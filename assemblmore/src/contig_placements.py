import numpy as np
import pandas as pd
from Bio import SeqIO
import networkx as nx

import os 


'''
TODO: 
In complex regions, primary contigs map to same regions.
Implement logic to handle these cases.

'''

def get_readlength_stats(reads_file):
    read_lengths = (len(record.seq) for record in SeqIO.parse(reads_file, "fasta"))
    return pd.Series(read_lengths).describe()



def filter_and_orient_contigs(mapped_contigs_path, min_contig_length=1000, phred_threshold=40):

    mappings_df = pd.read_csv(mapped_contigs_path, delimiter = '\t', header = None)
    mappings_df.columns = (
    ['contig', 'contig_length', 'contig_map_start', 'contig_map_end', 
    'strand', 'chr', 'chr_length', 'chr_map_start', 'chr_map_end', 'matched_bases', 'matched_total', 'mapping_quality'] + 
    [f"sup_{i}" for i in range(mappings_df.shape[1] - 12)]
    )


    organelle_options = ['MT', 'mitochondrion', 'chloroplast', 'plastid', 'organelle']
    organelle_contigs = mappings_df[mappings_df['chr'].isin(organelle_options)] # For now, ok, in future, filter out chloroplasts and use regex to account for multiple formats
    mappings_df = mappings_df[(~mappings_df['chr'].isin(organelle_options)) & (mappings_df['contig_length'] > min_contig_length) & (mappings_df['sup_0'] == 'tp:A:P')]


    #########################################
    prom_map = {item : {chr : 0 for chr in mappings_df['chr'].unique()} for item in mappings_df['contig'].unique()}
    for key, value in prom_map.items():
        for chr, count in value.items():
            prom_map[key][chr] += int(mappings_df[(mappings_df['contig'] == key) & (mappings_df['chr'] == chr)]['matched_total'].sum())


    contig_to_chr = pd.DataFrame(sorted([(item, next(key for key,value in prom_map[item].items() 
                                    if value == max(prom_map[item].values()))) for item in prom_map.keys()], 
                                    key = lambda x: x[1]), columns=['contig', 'chr'])
    ######################################### This really needs to be improved ^ 

 

    #####################################################
    # Calculate strand mode and max chr_map_end for each contig-chr pair and add as columns
    strand_modes, max_chr_map_ends, min_chr_map_starts, max_matched = [], [], [], [] #latest edit
    strand_map = {'+': 1, '-': -1}  # Map strand to numerical values for mode calculation
    for item in contig_to_chr.itertuples(index=False):

        poi_df = mappings_df[(mappings_df['contig'] == item[0]) & (mappings_df['chr'] == item[1])]

        ################################ Filter out extreme mappings
        most_contiguous = poi_df.loc[poi_df['matched_total'].idxmax()]
        most_contiguous_start = most_contiguous['chr_map_start']
        most_contiguous_end = most_contiguous['chr_map_end']
        matched_len = most_contiguous['contig_length']
        start_threshold = most_contiguous_start - matched_len
        end_threshold = most_contiguous_end + matched_len
        poi_df = poi_df[(poi_df['chr_map_start'] >= start_threshold) & (poi_df['chr_map_end'] <= end_threshold)]
        ################################


        ############# Calculate strand mode using a weighted approach

        matched_total = poi_df['matched_total'].sum()
        length_ratio_series = (poi_df['contig_map_end'] - poi_df['contig_map_start']) / matched_total
        criterion = (poi_df['strand'].apply(lambda x: strand_map[x] if x in strand_map else np.nan) @ length_ratio_series)
        if criterion >= 0:
            strand_modes.append('+')
        else:
            strand_modes.append('-')

        #mode_series = poi_df['strand'].mode() <- old method, but does not work well for complex regions with multiple mappings
        #############

        

        #max_end = mappings_df[(mappings_df['contig'] == item[0]) & (mappings_df['chr'] == item[1]) & (mappings_df['mapping_quality'] == 60)]['chr_map_end'].max()
        max_end = poi_df[poi_df['mapping_quality'] >= phred_threshold]['chr_map_end'].max() #latest edit
        #min_start = mappings_df[(mappings_df['contig'] == item[0]) & (mappings_df['chr'] == item[1]) & (mappings_df['mapping_quality'] == 60)]['chr_map_start'].min() #latest edit
        min_start = poi_df[poi_df['mapping_quality'] >= phred_threshold]['chr_map_start'].min() #latest edit
        #This will also have some bugs, check when more data available

        if pd.isnull(max_end):
            print(f"Warning: No valid chr_map_end found for contig {item[0]} on chr {item[1]} at phred threshold of {phred_threshold}.")

        max_chr_map_ends.append(max_end if pd.notnull(max_end) else None)
        min_chr_map_starts.append(min_start if pd.notnull(min_start) else None) #latest edit
        max_matched.append(poi_df['matched_total'].max() if not poi_df.empty else 0) #latest edit

    contig_to_chr['strand_mode'] = strand_modes
    contig_to_chr['min_chr_map_start'] = min_chr_map_starts  # latest edit
    contig_to_chr['max_chr_map_end'] = max_chr_map_ends
    contig_to_chr['max_matched'] = max_matched

    ####################################################

    G = nx.Graph()
    for i, item in contig_to_chr.iterrows():
        G.add_node(item['contig'], chr=item['chr'], strand_mode=item['strand_mode'], 
                   min_chr_map_start=item['min_chr_map_start'], max_chr_map_end=item['max_chr_map_end'],
                   max_matched=item['max_matched'])
    x = 0.5

    for i, item in contig_to_chr.iterrows():
        for j, other_item in contig_to_chr.iterrows():
            if i < j and item['chr'] == other_item['chr']:
                # Calculate overlap length
                overlap_start = max(item['min_chr_map_start'], other_item['min_chr_map_start'])
                overlap_end = min(item['max_chr_map_end'], other_item['max_chr_map_end'])
                overlap = max(0, overlap_end - overlap_start + 1)
                
                len1 = item['max_chr_map_end'] - item['min_chr_map_start'] + 1
                len2 = other_item['max_chr_map_end'] - other_item['min_chr_map_start'] + 1
                # Require overlap to be greater than fraction x of either contig's mapped region
                if (len1 > 0 and len2 > 0 and 
                    ((overlap / len1) > x or (overlap / len2) > x)):
                    # Add an edge if the overlap condition is met
                    G.add_edge(item['contig'], other_item['contig'])

    groups = list(nx.connected_components(G))

    contig_to_group = {label: i for i, group in enumerate(groups) for label in group}
    contig_to_chr['group'] = contig_to_chr['contig'].map(contig_to_group)
    ###################################################


    output_df = contig_to_chr.sort_values(by=['chr', 'min_chr_map_start'])

    # Optionally add organelle contigs if present
    if not organelle_contigs.empty:
        organelle_df = pd.DataFrame({
            'contig': organelle_contigs['contig'].unique(),
            'chr': organelle_contigs['chr'].unique(),
            'strand_mode': ['N/A'] * len(organelle_contigs['contig'].unique()),
            'max_chr_map_end': ['N/A'] * len(organelle_contigs['contig'].unique())
        })
        output_df = pd.concat([output_df, organelle_df], ignore_index=True)
    #^ This will fail if a contig is mapped to multiple organelles, but for now, this is ok
    return output_df


import sys

if len(sys.argv) < 3:
    print("Usage: python contig_placements.py <mapped_contigs_path> <reads_file>")
    sys.exit(1)

read_stats = get_readlength_stats(sys.argv[2])
threshold = read_stats['mean'] * 2
filtered_data = filter_and_orient_contigs(sys.argv[1], min_contig_length=threshold)


def removesuffix(s, suffix):
    if s.endswith(suffix):
        return s[:-len(suffix)]
    return s

rem_suffix = removesuffix(sys.argv[1], ".sorted.paf")

# Extract reference name after "to_"
if "to_" in rem_suffix:
    base_name_reference = rem_suffix.split("to_")[-1]
else:
    base_name_reference = rem_suffix  # fallback if "to_" not found

# Remove extensions if present
for ext in [".fasta", ".fa", ".fna"]:
    base_name_reference = removesuffix(base_name_reference, ext)


filtered_data.to_csv(f'filtered_by_{base_name_reference}_contigs.tsv', sep='\t', index=False)
print(f"Filtered contigs saved to 'filtered_by_{base_name_reference}_contigs.tsv'")
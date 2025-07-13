import pandas as pd 
import more_itertools as mit
from typing import List
from Bio import SeqIO
from Bio import SeqRecord
import sys
from copy import deepcopy
import os
import click

#import re

def get_readlength_stats(reads_file):
    print(f"[DEBUG] Reading read lengths from: {reads_file}")
    read_lengths = (len(record.seq) for record in SeqIO.parse(reads_file, "fasta"))
    stats = pd.Series(read_lengths).describe()
    print(f"[DEBUG] Read length stats:\n{stats}")
    return stats


def telomere_extension(alignments: pd.DataFrame, orderings: pd.DataFrame, expected_telomere_length: int = 8000, length_threshold: int = 0, phred_threshold: int = 20) -> List[List[str]]:
    """
    Extends the contigs at the telomeres using the best spanning reads.
    This function is a placeholder for future implementation.
    """
    print("[DEBUG] Running telomere_extension...")
    ###############################
    by_contig = [(i, item) for i, item in alignments.groupby(5, sort=False)]
    t2t_contigs = [item for item in orderings.iterrows() if orderings['chr'].tolist().count(item[1]['chr']) <= 1]
    both, left, right = [], [], []
    
    for item in by_contig:

        chr = orderings[orderings['contig'] == item[0]]['chr'].iloc[0]
        relevant_df = orderings[orderings['chr'] == chr]['contig'].to_list()
        num_contigs = len(relevant_df)

        #print(f"{item[0]} is on chromosome {chr} with {num_contigs} contigs in this chromosome.")
        if item[0] in [x[1]['contig'] for x in t2t_contigs]:

            if (orderings[orderings['contig'] == item[0]]['chr'] == 'MT').any(): #<- regex instead of hardcoding 'MT' would be better. 
                #print(f"Contig {item[0]} spans an organelle genome, skipping telomere extension.")
                pass
            else:
                #print(f"Contig {item[0]} is a telomere to telomere contig extend both ends.")
                both.append(item)

        else:
            item_idx = relevant_df.index(item[0])
            if  item_idx == 0:
                #print(f"Contig {item[0]} is not a t2t contig extend left end.")
                left.append(item)
            elif item_idx == num_contigs - 1:
                #print(f"Contig {item[0]} is not a t2t contig extend right end.")
                right.append(item)
            else:
                #print(f"Contig {item[0]} is not a t2t contig, but it is not at the end of the chromosome, skipping telomere extension.")
                continue
            
    ############################### Groups contigs into three categories: extend both ends, extend left end, extend right end, extend neither.
    print(f"[DEBUG] Found {len(both)} contigs to extend both ends, {len(left)} contigs to extend left end, and {len(right)} contigs to extend right end.")

    def extend(contig: pd.DataFrame, expected_telomere_length: int = 8000, phred_threshold: int = 20, length_threshold: int = 0, left: bool = True) -> pd.Series:
        """
        Extends the left end of the contig using the best spanning reads.
        """
        #print(f"[DEBUG] Extending {'left' if left else 'right'} end of contig...")
        # Placeholder for future implementation
        if left:
            #print("Extending left end of the contig.")
        
            points_of_interest = contig[(contig[1] >= length_threshold) 
                                & (contig[11] >= phred_threshold)

                                & (contig[4] == '+') # Not necessary but convenient. If needed take negative strand reads into account.
                                & (contig[7] < 100) #<- specific to left extension
                                & (contig[2] <= expected_telomere_length) # Make sure to tell users to overestimate the expected telomere length rather than underestimate it.
                                                                           # If they underestimate it, this will get rid of good reads.
                                ].copy()
            
            points_of_interest['criterion'] = points_of_interest[2]
            
        else:
            #print("Extending right end of the contig.")
            points_of_interest = contig[(contig[1] >= length_threshold) 
                                & (contig[11] >= phred_threshold)

                                & (contig[4] == '+') # Not necessary but convenient. If needed take negative strand reads into account.
                                & (contig[6] - contig[8] < 100) #<- specific to right extension (the specific value needs to be determined)
                                & (contig[1] - contig[3] <= expected_telomere_length) # Make sure to tell users to overestimate the expected telomere length rather than underestimate it.
                                                                        # If they underestimate it, this will get rid of good reads.
                                
                                ].copy()
            
            points_of_interest['criterion'] = points_of_interest[1] - points_of_interest[3]
        
            
        best_read = points_of_interest.loc[points_of_interest['criterion'].idxmax()] if not points_of_interest.empty else None
        
        return best_read
    
    # Extend both ends, left end, and right end of the contigs.
    #print("[DEBUG] Extending both ends, left end, and right end of the contigs...")
    both_t = [(item[0], extend(item[1], left = True, phred_threshold=phred_threshold, length_threshold = length_threshold), extend(item[1], left = False, phred_threshold=phred_threshold, length_threshold = length_threshold)) for item in both]
    left_t = [(item[0], extend(item[1], left = True, phred_threshold=phred_threshold, length_threshold = length_threshold)) for item in left]
    right_t = [(item[0], extend(item[1], left = False, phred_threshold=phred_threshold, length_threshold = length_threshold)) for item in right]
    #print(both_t)
    return [both_t, left_t, right_t]



failure_count = 0
def simple_contig_span(alignments: pd.DataFrame, orderings: pd.DataFrame, phred_threshold: int = 10, length_threshold: int = 0) -> List[str]:
    print("[DEBUG] Running simple_contig_span...") 
    def find_best_read(merged_alignments: pd.DataFrame, phred_threshold: int = 10, length_threshold: int = 0) -> pd.Series:
        
        #print(f"This is the input {merged_alignments}")

        if not isinstance(merged_alignments, pd.DataFrame):
            
            print(f"[DEBUG] Warning: Empty alignment found between {merged_alignments[0]} and {merged_alignments[1]}. Skipping.")
            #global failure_count
            #failure_count += 1
            return merged_alignments[1]
        
          # Return the right end of the alignment as a fallback
        # def get_clip_length(x: str) -> int:
        #     """
        #     Gets Hard/Soft clip length from the PAF file format.
        #     Returns the length of the clip if it exists, otherwise returns 0.
        #     """
        #     num = re.split('S|H', x)[0]
        #     return int(num) if num.isdigit() else 0
        
        print(f"[DEBUG] Finding best spanning read for {merged_alignments['5_x'].values[0]} and {merged_alignments['5_y'].values[0]}")
        points_of_interest = merged_alignments[(merged_alignments['1_x'] >= length_threshold)
                                               & (merged_alignments['6_x'] - merged_alignments['8_x'] < 100)
                                               & (merged_alignments['7_y'] < 100)
                                               #Above conditions are variable for potential branch and bound algorithm
                                               #Below conditions are fixed:
                                               & (merged_alignments['4_x'] == '+') 
                                               & (merged_alignments['4_y'] == '+')
                                               & (merged_alignments['11_x'] >= phred_threshold) #| (merged_alignments['11_y'] >= 60))
                                               & (merged_alignments['11_y'] >= phred_threshold)
                                                
                                               #& (merged_alignments['2_y'] > merged_alignments['3_x'])
                                               ].copy()

        #points_of_interest['criterion'] = points_of_interest['17_y'].apply(get_clip_length)
        points_of_interest['criterion'] = points_of_interest['2_y'] - points_of_interest['3_x']

        dup_chroms = [idx for idx, item in orderings['chr'].items() if orderings['chr'].tolist().count(item) > 1]
        chrom_windows = [list(item[1].itertuples()) for item in list(orderings.groupby('chr', sort=False))]



        #print(chrom_windows[i])
        def get_expected_gap_sizes(chr: List[pd.Series]) -> List[int]:
            res = []

            p1 = 0
            p2 = 1

            while p2 < len(chr):
                #print(f"Comparing {chr[p2][4]} and {chr[p1][5]}")
                res.append((chr[p1][1], chr[p2][1], chr[p2][4] - chr[p1][5]))
                p1 += 1
                p2 += 1
            return res
        
        expected_gap_sizes = [get_expected_gap_sizes(item) for item in chrom_windows]
        #for item in expected_gap_sizes:
            #print(f"[DEBUG] Expected gap sizes: {item}")

        #print(expected_gap_sizes)
        try:

            best_read = points_of_interest.loc[points_of_interest['criterion'].idxmax()]
        except:
            #global failure_count
            #failure_count += 1
            print(f"[DEBUG] Warning: No valid spanning read found between {merged_alignments['5_x'].values[0]} and {merged_alignments['5_y'].values[0]}. Failure count: {failure_count}")
            return merged_alignments['5_y'].values[0]

        # Return entire row as a Series (will need all data for merging of contigs later)
        return best_read
    

    by_contig = list(zip([item for _, item in alignments.groupby(5, sort=False)], orderings['chr']))

    # Group by contig and orderings to find pairs of alignments that are on the same chromosome
    print("[DEBUG] Merging alignments to find common reads that span the contigs...")
    merged = list(map(lambda x: pd.merge(x[0][0], x[1][0], how = 'inner', on = 0), (filter(lambda x: x[0][1] == x[1][1], mit.pairwise(by_contig)))))
    
    #print([(merged[i-1], item, merged[i+1]) for i, item in enumerate(merged) if item.empty])

    for i, item in enumerate(merged):
        if item.empty:
            print(f"[DEBUG] Warning: Empty alignment found between {merged[i-1]['5_y'].values[0]} and {merged[i+1]['5_x'].values[0]}. Skipping.")
            item = (merged[i-1]['5_y'].values[0], merged[i+1]['5_x'].values[0])
            continue

    merged = [item if not item.empty else (merged[i-1]['5_y'].values[0], merged[i+1]['5_x'].values[0]) for i, item in enumerate(merged)]
    
    # merge alignments to find common reads that span the contigs (probably no point in import more itertools just for the pairwise function, but it is more readable this way)
    #spanning_reads = [find_best_read(item) for item in merged if not item.empty]
    
    print("[DEBUG] Finding best spanning reads for merged alignments...")
    spanning_reads = [find_best_read(item, phred_threshold=phred_threshold, length_threshold=length_threshold) for item in merged]

    return spanning_reads

@click.command()
@click.argument('alignments_file', type=click.Path(exists=True))
@click.argument('orderings_file', type=click.Path(exists=True))
@click.argument('contigs_file', type=click.Path(exists=True))
@click.argument('reads_file', type=click.Path(exists=True))
@click.option('--expected_telomere_length', default=8000, help='Expected length of telomeres (default: 8000)')
@click.option('--length_threshold', default=0, help='Minimum length of reads to consider (default: 0)')
@click.option('--phred_threshold', default=20, help='Minimum Phred quality score of reads to consider (default: equivalent to 1/100 chance of mapping misplacement <20 for telomere extension, 10 for spanning reads>)')
def merge_contigs(alignments_file, orderings_file, contigs_file, reads_file, expected_telomere_length, length_threshold, phred_threshold) -> dict:


    click.echo(f"[DEBUG] Starting merge_contigs with parameters:\n")
    click.echo(f"alignments_file: {alignments_file}")
    click.echo(f"orderings_file: {orderings_file}")
    click.echo(f"contigs_file: {contigs_file}")
    click.echo(f"reads_file: {reads_file}")
    click.echo(f"expected_telomere_length: {expected_telomere_length}")
    click.echo(f"length_threshold: {length_threshold}")
    click.echo(f"phred_threshold: {phred_threshold}")

    print("[DEBUG] Reading input files...")
    alignments = pd.read_csv(alignments_file, delimiter='\t', header=None)
    orderings = pd.read_csv(orderings_file, delimiter = '\t')
    contigs = list(SeqIO.parse(contigs_file, "fasta"))
    reads = list(SeqIO.parse(reads_file, "fasta"))

    print("[DEBUG] Running merge_contigs...")
    print(orderings.head())
    groupings = list(orderings.groupby('chr', sort=False))
    chr_to_contigs = {chr: contigs['contig'].to_list() for chr, contigs in groupings}
    contigs_to_chr = {contig: chr for chr, contigs in groupings for contig in contigs['contig'].to_list()}

    print("[DEBUG] Running telomere_extension in merge_contigs...")
    both, left, right = telomere_extension(alignments, orderings, 
                                           expected_telomere_length=expected_telomere_length, 
                                           length_threshold=length_threshold, 
                                           phred_threshold=phred_threshold)
    
    print("[DEBUG] Running simple_contig_span in merge_contigs...")
    spanning_reads = simple_contig_span(alignments, orderings, 
                                        length_threshold=length_threshold, 
                                        phred_threshold=phred_threshold)  #<- phred_threshold is halved for spanning reads since the mapping is paired.
    
    final_hash = deepcopy(chr_to_contigs)
    print("[DEBUG] Inserting left, right, and both telomere extensions into final_hash...")
    for item in left:
        chr = contigs_to_chr[item[0]]
        final_hash[chr].insert(0, item[1]) if item[1] is not None else final_hash[chr].insert(0, None)

    for item in right:
        chr = contigs_to_chr[item[0]]
        final_hash[chr].append(item[1]) if item[1] is not None else final_hash[chr].append(None)

    for item in both:
        chr = contigs_to_chr[item[0]]
        final_hash[chr].insert(0, item[1]) if item[1] is not None else final_hash[chr].insert(0, None)
        final_hash[chr].append(item[2]) if item[2] is not None else final_hash[chr].append(None)

    def get_index(item, query):
        """
        Custom .index() pandas series break default
        """
        for i, contig in enumerate(item):
            if type(contig) != str: 
                pass
            elif contig == query:
                return i
        return(None)
                
              
    print("[DEBUG] Inserting spanning reads into final_hash...")
    for i, item in enumerate(spanning_reads):
        chr = contigs_to_chr[item['5_x']] if type(item) != str else contigs_to_chr[item]
        
        if type(item) != str:
            index = get_index(final_hash[chr], item['5_y'])
            final_hash[chr].insert(index, item)
        else:
            index = get_index(final_hash[chr], item)
            final_hash[chr].insert(index, None)
    
    read_ids = [item.id for item in reads]
    read_seqs = [item.seq for item in reads]
    contig_ids = [item.id for item in contigs]
    contig_seqs = [item.seq for item in contigs]

    #test = deepcopy(final_hash['III'])
    def combine_contigs(arr: list) -> list:
        print("[DEBUG] Combining contigs for a chromosome...")


        ################################## This has a logical error, if compound merges this breaks.
        def inter_contig_merge(str1, str2, series_rules):
            print("[DEBUG] Running inter_contig_merge...")
            if not isinstance(str1, str):
                l_str = str1
                base_x = len(l_str) - series_rules['7_x'] + series_rules['8_x'] #If str1 is SeqIO object, that means it has been merged to some extent and the indices stored in '8_x' are incorrect.
                
            else:
                l_str = contig_seqs[contig_ids.index(str1)]
                base_x = series_rules['8_x']
            #^^^^ If inter merge is used, str1 is a Bio.Seq.Seq object
            if not isinstance(str2, str):
                r_str = str2
            else:
                r_str = contig_seqs[contig_ids.index(str2)]
            #^^^^ If inter merge is used, str2 is a Bio.Seq.Seq object


            base_y = series_rules['7_y']
            #print(f"[DEBUG] base_x: {base_x}, base_y: {base_y}, str1: {str1}, str2: {str2}")

            x_prime = series_rules['3_x']
            y_prime = series_rules['2_y']
            

            if series_rules['criterion'] < 0:
                x = base_x - (series_rules['3_y'] - series_rules['2_x']) #<- this is the length of the left contig, not the right one.
                y = base_y + (series_rules['3_y'] - series_rules['2_x']) #<- this is the length of the right contig, not the left one.
                x_prime, y_prime = y_prime, x_prime 
            else:
                x = base_x
                y = base_y

            mid_str = read_seqs[read_ids.index(series_rules[0])]
            

            #print(f"[DEBUG] {series_rules[0]}_span_len {len(mid_str[x_prime + 1:y_prime])}, {str1}_len_{len(l_str[:x])}, {str2}_len_{len(r_str[y:])}")
            

            return(l_str[:x] + mid_str[x_prime + 1:y_prime] + r_str[y:])      #<----not sure +1 or not, but it is not important now. Also read mismatches at ends not accounted for yet.

        ##################################

        def telomere_merge(str1, series_rules, left = True):
            """
            Merges the contig with the telomere using the best spanning read.
            """
            print(f"[DEBUG] Running telomere_merge ({'left' if left else 'right'})...")
            if not isinstance(str1, str):
                rel_str = str1
            else:
                rel_str = contig_seqs[contig_ids.index(str1)]
            #^^^^ If inter merge is used, str1 is a Bio.Seq.Seq object, otherwise it is a string.

            telo_str = read_seqs[read_ids.index(series_rules[0])]

            if left:
                x = series_rules[2]
                #print(f"[DEBUG] {series_rules[0]}_span_len {len(telo_str[:x])}, {str1}_len_{len(rel_str)}")
                return telo_str[:x] + rel_str
            else:
                y = series_rules[3]
                #print(f"[DEBUG] {series_rules[0]}_span_len {len(telo_str[y+1:])}, {str1}_len_{len(rel_str)}")
                return rel_str + telo_str[y + 1:]

        
        i = 1
        while i < len(arr) - 2:
            #print(f"[DEBUG] Checking for inter-contig merge at position {i}...")
            s1 = arr[i]
            s2 = arr[i + 2]
            merge_rule = arr[i + 1]

            if isinstance(merge_rule, pd.Series):
                #print(f"[DEBUG] Performing inter-contig merge between {s1} and {s2}...")
                merged = inter_contig_merge(s1, s2, merge_rule)
                # Replace s1 and s2 with merged string
                arr[i] = merged
                # Remove merge rule and s2
                del arr[i + 1:i + 3]
                # Step back to try more merges with new merged result
                i = max(i - 2, 1)
            else:
                i += 2
        

        
        # Handle telomere merges
        if isinstance(arr[0], pd.Series):
            #print("[DEBUG] Performing left telomere merge...")
            arr[0] = telomere_merge(arr[1], arr[0], left=True)
            del arr[1]
        if isinstance(arr[-1], pd.Series):
            #print("[DEBUG] Performing right telomere merge...")
            arr[-1] = telomere_merge(arr[-2], arr[-1], left=False)
            del arr[-2]
        
        
        return(arr)
    


    def generate_descriptions(contigs: List[str]) -> List[str]:
        if isinstance(contigs, pd.Series):
            return(contigs[0])
        elif contigs is None:
            return('None')
        else:
            contig_length = len(contig_seqs[contig_ids.index(contigs)]) if contigs else 0
            return(f"{contigs}_len_{contig_length}")


    print("[DEBUG] Combining contigs for all chromosomes...")
    final_merged = {chr: combine_contigs(item) for chr, item in deepcopy(final_hash).items()}
    print("[DEBUG] Finished merge_contigs.")

    print("[DEBUG] Generating descriptions for contigs...")
    final_labels = {chr: ('-'.join([generate_descriptions(item) for item in final_hash[chr]])).split('-None-') for chr in final_hash}
    print(f"[DEBUG] Generated descriptions: {final_labels}")



    for chr, contigs in final_merged.items():
        final_merged[chr] = [item for item in final_merged[chr] if item is not None]  # Remove None values so the indices match the labels


    print("[DEBUG] Preparing output sequences...")
    fasta_record = []
    for chr, contigs in final_merged.items():
        for i, contig in enumerate(contigs):
            if isinstance(contig, str):
                #print(f"  Contig: {contig} (length: {len(contig)})")
                contig = contig_seqs[contig_ids.index(contig)]
                contig_record = SeqRecord.SeqRecord(seq=contig, id=f"{chr}_contig_{i+1}", description = final_labels[chr][i] if i < len(final_labels[chr]) else "")
                fasta_record.append(contig_record)
                i += 1
            elif contig is None:
                pass
            else:
                #print(f"  Contig: {contig} (length: {len(contig)})")
                contig_record = SeqRecord.SeqRecord(seq=contig, id=f"{chr}_contig_{i+1}", description = final_labels[chr][i] if i < len(final_labels[chr]) else "")
                fasta_record.append(contig_record)
                i += 1


    print("[DEBUG] Final merged result:")
    print(final_merged)

    cur_dir = os.getcwd()
    print(f"[DEBUG] Outputting FASTA file to: {cur_dir}")
    SeqIO.write(fasta_record, os.path.join(cur_dir, "final_assembly.fasta"), "fasta")
    print("[DEBUG] Finished writing output files.")

    return final_merged, fasta_record


if __name__ == "__main__":
    merge_contigs()
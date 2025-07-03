import pandas as pd 
import more_itertools as mit
import re
from typing import List
from Bio import SeqIO
from Bio import SeqRecord
import sys
from copy import deepcopy


def get_readlength_stats(reads_file):
    print(f"Reading read lengths from: {reads_file}")
    read_lengths = (len(record.seq) for record in SeqIO.parse(reads_file, "fasta"))
    stats = pd.Series(read_lengths).describe()
    print(f"Read length stats:\n{stats}")
    return stats


def telomere_extension(alignments: pd.DataFrame, orderings: pd.DataFrame, expected_telomere_length: int = 8000) -> List[List[str]]:
    """
    Extends the contigs at the telomeres using the best spanning reads.
    This function is a placeholder for future implementation.
    """
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
    print(f"Found {len(both)} contigs to extend both ends, {len(left)} contigs to extend left end, and {len(right)} contigs to extend right end.")

    def extend(contig: pd.DataFrame, expected_telomere_length: int = 8000, left: bool = True) -> pd.Series:
        """
        Extends the left end of the contig using the best spanning reads.
        """
        # Placeholder for future implementation
        if left:
            #print("Extending left end of the contig.")
        
            points_of_interest = contig[(contig[1] >= 100000) 
                                & (contig[11] >= 60)

                                & (contig[4] == '+') # Not necessary but convenient. If needed take negative strand reads into account.
                                & (contig[7] < 100) #<- specific to left extension
                                & (contig[2] <= expected_telomere_length) # Make sure to tell users to overestimate the expected telomere length rather than underestimate it.
                                                                           # If they underestimate it, this will get rid of good reads.
                                ].copy()
            
            points_of_interest['criterion'] = points_of_interest[2]
            
        else:
            #print("Extending right end of the contig.")
            points_of_interest = contig[(contig[1] >= 100000) 
                                & (contig[11] >= 60)

                                & (contig[4] == '+') # Not necessary but convenient. If needed take negative strand reads into account.
                                & (contig[6] - contig[8] < 100) #<- specific to right extension (the specific value needs to be determined)
                                & (contig[1] - contig[3] <= expected_telomere_length) # Make sure to tell users to overestimate the expected telomere length rather than underestimate it.
                                                                        # If they underestimate it, this will get rid of good reads.
                                
                                ].copy()
            
            points_of_interest['criterion'] = points_of_interest[1] - points_of_interest[3]
        
            
        best_read = points_of_interest.loc[points_of_interest['criterion'].idxmax()] if not points_of_interest.empty else None
        
        return best_read
    
    # Extend both ends, left end, and right end of the contigs.
    both_t = [(item[0], extend(item[1], left = True), extend(item[1], left = False)) for item in both]
    left_t = [(item[0], extend(item[1], left = True)) for item in left]
    right_t = [(item[0], extend(item[1], left = False)) for item in right]
    #print(both_t)
    return [both_t, left_t, right_t]



failure_count = 0
def simple_contig_span(alignments: pd.DataFrame, orderings: pd.DataFrame) -> List[str]:
     
    def find_best_read(merged_alignments: pd.DataFrame) -> pd.Series:
        
        #print(f"This is the input {merged_alignments}")

        if not isinstance(merged_alignments, pd.DataFrame):
            
            print(f"Warning: Empty alignment found between {merged_alignments[0]} and {merged_alignments[1]}. Skipping.")
            #global failure_count
            #failure_count += 1
            return merged_alignments[1]
        
          # Return the right end of the alignment as a fallback
        def get_clip_length(x: str) -> int:
            """
            Gets Hard/Soft clip length from the PAF file format.
            Returns the length of the clip if it exists, otherwise returns 0.
            """
            num = re.split('S|H', x)[0]
            return int(num) if num.isdigit() else 0
        
        print(f"Finding best spanning read for {merged_alignments['5_x'].values[0]} and {merged_alignments['5_y'].values[0]}")
        points_of_interest = merged_alignments[(merged_alignments['1_x'] >= 100000)
                                               & (merged_alignments['6_x'] - merged_alignments['8_x'] < 10)
                                               & (merged_alignments['7_y'] < 10)
                                               #Above conditions are variable for potential branch and bound algorithm
                                               #Below conditions are fixed:
                                               & (merged_alignments['4_x'] == '+') 
                                               & (merged_alignments['4_y'] == '+')
                                               & ((merged_alignments['11_x'] >= 60) | (merged_alignments['11_y'] >= 60))
                                               #& (merged_alignments['11_y'] >= 60)
                                                
                                               & (merged_alignments['2_y'] > merged_alignments['3_x'])
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
        for item in expected_gap_sizes:
            print(item)

        #print(expected_gap_sizes)
        try:

            best_read = points_of_interest.loc[points_of_interest['criterion'].idxmax()]
        except:
            #global failure_count
            #failure_count += 1
            print(f"Warning: No valid spanning read found between {merged_alignments['5_x'].values[0]} and {merged_alignments['5_y'].values[0]}. Failure count: {failure_count}")
            return merged_alignments['5_y'].values[0]

        # Return entire row as a Series (will need all data for merging of contigs later)
        return best_read
    

    by_contig = list(zip([item for _, item in alignments.groupby(5, sort=False)], orderings['chr']))

    # Group by contig and orderings to find pairs of alignments that are on the same chromosome
    merged = list(map(lambda x: pd.merge(x[0][0], x[1][0], how = 'inner', on = 0), (filter(lambda x: x[0][1] == x[1][1], mit.pairwise(by_contig)))))
    
    #print([(merged[i-1], item, merged[i+1]) for i, item in enumerate(merged) if item.empty])

    for i, item in enumerate(merged):
        if item.empty:
            print(f"Warning: Empty alignment found between {merged[i-1]['5_y'].values[0]} and {merged[i+1]['5_x'].values[0]}. Skipping.")
            item = (merged[i-1]['5_y'].values[0], merged[i+1]['5_x'].values[0])
            continue

    merged = [item if not item.empty else (merged[i-1]['5_y'].values[0], merged[i+1]['5_x'].values[0]) for i, item in enumerate(merged)]
    
    # merge alignments to find common reads that span the contigs (probably no point in import more itertools just for the pairwise function, but it is more readable this way)
    #spanning_reads = [find_best_read(item) for item in merged if not item.empty]
    
    spanning_reads = [find_best_read(item) for item in merged]

    return spanning_reads
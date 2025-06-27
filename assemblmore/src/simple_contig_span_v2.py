import pandas as pd 
import more_itertools as mit
import re
from typing import List
from Bio import SeqIO
from Bio import SeqRecord
import sys
"""
TODO:
There is a chance that selected spanning reads are chimeric. Need to implement a check for this.
Additionally, there is a chance that the spanning reads still contain adapter sequences. Need to implement a check for this as well.

The script also needs to be able to span the gap between contigs where multiple reads are needed to connect the contigs.
---Base Case: the contigs can be spanned by a single read.
---Recursive Case: extend the left contig as much as possible and check for base case again.

---- Procedure Needs: python subprocess module that aligns reads to extended contig via minimap2.


The base case is that there is a single read that spans the gap between two contigs.

"""

#I'm pretty sure its a moot point since the dna extraction will sample DNA from multiple cells which will have variable lengths of telomeres.
#Update: This is not a moot point, Ryan mentioned that every replication cycle shortens the telomeres by around 2~6 bp (size of primase); this is results in a negliglbe difference between cells.

#NOTE: This function produces good results for the C. briggsae AF16 assembly, however, the section where a bound is set on edge mismatches needs to be variably adjusted. 


def telomere_extension(alignments: pd.DataFrame, orderings: pd.DataFrame, expected_telomere_length: int = 8000) -> List[List[str]]:
    """
    Extends the contigs at the telomeres using the best spanning reads.
    This function is a placeholder for future implementation.
    """
    ###############################
    by_contig = [(i, item) for i, item in alignments.groupby(5, sort=False)]
    t2t_contigs = [item for item in orderings.iterrows() if orderings['chr'].tolist().count(item[1]['chr']) <= 1]
    i = 0
    both, left, right = [], [], []
    for item in by_contig:
        if item[0] in [x[1]['contig'] for x in t2t_contigs]:

            if (orderings[orderings['contig'] == item[0]]['chr'] == 'MT').any(): #<- regex instead of hardcoding 'MT' would be better. 
                 print(f"Contig {item[0]} spans an organelle genome, skipping telomere extension.")
            else:
                print(f"Contig {item[0]} is a telomere to telomere contig extend both ends.")
                both.append(item)
        else:
            if i == 0:
                print(f"Contig {item[0]} is not a t2t contig extend left end.")
                left.append(item)
            elif i == 1:
                print(f"Contig {item[0]} is not a t2t contig extend right end.")
                right.append(item)
            i += 1
            i %= 2
    ############################### Disgusting implementation. Groups contigs into three categories: extend both ends, extend left end, and extend right end.
    
    def extend(contig: pd.DataFrame, expected_telomere_length: int = 8000, left: bool = True) -> pd.Series:
        """
        Extends the left end of the contig using the best spanning reads.
        """
        # Placeholder for future implementation
        if left:
            print("Extending left end of the contig.")
        
            points_of_interest = contig[(contig[1] >= 100000) 
                                & (contig[11] >= 60)

                                & (contig[4] == '+') # Not necessary but convenient. If needed take negative strand reads into account.
                                & (contig[7] < 100) #<- specific to left extension
                                & (contig[2] <= expected_telomere_length) # Make sure to tell users to overestimate the expected telomere length rather than underestimate it.
                                                                           # If they underestimate it, this will get rid of good reads.
                                ].copy()
            
        else:
            print("Extending right end of the contig.")
            points_of_interest = contig[(contig[1] >= 100000) 
                                & (contig[11] >= 60)

                                & (contig[4] == '+') # Not necessary but convenient. If needed take negative strand reads into account.
                                & (contig[6] - contig[8] < 100) #<- specific to right extension (the specific value needs to be determined)
                                & (contig[1] - contig[3] <= expected_telomere_length) # Make sure to tell users to overestimate the expected telomere length rather than underestimate it.
                                                                        # If they underestimate it, this will get rid of good reads.
                                
                                ].copy()
        
            
        if left:
            points_of_interest['criterion'] = points_of_interest[2]
        else:
            points_of_interest['criterion'] = points_of_interest[1] - points_of_interest[3]

        
        best_read = points_of_interest.loc[points_of_interest['criterion'].idxmax()]
        
        # print(best_read)
        return best_read
    
    print(len(both), len(left), len(right))
    
    

    both_t = [(item[0], extend(item[1], left = True), extend(item[1], left = False)) for item in both]
    left_t = [extend(item[1], left = True) for item in left]
    right_t = [extend(item[1], left = False) for item in right]
    #print(both_t)
    return [both_t, left_t, right_t]


def simple_contig_span(alignments: pd.DataFrame, orderings: pd.DataFrame) -> List[str]:
     
    def find_best_read(merged_alignments: pd.DataFrame) -> pd.Series:

        def get_clip_length(x: str) -> int:
            """
            Gets Hard/Soft clip length from the PAF file format.
            Returns the length of the clip if it exists, otherwise returns 0.
            """
            num = re.split('S|H', x)[0]
            return int(num) if num.isdigit() else 0
        
        points_of_interest = merged_alignments[(merged_alignments['1_x'] >= 100000)
                                               & (merged_alignments['6_x'] - merged_alignments['8_x'] < 10)
                                               & (merged_alignments['7_y'] < 10)
                                               #Above conditions are variable for potential branch and bound algorithm
                                               #Below conditions are fixed:
                                               & (merged_alignments['4_x'] == '+') 
                                               & (merged_alignments['4_y'] == '+')
                                               & (merged_alignments['11_x'] >= 60)
                                               & (merged_alignments['11_y'] >= 60)
                                               & (merged_alignments['2_y'] > merged_alignments['3_x'])
                                               ].copy()

        points_of_interest['criterion'] = points_of_interest['17_y'].apply(get_clip_length)

        best_read = points_of_interest.loc[points_of_interest['criterion'].idxmax()]

        # Return entire row as a Series (will need all data for merging of contigs later)
        return best_read
    

    by_contig = list(zip([item for _, item in alignments.groupby(5, sort=False)], orderings['chr']))

    # Group by contig and orderings to find pairs of alignments that are on the same chromosome
    merged = list(map(lambda x: pd.merge(x[0][0], x[1][0], how = 'inner', on = 0), (filter(lambda x: x[0][1] == x[1][1], mit.pairwise(by_contig)))))
    

    # merge alignments to find common reads that span the contigs (probably no point in import more itertools just for the pairwise function, but it is more readable this way)
    spanning_reads = [find_best_read(item) for item in merged if not item.empty]

    return spanning_reads

def remap_reads(relevant_contigs: List[SeqRecord.SeqRecord], all_reads: SeqIO.FastaIO.FastaIterator) -> pd.DataFrame:
    """
    Remaps reads to relevant contigs.
    This function is a placeholder for future implementation.
    """
    # Placeholder for future implementation
    return []


def contig_span(alignments: pd.DataFrame, orderings: pd.DataFrame) -> List[str]:
    """
    Spans the gap between contigs using the best spanning reads.
    This function is a placeholder for future implementation.
    """
    # Placeholder for future implementation
    
    return []


def connect_contigs(spanning_reads: List[pd.Series], contigs: SeqIO.FastaIO.FastaIterator, all_reads: SeqIO.FastaIO.FastaIterator, orderings: pd.DataFrame) -> List[str]:
    """
    Connects contigs based on spanning reads.
    This function is a placeholder for future implementation.
    """
    contigs = list(contigs)
    contig_ids = [item.id for item in contigs]
    contig_seqs = [item.seq for item in contigs]
    all_reads = list(all_reads)

    connected_contigs = []
    for df in spanning_reads:
        left_contig = contig_seqs[contig_ids.index(df['5_x'])]
        right_contig = contig_seqs[contig_ids.index(df['5_y'])]
        l_index = df['3_x']
        r_index = df['2_y']
        l_contig_l_index, r_contig_l_index = df['7_x'], df['7_y']
        l_contig_r_index, r_contig_r_index = df['8_x'], df['8_y']
        spanning_read = [item for item in all_reads if item.id == df[0]][0]
        chr = orderings[orderings['contig'] == df['5_x']]['chr'].values[0]
        connected_contigs.append(((left_contig + spanning_read.seq[l_index:r_index] + right_contig), chr, f"{chr}-{df['5_x']}-{df[0]}-{df['5_y']}-connected_at-{l_contig_r_index}"))
        # Does not handle the case where multiple reads are needed to connect the contigs.
        # Or in the case of disagreement between the overlapping part of the reads and the contigs.
    
    return (connected_contigs)





if __name__ == "__main__":


    if len(sys.argv) < 5:
        print("Usage: python simple_contig_span_v2.py <alignments_file> <orderings_file> <contigs_file> <reads_file>")
        sys.exit(1)

    #data = pd.read_csv(r"C:\Users\tkoti\Desktop\Genomes\gap_fill_automation\C_briggsae_AF16.fasta_mapped_to_ordered_and_oriented_assembly.sorted.sed.paf", delimiter='\t', header=None)
    data = pd.read_csv(sys.argv[1], delimiter='\t', header=None)
    ordering = pd.read_csv(sys.argv[2], delimiter='\t')
    #ordering = pd.read_csv(r"C:\Users\tkoti\Desktop\Genomes\gap_fill_automation\filtered_contigs.tsv", delimiter='\t')
    telomere_reads = telomere_extension(data, ordering, expected_telomere_length=8000)
    spanning_reads = simple_contig_span(data, ordering)


    # print(telomere_reads)
    # print(spanning_reads)
    #global_hash = {item : [] for item in ordering['chr'].unique()}



    
    #contigs = list(SeqIO.parse(r"C:\Users\tkoti\Desktop\Genomes\gap_fill_automation\ordered_and_oriented_assembly.fasta", 'fasta'))
    contigs = list(SeqIO.parse(sys.argv[3], 'fasta'))
    #all_reads = list(SeqIO.parse(r"C:\Users\tkoti\Desktop\Genomes\raw_reads\C_briggsae_AF16.fasta", 'fasta'))
    all_reads = list(SeqIO.parse(sys.argv[4], 'fasta'))


    #meh = (connect_contigs(spanning_reads, contigs, all_reads, ordering))
    dup_chroms = [item for item in ordering['chr'] if ordering['chr'].tolist().count(item) > 1]
    t2t_contigs = [item for item in ordering.iterrows() if item[1]['chr'] not in dup_chroms and not (ordering[ordering['contig'] == item[1]['contig']]['chr'] == 'MT').any()]

    filled_assembly = []
    print(len(telomere_reads[0][0]))
    print(telomere_reads[0][0][2])
    for i, item in enumerate(t2t_contigs):
        contig_id = item[1]['contig']
        contig_seq = [item for item in contigs if item.id == contig_id][0]
        chr = item[1]['chr']
        l_telomere = telomere_reads[0][i][1]
        l_telomere_seq = [item for item in all_reads if item.id == l_telomere[0]][0]
        r_telomere = telomere_reads[0][i][2]
        r_telomere_seq = [item for item in all_reads if item.id == r_telomere[0]][0]
        print('length of contig before filling:', len(contig_seq.seq))
        contig_seq.seq = l_telomere_seq.seq[:l_telomere[2]] + contig_seq.seq + r_telomere_seq.seq[r_telomere[3]:]
        print('length of contig after filling:', len(contig_seq.seq))
        filled_assembly.append(SeqRecord.SeqRecord(seq=contig_seq.seq, id=f"chr{chr}", description=f"{chr}-{contig_id}"))

    # for item in meh:
    #     filled_assembly.append(SeqRecord.SeqRecord(seq=item[0], id=f"chr{item[1]}", description=item[2]))





    #SeqIO.write(filled_assembly, 'automated_filled_assembly.fasta', 'fasta')
    # print(spanning_reads)
    
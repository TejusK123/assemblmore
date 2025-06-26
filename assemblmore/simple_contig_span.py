import pandas as pd 
import more_itertools as mit
import re
from typing import List
from Bio import SeqIO
from Bio import SeqRecord

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
                                               & (merged_alignments['4_x'] == '+') 
                                               & (merged_alignments['4_y'] == '+')
                                               & (merged_alignments['11_x'] >= 60)
                                               & (merged_alignments['11_y'] >= 60)
                                               ].copy()

        points_of_interest['clipping'] = points_of_interest['17_y'].apply(get_clip_length)

        best_read = points_of_interest.loc[points_of_interest['clipping'].idxmax()]

        # Return entire row as a Series (will need all data for merging of contigs later)
        return best_read
    

    by_contig = list(zip([item for _, item in alignments.groupby(5, sort=False)], orderings['chr']))

    # Group by contig and orderings to find pairs of alignments that are on the same chromosome
    merged = list(map(lambda x: pd.merge(x[0][0], x[1][0], how = 'inner', on = 0), (filter(lambda x: x[0][1] == x[1][1], mit.pairwise(by_contig)))))
    
    # merge alignments to find common reads that span the contigs (probably no point in import more itertools just for the pairwise function, but it is more readable this way)
    spanning_reads = [find_best_read(item) for item in merged if not item.empty]

    return spanning_reads


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

def telomere_extension(alignments: pd.DataFrame, orderings: pd.DataFrame) -> List[str]:
    """
    Extends telomeres based on alignments.
    This function is a placeholder for future implementation.
    """
    # Placeholder for future implementation
    return []


def remap_reads(relevant_contigs: List[SeqRecord.SeqRecord], all_reads: SeqIO.FastaIO.FastaIterator) -> pd.DataFrame:
    """
    Remaps reads to relevant contigs.
    This function is a placeholder for future implementation.
    """
    # Placeholder for future implementation
    return []

if __name__ == "__main__":
    data = pd.read_csv('C_briggsae_AF16.fasta_mapped_to_ordered_and_oriented_assembly.sorted.sed.paf', delimiter='\t', header=None)
    ordering = pd.read_csv('filtered_contigs.tsv', delimiter='\t')
    spanning_reads = simple_contig_span(data, ordering)


    contigs = list(SeqIO.parse('ordered_and_oriented_assembly.fasta', 'fasta'))
    all_reads = list(SeqIO.parse(r"C:\Users\tkoti\Desktop\Genomes\raw_reads\C_briggsae_AF16.fasta", 'fasta'))


    meh = (connect_contigs(spanning_reads, contigs, all_reads, ordering))
    dup_chroms = [item for item in ordering['chr'] if ordering['chr'].tolist().count(item) > 1]
    t2t_contigs = [item for item in ordering.iterrows() if item[1]['chr'] not in dup_chroms]
    
    filled_assembly = []

    for item in t2t_contigs:
        contig_id = item[1]['contig']
        contig_seq = [item for item in contigs if item.id == contig_id][0]
        chr = item[1]['chr']
        filled_assembly.append(SeqRecord.SeqRecord(seq=contig_seq.seq, id=f"chr{chr}", description=f"{chr}-{contig_id}"))

    for item in meh:
        filled_assembly.append(SeqRecord.SeqRecord(seq=item[0], id=f"chr{item[1]}", description=item[2]))

    #SeqIO.write(filled_assembly, 'automated_filled_assembly.fasta', 'fasta')
    # print(spanning_reads)
    
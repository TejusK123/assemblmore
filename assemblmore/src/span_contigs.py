import pandas as pd 
import more_itertools as mit
from typing import List
from Bio import SeqIO
from Bio import SeqRecord
import sys
from copy import deepcopy
import os
import click
import os
import subprocess
#import re



def combine_contigs(arr: list, reads, contigs, extended_contigs: dict = None) -> list:
    print("[DEBUG] Combining contigs for a chromosome...")
    
    if extended_contigs is None:
        extended_contigs = {}

    read_ids = [item.id for item in reads]
    read_seqs = [item.seq for item in reads]
    contig_ids = [item.id for item in contigs]
    contig_seqs = [item.seq for item in contigs]

    ################################## This has a logical error, if compound merges this breaks.
    def inter_contig_merge(str1, str2, series_rules):
        print("[DEBUG] Running inter_contig_merge...")
        if not isinstance(str1, str):
            l_str = str1
            base_x = len(l_str) - series_rules['7_x'] + series_rules['8_x'] #If str1 is SeqIO object, that means it has been merged to some extent and the indices stored in '8_x' are incorrect.
            
        else:
            print(f"[DEBUG] str1 is a string, using contig_seqs for {str1}...")
            # Check if this contig was extended, use extended sequence if available
            if str1 in extended_contigs:
                l_str = extended_contigs[str1]
                print(f"[DEBUG] Using extended sequence for {str1} (length: {len(l_str)})")
            else:
                l_str = contig_seqs[contig_ids.index(str1)]
            base_x = series_rules['8_x']
        #^^^^ If inter merge is used, str1 is a Bio.Seq.Seq object
        if not isinstance(str2, str):
            r_str = str2
        else:
            print(f"[DEBUG] str2 is a string, using contig_seqs for {str2}...")
            # Check if this contig was extended, use extended sequence if available  
            if str2 in extended_contigs:
                r_str = extended_contigs[str2]
                print(f"[DEBUG] Using extended sequence for {str2} (length: {len(r_str)})")
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
            # Check if this contig was extended, use extended sequence if available
            if str1 in extended_contigs:
                rel_str = extended_contigs[str1]  
                print(f"[DEBUG] Using extended sequence for telomere merge {str1} (length: {len(rel_str)})")
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
def simple_contig_span(alignments: pd.DataFrame, orderings: pd.DataFrame, reads, contigs, reads_file, phred_threshold: int = 10, length_threshold: int = 0, enable_extension: bool = True) -> tuple:
    print("[DEBUG] Running simple_contig_span...") 
    read_ids = [item.id for item in reads]
    read_seqs = [item.seq for item in reads]
    contig_ids = [item.id for item in contigs]
    contig_seqs = [item.seq for item in contigs]
    
    # Track extended contigs - maps original contig name to extended sequence
    extended_contigs = {}


    def find_best_read(merged_alignments: pd.DataFrame, whole_alignments: pd.DataFrame, phred_threshold: int = 10, length_threshold: int = 0, max_iterations: int = 3, enable_extension: bool = True) -> pd.Series:  
        """
        Find the best spanning read with iterative contig extension if needed.
        
        Args:
            merged_alignments: DataFrame with merged alignment data
            whole_alignments: Full alignments DataFrame for extension search
            phred_threshold: Minimum quality threshold
            length_threshold: Minimum read length threshold  
            max_iterations: Maximum number of extension iterations to attempt
            enable_extension: Whether to perform computationally expensive contig extensions
            
        Returns:
            Best spanning read as pandas Series or fallback contig name
        """
        global failure_count
        nonlocal extended_contigs  # Access the extended_contigs from outer scope
        
        if not isinstance(merged_alignments, pd.DataFrame):
            print(f"[DEBUG] Warning: Empty alignment found between {merged_alignments[0]} and {merged_alignments[1]}. Skipping.")
            return merged_alignments[1]
        
        # Store original contigs for reference
        original_left_contig = merged_alignments['5_x'].values[0]
        original_right_contig = merged_alignments['5_y'].values[0]
        current_merged_alignments = merged_alignments
        
        print(f"[DEBUG] Finding best spanning read for {original_left_contig} and {original_right_contig}")
        
        # Initialize temp extended contig for this function call
        temp_extended_contig = None
        
        # Iterative approach instead of recursion
        for iteration in range(max_iterations):
            print(f"[DEBUG] Iteration {iteration + 1}/{max_iterations}: Finding best read...")
            
            # Define relaxation parameters for inner iterations
            relaxation_steps = [
                # Step 0: Original strict conditions
                {"gap_threshold_x": 100, "gap_threshold_y": 100, "phred_factor": 1.0},
                # Step 1: Relax gap thresholds
                {"gap_threshold_x": 1000, "gap_threshold_y": 1000, "phred_factor": 1.0},
                # Step 2: Further relax gap thresholds and slightly lower phred
                {"gap_threshold_x": 5000, "gap_threshold_y": 5000, "phred_factor": 1.0},
                # Step 3: Most relaxed conditions before extension
                {"gap_threshold_x": 10000, "gap_threshold_y": 10000, "phred_factor": 1.0}
            ]
            
            # Try progressively relaxed conditions before attempting extension
            for relax_step, params in enumerate(relaxation_steps):
                print(f"[DEBUG] Iteration {iteration + 1}, relaxation step {relax_step + 1}: gap_x<{params['gap_threshold_x']}, gap_y<{params['gap_threshold_y']}, phred*{params['phred_factor']}")
                
                adjusted_phred = int(phred_threshold * params['phred_factor'])
                
                # Try to find spanning read with current relaxed conditions
                points_of_interest = current_merged_alignments[(current_merged_alignments['1_x'] >= length_threshold)
                                                       & (current_merged_alignments['6_x'] - current_merged_alignments['8_x'] < params['gap_threshold_x'])
                                                       & (current_merged_alignments['7_y'] < params['gap_threshold_y'])
                                                       #Above conditions are variable for potential branch and bound algorithm
                                                       #Below conditions are fixed:
                                                       & (current_merged_alignments['4_x'] == '+') 
                                                       & (current_merged_alignments['4_y'] == '+')
                                                       & (current_merged_alignments['11_x'] >= adjusted_phred)
                                                       & (current_merged_alignments['11_y'] >= adjusted_phred)
                                                        
                                                       #& (merged_alignments['2_y'] > merged_alignments['3_x'])
                                                       ].copy()

                #points_of_interest['criterion'] = points_of_interest['17_y'].apply(get_clip_length)
                points_of_interest['criterion'] = points_of_interest['2_y'] - points_of_interest['3_x']

                if iteration == max_iterations - 1 and relax_step == len(relaxation_steps) - 1:  # verbose on last iteration and step
                    print("ULTRAMEGADEBUG", points_of_interest)

                # Try to find spanning read with current relaxed conditions
                try:
                    best_read = points_of_interest.loc[points_of_interest['criterion'].idxmax()]
                    print(f"[DEBUG] Found spanning read {best_read[0]} on iteration {iteration + 1}, relaxation step {relax_step + 1}")
                    print("THIS IS EXTENDED CONTIGS", extended_contigs)
                    
                    # If we found a spanning read and we have a temp extended contig from this iteration, store it
                    print(f"[DEBUG] Checking temp_extended_contig: exists={temp_extended_contig is not None}, value={temp_extended_contig}")
                    if 'temp_extended_contig' in locals() and temp_extended_contig is not None:
                        print(f"[DEBUG] Adding extended contig for {original_left_contig}")
                        extended_contigs[original_left_contig] = temp_extended_contig
                        print(f"[DEBUG] Successfully stored extended sequence for {original_left_contig} (length: {len(temp_extended_contig)}) - spanning read found!")
                    else:
                        print(f"[DEBUG] No temp extended contig to store for {original_left_contig}")
                    
                    return best_read
                except Exception as e:
                    print(f"[DEBUG] Exception in spanning read search: {e}")
                    print(f"[DEBUG] No spanning read found with relaxation step {relax_step + 1}")
                    continue  # Try next relaxation step
            
            ############## Expected gap size not used considerremove
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
            ##################

            # If all relaxation steps failed, try extension (only if not last iteration and extension is enabled)
            if iteration < max_iterations - 1 and enable_extension:
                
                # Try to extend the left contig
                print(f"[DEBUG] Iteration {iteration + 1}: No valid spanning read found between {current_merged_alignments['5_x'].values[0]} and {current_merged_alignments['5_y'].values[0]}. Attempting contig extension...")
                
                # Get alignments for the left contig to attempt extension
                left_contig_name = current_merged_alignments['5_x'].values[0]
                right_contig_name = current_merged_alignments['5_y'].values[0]


                if iteration == 0:
                    left_contig_alignments = whole_alignments[whole_alignments[5] == original_left_contig]  # Always use original left contig for extension search
                else:
                    left_contig_alignments = whole_alignments[whole_alignments[5] == left_contig_name]
                
                #^ This is iffy

                # Apply extend(left=False) logic to find extension reads
                extension_candidates = left_contig_alignments[
                    (left_contig_alignments[1] >= length_threshold) &
                    (left_contig_alignments[11] >= phred_threshold) &
                    (left_contig_alignments[4] == '+') &  # Positive strand
                    (left_contig_alignments[6] - left_contig_alignments[8] < 100)  # Close to end of contig
                ].copy()

                print(f"[DEBUG] Found {len(extension_candidates)} extension candidates for contig {left_contig_name}")
                print(extension_candidates.head())
                
                if not extension_candidates.empty:
                    # Find best extension read
                    extension_candidates['extension_criterion'] = extension_candidates[1] - extension_candidates[3]
                    best_extension = extension_candidates.loc[extension_candidates['extension_criterion'].idxmax()]
                    print(f"[DEBUG] Found best extension read {best_extension[0]} for contig {left_contig_name} with extension criterion {best_extension['extension_criterion']}")
                    
                    xprime = best_extension[2]
                    yprime = best_extension[3]
                    x = best_extension[7]
                    y = best_extension[8]

                    new_contig = contig_seqs[contig_ids.index(original_left_contig)][:y] + read_seqs[read_ids.index(best_extension[0])][yprime:]
                    
                    # Don't store extended contig yet - only store if it leads to successful spanning read
                    print(f"[DEBUG] Created extended sequence for {original_left_contig} (length: {len(new_contig)}) - will store only if successful")

                    mapping_record = [
                        SeqRecord.SeqRecord(new_contig, id=f"{original_left_contig}_iter_{iteration + 1}", description=f"Extended {original_left_contig} with read {best_extension[0]}"),
                        contigs[contig_ids.index(right_contig_name)]
                    ]

                    cur_dir = os.getcwd()
                    if not os.path.exists(os.path.join(cur_dir, "temp_contigs")):
                        os.makedirs("temp_contigs")

                    SeqIO.write(mapping_record, os.path.join(cur_dir, "temp_contigs", f"{original_left_contig}_iter_{iteration + 1}.fasta"), "fasta")

                    # Call fill_gaps.sh to fill gaps in the new contig
                    cleaned_reads_file = os.path.basename(reads_file).removesuffix('.fasta')
                    cleaned_reads_file = cleaned_reads_file.removesuffix('.fna')  # Remove .fna extension for consistency
                    cleaned_reads_file = cleaned_reads_file.removesuffix('.fastq')  # Remove .fastq extension for consistency
                    cleaned_reads_file = cleaned_reads_file.removesuffix('.fa')  # Remove .fa extension for consistency
                    
                    paf_filename = f"{cleaned_reads_file}_mapped_to_{original_left_contig}_iter_{iteration + 1}.sorted.paf"
                    if not os.path.exists(os.path.join(cur_dir, "temp_contigs", paf_filename)):
                        print(f"[DEBUG] {paf_filename} not found. Creating it now.")
                        subprocess.call(["fill_gaps.sh",
                                      os.path.join(cur_dir, "temp_contigs", f"{original_left_contig}_iter_{iteration + 1}.fasta"),
                                       f"{reads_file}",
                                        "--no_bam", "-o", "temp_contigs"
                                      ])
                    
                    # Load new alignments and prepare for next iteration
                    new_alignments = pd.read_csv(os.path.join(cur_dir, "temp_contigs", paf_filename), delimiter='\t', header=None)
                    by_contig2 = list(zip([item for _, item in new_alignments.groupby(5, sort=False)], orderings['chr']))

                    merged2 = list(map(lambda x: pd.merge(x[0][0], x[1][0], how = 'inner', on = 0), (filter(lambda x: x[0][1] == x[1][1], mit.pairwise(by_contig2)))))
                    
                    if merged2:  # If we have new merged alignments, use them for next iteration
                        current_merged_alignments = merged2[0]
                        whole_alignments = new_alignments
                        print(f"[DEBUG] Updated alignments for iteration {iteration + 2}")
                        # Store the extended contig and new_contig for potential success
                        temp_extended_contig = new_contig
                    else:
                        print(f"[DEBUG] No new merged alignments found, will retry with same data")
                        
                    print(f"[DEBUG] Found potential extension read {best_extension[0]} for contig {left_contig_name}")
                else:
                    print(f"[DEBUG] No extension candidates found for iteration {iteration + 1}")
                    # Continue to next iteration anyway, maybe with relaxed criteria
                    continue
            else:
                # This is the last iteration, all relaxation steps failed, or extension is disabled
                failure_count += 1
                extension_msg = "and extension disabled" if not enable_extension else ""
                print(f"[DEBUG] All {max_iterations} iterations and relaxation steps failed {extension_msg}. No valid spanning read found between {original_left_contig} and {original_right_contig}. Failure count: {failure_count}")
                return original_right_contig
        
        # This should never be reached due to the max_iterations check above, but just in case
        failure_count += 1
        print(f"[DEBUG] Unexpected end of iterations. Failure count: {failure_count}")
        return original_right_contig


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
    #spanning_reads = [find_best_read(item, phred_threshold=phred_threshold, length_threshold=length_threshold) for item in merged]

    spanning_reads = []
    for item in merged:
        chosen_read = find_best_read(item, whole_alignments=alignments, phred_threshold=phred_threshold, length_threshold=length_threshold, max_iterations=2, enable_extension=enable_extension)
        try:
            #print('AAAAAAAAAAAAA', item['5_x'].values[0] if not isinstance(item, str) else item, item['5_y'].values[0], chosen_read[0] if isinstance(chosen_read, pd.Series) else chosen_read)
            pass
        except:
            print(f"[DEBUGGGGGGGG] Warning: No valid spanning read found. Skipping.")
            chosen_read = item['5_y'].values[0]
        spanning_reads.append(chosen_read)

    return spanning_reads, extended_contigs




def generate_descriptions(contigs: List[str], reads, contigos, extended_contigs: dict = None) -> List[str]:
    if extended_contigs is None:
        extended_contigs = {}
        
    read_ids = [item.id for item in reads]
    read_seqs = [item.seq for item in reads]
    contig_ids = [item.id for item in contigos]
    contig_seqs = [item.seq for item in contigos]
    if isinstance(contigs, pd.Series):
        return(contigs[0])
    elif contigs is None:
        return('None')
    else:
        # Check if this contig was extended
        if contigs in extended_contigs:
            contig_length = len(extended_contigs[contigs])
            return(f"{contigs}_extended_len_{contig_length}")
        else:
            contig_length = len(contig_seqs[contig_ids.index(contigs)]) if contigs else 0
            return(f"{contigs}_len_{contig_length}")
    

@click.command()
@click.argument('alignments_file', type=click.Path(exists=True))
@click.argument('orderings_file', type=click.Path(exists=True))
@click.argument('contigs_file', type=click.Path(exists=True))
@click.argument('reads_file', type=click.Path(exists=True))
@click.option('--expected_telomere_length', default=8000, help='Expected length of telomeres (default: 8000)')
@click.option('--length_threshold', default=0, help='Minimum length of reads to consider (default: 0)')
@click.option('--phred_threshold', default=20, help='Minimum Phred quality score of reads to consider (default: equivalent to 1/100 chance of mapping misplacement <20 for telomere extension, 10 for spanning reads>)')
@click.option('--disable_extension', is_flag=True, help='Disable computationally expensive contig extension iterations (default: extension enabled)')
def merge_contigs(alignments_file, orderings_file, contigs_file, reads_file, expected_telomere_length, length_threshold, phred_threshold, disable_extension) -> dict:


    click.echo(f"[DEBUG] Starting merge_contigs with parameters:\n")
    click.echo(f"alignments_file: {alignments_file}")
    click.echo(f"orderings_file: {orderings_file}")
    click.echo(f"contigs_file: {contigs_file}")
    click.echo(f"reads_file: {reads_file}")
    click.echo(f"expected_telomere_length: {expected_telomere_length}")
    click.echo(f"length_threshold: {length_threshold}")
    click.echo(f"phred_threshold: {phred_threshold}")
    click.echo(f"extension_enabled: {not disable_extension}")

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
    spanning_reads, extended_contigs = simple_contig_span(alignments, orderings, reads, contigs, reads_file=reads_file,
                                        length_threshold=length_threshold, 
                                        phred_threshold=phred_threshold,
                                        enable_extension=not disable_extension)  #<- phred_threshold is halved for spanning reads since the mapping is paired.
    
    print(f"[DEBUG] Found {len(extended_contigs)} extended contigs: {list(extended_contigs.keys())}")
    
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
        # Extract original contig name if this is an extended contig
        if type(item) != str:
            left_contig = item['5_x']
            # If this is an extended contig (contains _iter_), extract original name
            if '_iter_' in left_contig:
                original_left_contig = left_contig.split('_iter_')[0]
            else:
                original_left_contig = left_contig
            chr = contigs_to_chr[original_left_contig]
        else:
            # If this is an extended contig name (contains _iter_), extract original name
            if '_iter_' in item:
                original_contig = item.split('_iter_')[0]
            else:
                original_contig = item
            chr = contigs_to_chr[original_contig]
        
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

    print("[DEBUG] Combining contigs for all chromosomes...")
    final_merged = {chr: combine_contigs(item, reads, contigs, extended_contigs) for chr, item in deepcopy(final_hash).items()}
    print("[DEBUG] Finished merge_contigs.")

    print("[DEBUG] Generating descriptions for contigs...")
    final_labels = {chr: ('-'.join([generate_descriptions(item, reads, contigs, extended_contigs) for item in final_hash[chr]])).split('-None-') for chr in final_hash}
    print(f"[DEBUG] Generated descriptions: {final_labels}")



    for chr, contigs in final_merged.items():
        final_merged[chr] = [item for item in final_merged[chr] if item is not None]  # Remove None values so the indices match the labels


    print("[DEBUG] Preparing output sequences...")
    fasta_record = []
    for chr, contigs in final_merged.items():
        for i, contig in enumerate(contigs):
            if isinstance(contig, str):
                # Check if this contig was extended, use extended sequence if available
                if contig in extended_contigs:
                    contig_seq = extended_contigs[contig]
                    print(f"[DEBUG] Using extended sequence for output {contig} (length: {len(contig_seq)})")
                else:
                    contig_seq = contig_seqs[contig_ids.index(contig)]
                contig_record = SeqRecord.SeqRecord(seq=contig_seq, id=f"{chr}_contig_{i+1}", description = final_labels[chr][i] if i < len(final_labels[chr]) else "")
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
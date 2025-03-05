#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
optimized_te_insertion.py - Optimized version for inserting transposable elements

This script is an optimized version of te_insertion.py that significantly improves
performance for large genomes by:
1. Using batch insertion instead of sequential insertion
2. Pre-calculating insertion positions
3. Optimizing memory usage and algorithm complexity
4. Making weight calculation optional (disabled by default)
"""

import random
import argparse
import os
import sys
import json
from collections import defaultdict
import datetime
import re
import numpy as np
# Import parameter validation and parallel processing modules
from parameter_validation import validate_te_insertion_params
from parallel_processing import get_optimal_processes

def parse_fasta(fasta_file):
    """Parse FASTA file, return sequence dictionary"""
    sequences = {}
    current_id = None
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Extract sequence ID (remove '>' symbol and other annotations)
                current_id = line[1:].split()[0]
                sequences[current_id] = []
            elif current_id:
                # Add sequence line
                sequences[current_id].append(line)
    
    # Merge all lines for each sequence
    for seq_id in sequences:
        sequences[seq_id] = ''.join(sequences[seq_id])
    
    return sequences

def parse_te_metadata(te_library):
    """Extract TE family metadata from TE library (such as preferences)"""
    metadata = {}
    
    # Default metadata
    default_metadata = {
        'gc_preference': 'neutral',  # 'high', 'low', 'neutral'
        'region_preference': []      # Preferred regions, e.g., ['telomere', 'centromere', 'arm']
    }
    
    for te_id, sequence in te_library.items():
        # Create a new metadata copy to avoid sharing references
        metadata[te_id] = {
            'gc_preference': default_metadata['gc_preference'],
            'region_preference': list(default_metadata['region_preference'])
        }
        
        # Analyze TE sequence features to infer preferences
        try:
            gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
            
            # Set preference based on GC content
            if gc_content > 0.55:
                metadata[te_id]['gc_preference'] = 'high'
            elif gc_content < 0.45:
                metadata[te_id]['gc_preference'] = 'low'
            
            # Infer region preference based on TE name (based on common naming conventions)
            te_name = te_id.lower()
            if 'telo' in te_name or re.search(r'(line|l1)', te_name):
                metadata[te_id]['region_preference'].append('telomere')
            if 'centro' in te_name or re.search(r'(alpha|sat)', te_name):
                metadata[te_id]['region_preference'].append('centromere')
            if re.search(r'(copia|gypsy|retro)', te_name):
                metadata[te_id]['region_preference'].append('arm')
        except Exception as e:
            print(f"Warning: Error processing metadata for TE '{te_id}': {e}")
    
    return metadata

def load_genome_regions(config_file):
    """Load genome region definitions"""
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
        return config.get('regions', {})
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Warning: Error loading config file: {e}")
        return {}

def calculate_gc_content_windows(genome, window_size=1000, step_size=500):
    """Calculate GC content in sliding windows across the genome"""
    gc_windows = []
    genome_length = len(genome)
    
    for i in range(0, genome_length - window_size + 1, step_size):
        window = genome[i:i+window_size]
        gc_content = (window.count('G') + window.count('C')) / window_size
        gc_windows.append((i, i+window_size, gc_content))
    
    return gc_windows

def calculate_region_map(regions, genome_length):
    """Create a region mapping for each position in the genome"""
    region_map = {}
    
    for region_name, (start, end) in regions.items():
        region_type = region_name.split('_')[0]  # Extract region type (without number)
        for pos in range(start, end + 1):
            region_map[pos] = region_type
    
    return region_map

# 修改batch_calculate_insertion_weights函数
def batch_calculate_insertion_weights(genome, regions, te_metadata, window_size=1000):
    """预计算所有TE类型的插入权重 - 优化版"""
    # 计算GC含量窗口，减少生成窗口数量
    window_size = min(5000, max(1000, len(genome) // 10000))
    step_size = window_size // 2
    
    print(f"使用动态窗口大小: {window_size}bp, 步长: {step_size}bp")
    gc_windows = calculate_gc_content_windows(genome, window_size, step_size)
    
    # 创建区域映射以快速查找
    genome_length = len(genome)
    region_map = calculate_region_map(regions, genome_length) if regions else {}
    
    # 并行计算不同TE类型的权重
    te_weights = {}
    te_ids = list(te_metadata.keys())
    
    # 创建并行执行的任务
    from functools import partial
    import multiprocessing as mp
    
    # 定义并行处理函数
    def process_te_weights(te_id, genome_length, gc_windows, region_map, te_metadata):
        weights = np.ones(genome_length, dtype=np.float32)
        metadata = te_metadata[te_id]
        gc_pref = metadata.get('gc_preference', 'neutral')
        region_prefs = metadata.get('region_preference', [])
        
        # 应用GC含量偏好 (优化模式)
        if gc_pref != 'neutral':
            for start, end, gc_content in gc_windows:
                # 简化GC乘数计算
                if gc_pref == 'high':
                    gc_multiplier = 3.0 if gc_content > 0.55 else (0.3 if gc_content < 0.45 else 1.0)
                elif gc_pref == 'low':
                    gc_multiplier = 3.0 if gc_content < 0.45 else (0.3 if gc_content > 0.55 else 1.0)
                else:
                    gc_multiplier = 1.0
                
                weights[start:end] *= gc_multiplier
        
        # 应用区域偏好 (优化为批量处理)
        if region_map and region_prefs:
            for region_name, (start, end) in regions.items():
                region_type = region_name.split('_')[0]
                weight_multiplier = 5.0 if region_type in region_prefs else 0.2
                weights[start:end+1] *= weight_multiplier
        
        return (te_id, weights)
    
    print("并行预计算TE插入权重...")
    # 使用并行处理计算所有TE类型的权重
    num_processes = get_optimal_processes(task_type='cpu_bound')
    with mp.Pool(processes=num_processes) as pool:
        process_func = partial(process_te_weights, 
                              genome_length=genome_length, 
                              gc_windows=gc_windows,
                              region_map=region_map,
                              te_metadata=te_metadata)
        
        # 使用imap_unordered处理不同的TE类型
        results = pool.imap_unordered(process_func, te_ids)
        
        # 收集结果
        for te_id, weights in results:
            te_weights[te_id] = weights
    
    return te_weights

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

def optimized_insert_te(genome, te_library, target_percentage, regions=None, te_metadata=None, num_processes=None):
    """
    Optimized version of TE insertion that uses batch processing
    
    Parameters:
    genome -- Genome sequence string
    te_library -- TE sequence dictionary {id: sequence}
    target_percentage -- Target TE percentage (0-100)
    regions -- Chromosome region definition {region_name: (start, end), ...}
    te_metadata -- TE metadata {te_id: {'gc_preference': '...', 'region_preference': [...]}, ...}
    num_processes -- Number of processes for parallel processing
    
    Returns:
    modified_genome -- Genome after TE insertion
    te_annotations -- TE annotation list [(te_id, start, end, strand, nested_in, te_seq), ...]
    """
    genome_length = len(genome)
    
    # Calculate the total TE length needed to reach target percentage
    # Formula: x / (original_length + x) = target_percentage / 100
    # where x is the total TE length to insert
    total_te_length_needed = (target_percentage * genome_length)/(100-target_percentage)
    total_te_length_needed = int(total_te_length_needed)
    
    print(f"Need to insert approximately {total_te_length_needed:,} bp of TEs to reach {target_percentage}% of genome")
    
    # Pre-calculate weights for all TE types if regions and metadata are provided
    # This is now optional and disabled by default
    te_weights = None
    if regions and te_metadata:
        print("Calculating insertion weights based on GC content and region preferences...")
        te_weights = batch_calculate_insertion_weights(genome, regions, te_metadata)
    
    # Generate all TE insertions
    insertions = []
    te_counter = defaultdict(int)
    current_te_length = 0
    
    print("Planning TE insertions...")
    
    # Continue until we reach the target length
    while current_te_length < total_te_length_needed:
        # Randomly select a TE from the library
        te_id = random.choice(list(te_library.keys()))
        te_seq = te_library[te_id]
        te_length = len(te_seq)
        
        # Choose insertion position
        if te_weights is not None and te_id in te_weights:
            # Weighted random selection based on pre-calculated weights
            weights = te_weights[te_id]
            insert_pos = random.choices(range(genome_length), weights=weights, k=1)[0]
        else:
            # Uniform random selection
            insert_pos = random.randint(0, genome_length - 1)
        
        # Choose strand
        strand = random.choice(['+', '-'])
        
        # Count instances of this TE
        te_counter[te_id] += 1
        unique_te_id = f"{te_id}_{te_counter[te_id]}"
        
        # Store insertion information
        insertions.append((unique_te_id, insert_pos, te_length, strand, te_seq))
        
        current_te_length += te_length
        
        # Show progress for every 1000 TEs
        if len(insertions) % 1000 == 0:
            progress = (current_te_length / total_te_length_needed) * 100
            print(f"Planned {len(insertions)} TE insertions ({progress:.2f}% of target)")
    
    print(f"Total planned insertions: {len(insertions)}")
    
    # Sort insertions by position (ascending)
    insertions.sort(key=lambda x: x[1])
    
    # Apply insertions in order
    # Using a list for genome manipulation is faster than string concatenation
    genome_parts = []
    current_pos = 0
    te_annotations = []
    
    print("Applying TE insertions to genome...")
    
    # Track nested TEs
    position_to_te = {}  # Maps genome positions to TE IDs for detecting nesting
    
    for i, (te_id, insert_pos, te_length, strand, te_seq) in enumerate(insertions):
        # Adjust insertion position based on previous insertions
        adjusted_pos = insert_pos
        
        # Add the genome segment before this insertion
        genome_parts.append(genome[current_pos:insert_pos])
        
        # Apply reverse complement if needed
        if strand == '-':
            inserted_seq = reverse_complement(te_seq)
        else:
            inserted_seq = te_seq
        
        # Add the TE sequence
        genome_parts.append(inserted_seq)
        
        # Find if this insertion is nested in another TE
        nested_in = None
        for pos in range(adjusted_pos, adjusted_pos + 1):  # Just check insertion point
            if pos in position_to_te:
                nested_in = position_to_te[pos]
                break
        
        # Record the annotation
        end_pos = adjusted_pos + te_length - 1
        te_annotations.append((te_id, adjusted_pos, end_pos, strand, nested_in, inserted_seq))
        
        # Mark the positions this TE occupies
        for pos in range(adjusted_pos, end_pos + 1):
            position_to_te[pos] = te_id
        
        # Update current position
        current_pos = insert_pos
        
        # Show progress for every 1000 TEs
        if (i + 1) % 1000 == 0:
            print(f"Applied {i+1}/{len(insertions)} insertions")
    
    # Add the remaining genome
    genome_parts.append(genome[current_pos:])
    
    # Join all parts to create the modified genome
    modified_genome = ''.join(genome_parts)
    
    # Calculate statistics
    inserted_te_content = sum(len(te[5]) for te in te_annotations)
    actual_te_percentage = (inserted_te_content / len(modified_genome)) * 100
    
    print(f"Finished applying {len(insertions)} TE insertions")
    print(f"Original genome length: {genome_length:,} bp")
    print(f"New genome length: {len(modified_genome):,} bp")
    print(f"Inserted TE content: {inserted_te_content:,} bp")
    print(f"Actual TE percentage: {actual_te_percentage:.2f}%")
    
    return modified_genome, te_annotations

def write_te_annotations(annotations, output_file, te_library):
    """Write TE annotations to an enhanced BED format file"""
    with open(output_file, 'w') as f:
        f.write("# TE_ID\tFamily\tStart\tEnd\tStrand\tLength\tNested_In\tOriginal_Seq\tInserted_Seq\tInsertion_Time\n")
        
        # Get current time as "insertion time"
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        for i, (te_id, start, end, strand, nested_in, te_seq) in enumerate(annotations):
            # Extract TE family (original ID, without counter)
            original_id = te_id.rsplit('_', 1)[0]
            length = end - start + 1
            
            # Get original TE sequence
            original_seq = te_library[original_id]
            if strand == '-':
                original_seq = reverse_complement(original_seq)
            
            # Write basic information
            f.write(f"{te_id}\t{original_id}\t{start}\t{end}\t{strand}\t{length}\t"
                   f"{nested_in if nested_in else 'None'}\t")
            
            # Write sequence hash values to save space
            orig_hash = hash(original_seq) % 10000000
            seq_hash = hash(te_seq) % 10000000
            f.write(f"{orig_hash}\t{seq_hash}\t{timestamp}\n")

def write_te_sequences(annotations, output_file):
    """Write inserted TE sequences to a FASTA file for subsequent analysis"""
    with open(output_file, 'w') as f:
        for te_id, start, end, strand, nested_in, te_seq in annotations:
            f.write(f">{te_id} start={start} end={end} strand={strand} nested_in={nested_in if nested_in else 'None'}\n")
            
            # Write 60 bases per line
            for i in range(0, len(te_seq), 60):
                f.write(te_seq[i:i+60] + '\n')

def write_gff3(annotations, output_file, te_library):
    """Write TE annotations to GFF3 format file"""
    with open(output_file, 'w') as f:
        f.write("##gff-version 3\n")
        
        for te_id, start, end, strand, nested_in, te_seq in annotations:
            original_id = te_id.rsplit('_', 1)[0]
            
            # Get original TE sequence
            original_seq = te_library[original_id]
            if strand == '-':
                original_seq = reverse_complement(original_seq)
            
            # Calculate similarity (simple implementation)
            min_len = min(len(original_seq), len(te_seq))
            similarity = sum(a == b for a, b in zip(original_seq[:min_len], te_seq[:min_len])) / min_len * 100
            
            attributes = f"ID={te_id};Name={original_id};similarity={similarity:.2f}"
            if nested_in:
                attributes += f";Parent={nested_in}"
            
            f.write(f"chr1\tTEsimulation\ttransposable_element\t{start+1}\t{end+1}\t.\t{strand}\t.\t{attributes}\n")

def write_fasta(sequence, output_file, chunk_size=60):
    """Write sequence to FASTA format file"""
    with open(output_file, 'w') as f:
        f.write(">chr1 with inserted TEs\n")
        # Write chunk_size bases per line
        for i in range(0, len(sequence), chunk_size):
            f.write(sequence[i:i+chunk_size] + '\n')

def main():
    parser = argparse.ArgumentParser(description='Insert transposable elements into a genome considering chromosome structure and TE preferences')
    parser.add_argument('-g', '--genome', type=str, required=True,
                        help='Input genome FASTA file')
    parser.add_argument('-t', '--te_library', type=str, required=True,
                        help='TE library FASTA file')
    parser.add_argument('-p', '--percentage', type=float, default=30.0,
                        help='Target TE percentage (default: 30.0%)')
    parser.add_argument('-o', '--output_prefix', type=str, default='te_inserted',
                        help='Output file prefix (default: te_inserted)')
    parser.add_argument('-c', '--config', type=str, default=None,
                        help='Genome configuration file (JSON format, with region definitions)')
    parser.add_argument('--use_weights', action='store_true',
                        help='Enable weight calculation for biased TE insertion (slower but more realistic)')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed for reproducibility (default: random)')
    parser.add_argument('--num_processes', type=int, default=None,
                        help='Number of processes for parallel processing (default: 75% of CPU cores)')
    
    args = parser.parse_args()
    
    # Parameter validation
    validation_result = validate_te_insertion_params(args)
    validation_result.print_messages()
    validation_result.exit_if_invalid()
    
    # Set random seed
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)
    
    # Set number of processes
    num_processes = args.num_processes if args.num_processes else get_optimal_processes()
    print(f"Using {num_processes} processes for parallel computation")
    
    # Parse input files
    print("Parsing genome and TE library...")
    genome_dict = parse_fasta(args.genome)
    te_library = parse_fasta(args.te_library)
    
    if not genome_dict:
        sys.exit("Error: Genome file is empty or format is incorrect")
    if not te_library:
        sys.exit("Error: TE library is empty or format is incorrect")
    
    # Get the first chromosome sequence
    chromosome_id = list(genome_dict.keys())[0]
    genome_sequence = genome_dict[chromosome_id]
    
    # Load genome region definitions (if provided)
    regions = None
    if args.config:
        if os.path.exists(args.config):
            regions = load_genome_regions(args.config)
            print(f"Loaded genome region definitions, containing {len(regions)} regions")
        else:
            print(f"Warning: Configuration file {args.config} does not exist, will use uniform insertion model")
    
    # Analyze TE metadata only if using weights
    te_metadata = None
    if args.use_weights and regions:
        te_metadata = parse_te_metadata(te_library)
        print(f"Analyzed TE metadata, total {len(te_metadata)} TE families")
    else:
        print("Using fast uniform insertion (no weight calculation)")
        # Set regions to None to ensure uniform insertion
        regions = None
    
    print(f"Genome length: {len(genome_sequence):,} bp")
    print(f"Number of TEs in the library: {len(te_library)}")
    print(f"Target TE percentage: {args.percentage}%")
    
    # Create output directory
    os.makedirs(os.path.dirname(os.path.abspath(args.output_prefix)), exist_ok=True)
    
    # Insert TEs using the optimized method
    print("Starting TE insertion (optimized method)...")
    start_time = datetime.datetime.now()
    modified_genome, te_annotations = optimized_insert_te(
        genome_sequence, te_library, args.percentage, regions, te_metadata, num_processes)
    end_time = datetime.datetime.now()
    
    print(f"Insertion completed in {(end_time - start_time).total_seconds():.2f} seconds")
    
    # Save results
    genome_output = f"{args.output_prefix}_genome.fa"
    bed_output = f"{args.output_prefix}_te_annotations.bed"
    gff_output = f"{args.output_prefix}_te_annotations.gff3"
    seq_output = f"{args.output_prefix}_te_sequences.fa"
    
    print("Saving results...")
    write_fasta(modified_genome, genome_output)
    write_te_annotations(te_annotations, bed_output, te_library)
    write_gff3(te_annotations, gff_output, te_library)
    write_te_sequences(te_annotations, seq_output)
    
    print(f"Done! Inserted {len(te_annotations)} TEs")
    print(f"Modified genome saved to: {genome_output}")
    print(f"TE annotations (BED format) saved to: {bed_output}")
    print(f"TE annotations (GFF3 format) saved to: {gff_output}")
    print(f"TE sequences saved to: {seq_output}")

if __name__ == "__main__":
    main()

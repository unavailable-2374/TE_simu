#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
optimized_te_insertion.py - 在基因组中插入转座元件的优化版本

此脚本是te_insertion.py的优化版本，主要改进：
1. 使用批量插入而非顺序插入
2. 预先计算插入位置
3. 优化内存使用和算法复杂度
4. 权重计算设为可选（默认禁用以提高性能）
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
# 导入参数验证和并行处理模块
from parameter_validation import validate_te_insertion_params
from parallel_processing import get_optimal_processes

def parse_fasta(fasta_file):
    """解析FASTA文件，返回序列字典"""
    sequences = {}
    current_id = None
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # 提取序列ID（去除'>'符号和其他注释）
                current_id = line[1:].split()[0]
                sequences[current_id] = []
            elif current_id:
                # 添加序列行
                sequences[current_id].append(line)
    
    # 合并每个序列的所有行
    for seq_id in sequences:
        sequences[seq_id] = ''.join(sequences[seq_id])
    
    return sequences

def parse_te_metadata(te_library):
    """从TE文库中提取TE家族元数据（如倾向性）"""
    metadata = {}
    
    # 默认元数据
    default_metadata = {
        'gc_preference': 'neutral',  # 'high', 'low', 'neutral'
        'region_preference': []      # 偏好的区域，如 ['telomere', 'centromere', 'arm']
    }
    
    for te_id, sequence in te_library.items():
        # 创建一个新的元数据副本以避免共享引用
        metadata[te_id] = {
            'gc_preference': default_metadata['gc_preference'],
            'region_preference': list(default_metadata['region_preference'])
        }
        
        # 分析TE序列特征来推断其偏好
        try:
            gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
            
            # 根据GC含量设置偏好
            if gc_content > 0.55:
                metadata[te_id]['gc_preference'] = 'high'
            elif gc_content < 0.45:
                metadata[te_id]['gc_preference'] = 'low'
            
            # 根据TE名称推断区域偏好（基于常见命名规则）
            te_name = te_id.lower()
            if 'telo' in te_name or re.search(r'(line|l1)', te_name):
                metadata[te_id]['region_preference'].append('telomere')
            if 'centro' in te_name or re.search(r'(alpha|sat)', te_name):
                metadata[te_id]['region_preference'].append('centromere')
            if re.search(r'(copia|gypsy|retro)', te_name):
                metadata[te_id]['region_preference'].append('arm')
        except Exception as e:
            print(f"警告: 处理TE '{te_id}' 元数据时出错: {e}")
    
    return metadata

def load_genome_regions(config_file):
    """加载基因组区域定义"""
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
        return config.get('regions', {})
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"警告: 加载配置文件出错: {e}")
        return {}

def calculate_gc_content_windows(genome, window_size=1000, step_size=500):
    """计算基因组中滑动窗口的GC含量"""
    gc_windows = []
    genome_length = len(genome)
    
    for i in range(0, genome_length - window_size + 1, step_size):
        window = genome[i:i+window_size]
        gc_content = (window.count('G') + window.count('C')) / window_size
        gc_windows.append((i, i+window_size, gc_content))
    
    return gc_windows

def calculate_region_map(regions, genome_length):
    """为基因组中的每个位置创建区域映射"""
    region_map = {}
    
    for region_name, (start, end) in regions.items():
        region_type = region_name.split('_')[0]  # 提取区域类型（不含编号）
        for pos in range(start, end + 1):
            region_map[pos] = region_type
    
    return region_map

def batch_calculate_insertion_weights(genome, regions, te_metadata, window_size=1000):
    """
    预先计算所有TE类型在基因组所有位置的权重
    
    返回:
    te_weights -- 字典: {te_id: [所有位置的权重]}
    """
    # 计算GC含量窗口
    gc_windows = calculate_gc_content_windows(genome, window_size)
    
    # 创建区域映射以快速查找
    genome_length = len(genome)
    region_map = calculate_region_map(regions, genome_length) if regions else {}
    
    # 初始化所有TE的权重
    te_weights = {}
    
    print("正在预计算所有TE类型的插入权重...")
    
    for te_id, metadata in te_metadata.items():
        # 创建一个初始化为1.0的权重数组
        weights = np.ones(genome_length, dtype=np.float32)
        
        gc_pref = metadata.get('gc_preference', 'neutral')
        region_prefs = metadata.get('region_preference', [])
        
        # 使用窗口应用GC含量偏好
        if gc_pref != 'neutral':
            for start, end, gc_content in gc_windows:
                gc_multiplier = 1.0
                
                if gc_pref == 'high' and gc_content > 0.55:
                    gc_multiplier = 3.0
                elif gc_pref == 'high' and gc_content < 0.45:
                    gc_multiplier = 0.3
                elif gc_pref == 'low' and gc_content < 0.45:
                    gc_multiplier = 3.0
                elif gc_pref == 'low' and gc_content > 0.55:
                    gc_multiplier = 0.3
                
                weights[start:end] *= gc_multiplier
        
        # 应用区域偏好
        if region_map and region_prefs:
            for pos in range(genome_length):
                if pos in region_map:
                    region_type = region_map[pos]
                    weight_multiplier = 5.0 if region_type in region_prefs else 0.2
                    weights[pos] *= weight_multiplier
        
        te_weights[te_id] = weights
    
    return te_weights

def reverse_complement(sequence):
    """返回DNA序列的反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

def optimized_insert_te(genome, te_library, target_percentage, regions=None, te_metadata=None, num_processes=None):
    """
    使用批处理的TE插入优化版本
    
    参数:
    genome -- 基因组序列字符串
    te_library -- TE序列字典 {id: sequence}
    target_percentage -- 目标TE占比（0-100）
    regions -- 染色体区域定义 {region_name: (start, end), ...}
    te_metadata -- TE元数据 {te_id: {'gc_preference': '...', 'region_preference': [...]}, ...}
    num_processes -- 并行处理的进程数量
    
    返回:
    modified_genome -- 插入TE后的基因组
    te_annotations -- TE注释列表 [(te_id, start, end, strand, nested_in, te_seq), ...]
    """
    genome_length = len(genome)
    
    # 计算达到目标百分比所需的总TE长度
    # 公式: x / (original_length + x) = target_percentage / 100
    # 其中x是要插入的总TE长度
    total_te_length_needed = (target_percentage * genome_length) / (100 - target_percentage)
    total_te_length_needed = int(total_te_length_needed)
    
    print(f"需要插入约 {total_te_length_needed:,} bp的TE以达到基因组的{target_percentage}%")
    
    # 如果提供了区域和元数据，预先计算所有TE类型的权重
    # 现在这是可选的，默认禁用
    te_weights = None
    if regions and te_metadata:
        print("正在根据GC含量和区域偏好计算插入权重...")
        te_weights = batch_calculate_insertion_weights(genome, regions, te_metadata)
    
    # 生成所有TE插入
    insertions = []
    te_counter = defaultdict(int)
    current_te_length = 0
    
    print("正在规划TE插入...")
    
    # 继续，直到达到目标长度
    while current_te_length < total_te_length_needed:
        # 从TE文库中随机选择一个TE
        te_id = random.choice(list(te_library.keys()))
        te_seq = te_library[te_id]
        te_length = len(te_seq)
        
        # 选择插入位置
        if te_weights is not None and te_id in te_weights:
            # 基于预先计算的权重进行加权随机选择
            weights = te_weights[te_id]
            insert_pos = random.choices(range(genome_length), weights=weights, k=1)[0]
        else:
            # 均匀随机选择
            insert_pos = random.randint(0, genome_length - 1)
        
        # 选择链方向
        strand = random.choice(['+', '-'])
        
        # 计数该TE的实例
        te_counter[te_id] += 1
        unique_te_id = f"{te_id}_{te_counter[te_id]}"
        
        # 存储插入信息
        insertions.append((unique_te_id, insert_pos, te_length, strand, te_seq))
        
        current_te_length += te_length
        
        # 每1000个TE显示一次进度
        if len(insertions) % 1000 == 0:
            progress = (current_te_length / total_te_length_needed) * 100
            print(f"已规划 {len(insertions)} 个TE插入 ({progress:.2f}% 的目标)")
    
    print(f"总计划插入数量: {len(insertions)}")
    
    # 按位置排序插入（升序）
    insertions.sort(key=lambda x: x[1])
    
    # 按顺序应用插入
    # 使用列表进行基因组操作比字符串连接更快
    genome_parts = []
    current_pos = 0
    te_annotations = []
    
    print("正在应用TE插入到基因组...")
    
    # 跟踪嵌套的TE
    position_to_te = {}  # 将基因组位置映射到TE ID，用于检测嵌套
    
    for i, (te_id, insert_pos, te_length, strand, te_seq) in enumerate(insertions):
        # 根据之前的插入调整插入位置
        adjusted_pos = insert_pos
        
        # 添加此插入前的基因组段
        genome_parts.append(genome[current_pos:insert_pos])
        
        # 如果需要，应用反向互补
        if strand == '-':
            inserted_seq = reverse_complement(te_seq)
        else:
            inserted_seq = te_seq
        
        # 添加TE序列
        genome_parts.append(inserted_seq)
        
        # 检查该插入是否嵌套在另一个TE中
        nested_in = None
        for pos in range(adjusted_pos, adjusted_pos + 1):  # 只检查插入点
            if pos in position_to_te:
                nested_in = position_to_te[pos]
                break
        
        # 记录注释
        end_pos = adjusted_pos + te_length - 1
        te_annotations.append((te_id, adjusted_pos, end_pos, strand, nested_in, inserted_seq))
        
        # 标记该TE占据的位置
        for pos in range(adjusted_pos, end_pos + 1):
            position_to_te[pos] = te_id
        
        # 更新当前位置
        current_pos = insert_pos
        
        # 每1000个TE显示一次进度
        if (i + 1) % 1000 == 0:
            print(f"已应用 {i+1}/{len(insertions)} 个插入")
    
    # 添加剩余的基因组
    genome_parts.append(genome[current_pos:])
    
    # 连接所有部分以创建修改后的基因组
    modified_genome = ''.join(genome_parts)
    
    # 计算统计数据
    inserted_te_content = sum(len(te[5]) for te in te_annotations)
    actual_te_percentage = (inserted_te_content / len(modified_genome)) * 100
    
    print(f"已完成应用 {len(insertions)} 个TE插入")
    print(f"原始基因组长度: {genome_length:,} bp")
    print(f"新基因组长度: {len(modified_genome):,} bp")
    print(f"插入的TE内容: {inserted_te_content:,} bp")
    print(f"实际TE百分比: {actual_te_percentage:.2f}%")
    
    return modified_genome, te_annotations

def write_te_annotations(annotations, output_file, te_library):
    """将TE注释写入增强的BED格式文件"""
    with open(output_file, 'w') as f:
        f.write("# TE_ID\tFamily\tStart\tEnd\tStrand\tLength\tNested_In\tOriginal_Seq\tInserted_Seq\tInsertion_Time\n")
        
        # 获取当前时间作为"插入时间"
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        for i, (te_id, start, end, strand, nested_in, te_seq) in enumerate(annotations):
            # 提取TE家族（原始ID，不包含计数器）
            original_id = te_id.rsplit('_', 1)[0]
            length = end - start + 1
            
            # 获取原始TE序列
            original_seq = te_library[original_id]
            if strand == '-':
                original_seq = reverse_complement(original_seq)
            
            # 写入基本信息
            f.write(f"{te_id}\t{original_id}\t{start}\t{end}\t{strand}\t{length}\t"
                   f"{nested_in if nested_in else 'None'}\t")
            
            # 写入序列哈希值，以节省空间
            orig_hash = hash(original_seq) % 10000000
            seq_hash = hash(te_seq) % 10000000
            f.write(f"{orig_hash}\t{seq_hash}\t{timestamp}\n")

def write_te_sequences(annotations, output_file):
    """将插入的TE序列写入FASTA文件，用于后续分析"""
    with open(output_file, 'w') as f:
        for te_id, start, end, strand, nested_in, te_seq in annotations:
            f.write(f">{te_id} start={start} end={end} strand={strand} nested_in={nested_in if nested_in else 'None'}\n")
            
            # 每行写入60个碱基
            for i in range(0, len(te_seq), 60):
                f.write(te_seq[i:i+60] + '\n')

def write_gff3(annotations, output_file, te_library):
    """将TE注释写入GFF3格式文件"""
    with open(output_file, 'w') as f:
        f.write("##gff-version 3\n")
        
        for te_id, start, end, strand, nested_in, te_seq in annotations:
            original_id = te_id.rsplit('_', 1)[0]
            
            # 获取原始TE序列
            original_seq = te_library[original_id]
            if strand == '-':
                original_seq = reverse_complement(original_seq)
            
            # 计算相似度（简单实现）
            min_len = min(len(original_seq), len(te_seq))
            similarity = sum(a == b for a, b in zip(original_seq[:min_len], te_seq[:min_len])) / min_len * 100
            
            attributes = f"ID={te_id};Name={original_id};similarity={similarity:.2f}"
            if nested_in:
                attributes += f";Parent={nested_in}"
            
            f.write(f"chr1\tTEsimulation\ttransposable_element\t{start+1}\t{end+1}\t.\t{strand}\t.\t{attributes}\n")

def write_fasta(sequence, output_file, chunk_size=60):
    """将序列写入FASTA格式文件"""
    with open(output_file, 'w') as f:
        f.write(">chr1 with inserted TEs\n")
        # 每行写入chunk_size个碱基
        for i in range(0, len(sequence), chunk_size):
            f.write(sequence[i:i+chunk_size] + '\n')

def main():
    parser = argparse.ArgumentParser(description='在基因组中插入转座元件，考虑染色体结构和TE偏好')
    parser.add_argument('-g', '--genome', type=str, required=True,
                        help='输入基因组FASTA文件')
    parser.add_argument('-t', '--te_library', type=str, required=True,
                        help='TE文库FASTA文件')
    parser.add_argument('-p', '--percentage', type=float, default=30.0,
                        help='目标TE占比 (默认: 30.0%)')
    parser.add_argument('-o', '--output_prefix', type=str, default='te_inserted',
                        help='输出文件前缀 (默认: te_inserted)')
    parser.add_argument('-c', '--config', type=str, default=None,
                        help='基因组配置文件 (JSON格式，包含区域定义)')
    parser.add_argument('--use_weights', action='store_true',
                        help='启用权重计算以进行有偏TE插入 (速度较慢但更真实)')
    parser.add_argument('--seed', type=int, default=None,
                        help='随机数种子，用于结果复现 (默认: 随机)')
    parser.add_argument('--num_processes', type=int, default=None,
                        help='并行处理的进程数量 (默认: CPU核心数的75%)')
    
    args = parser.parse_args()
    
    # 参数验证
    validation_result = validate_te_insertion_params(args)
    validation_result.print_messages()
    validation_result.exit_if_invalid()
    
    # 设置随机种子
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)
    
    # 设置并行处理进程数
    num_processes = args.num_processes if args.num_processes else get_optimal_processes()
    print(f"使用 {num_processes} 个进程进行并行计算")
    
    # 解析输入文件
    print("解析基因组和TE文库...")
    genome_dict = parse_fasta(args.genome)
    te_library = parse_fasta(args.te_library)
    
    if not genome_dict:
        sys.exit("错误: 基因组文件为空或格式不正确")
    if not te_library:
        sys.exit("错误: TE文库为空或格式不正确")
    
    # 获取第一条染色体序列
    chromosome_id = list(genome_dict.keys())[0]
    genome_sequence = genome_dict[chromosome_id]
    
    # 加载基因组区域定义（如果提供）
    regions = None
    if args.config:
        if os.path.exists(args.config):
            regions = load_genome_regions(args.config)
            print(f"已加载基因组区域定义，包含 {len(regions)} 个区域")
        else:
            print(f"警告: 配置文件 {args.config} 不存在，将使用均匀插入模型")
    
    # 仅在使用权重时分析TE元数据
    te_metadata = None
    if args.use_weights and regions:
        te_metadata = parse_te_metadata(te_library)
        print(f"已分析TE元数据，共 {len(te_metadata)} 个TE家族")
    else:
        print("使用快速均匀插入（无权重计算）")
        # 将regions设为None以确保均匀插入
        regions = None
    
    print(f"基因组长度: {len(genome_sequence):,} bp")
    print(f"TE文库中TE数量: {len(te_library)}")
    print(f"目标TE占比: {args.percentage}%")
    
    # 创建输出目录
    os.makedirs(os.path.dirname(os.path.abspath(args.output_prefix)), exist_ok=True)
    
    # 使用优化方法插入TE
    print("开始TE插入（优化方法）...")
    start_time = datetime.datetime.now()
    modified_genome, te_annotations = optimized_insert_te(
        genome_sequence, te_library, args.percentage, regions, te_metadata, num_processes)
    end_time = datetime.datetime.now()
    
    print(f"插入完成，用时 {(end_time - start_time).total_seconds():.2f} 秒")
    
    # 保存结果
    genome_output = f"{args.output_prefix}_genome.fa"
    bed_output = f"{args.output_prefix}_te_annotations.bed"
    gff_output = f"{args.output_prefix}_te_annotations.gff3"
    seq_output = f"{args.output_prefix}_te_sequences.fa"
    
    print("保存结果...")
    write_fasta(modified_genome, genome_output)
    write_te_annotations(te_annotations, bed_output, te_library)
    write_gff3(te_annotations, gff_output, te_library)
    write_te_sequences(te_annotations, seq_output)
    
    print(f"完成! 插入了 {len(te_annotations)} 个TE")
    print(f"修改后的基因组保存至: {genome_output}")
    print(f"TE注释 (BED格式) 保存至: {bed_output}")
    print(f"TE注释 (GFF3格式) 保存至: {gff_output}")
    print(f"TE序列保存至: {seq_output}")

if __name__ == "__main__":
    main()

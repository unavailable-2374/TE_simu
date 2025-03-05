#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_base_genome.py - 生成初始随机基因组序列

此脚本生成一个具有特定GC含量和染色体结构特征的序列，作为模拟TE扩增的初始基因组。
支持GC含量梯度和特殊染色体区域定义。
"""

import random
import argparse
import os
import json
import numpy as np
from collections import OrderedDict
# 导入参数验证模块
from parameter_validation import validate_genome_generation_params

def generate_sequence_with_gc(length, gc_content):
    """生成指定GC含量的随机DNA序列"""
    g_c_count = int(length * gc_content)
    a_t_count = length - g_c_count
    
    # 生成G和C，数量由gc_content决定
    g_count = g_c_count // 2
    c_count = g_c_count - g_count
    
    # 生成A和T，数量由剩余决定
    a_count = a_t_count // 2
    t_count = a_t_count - a_count
    
    # 合并所有碱基并随机打乱
    bases = ['G'] * g_count + ['C'] * c_count + ['A'] * a_count + ['T'] * t_count
    random.shuffle(bases)
    
    return ''.join(bases)

def generate_gc_gradient_sequence(length, min_gc, max_gc, window_size=1000):
    """生成具有GC含量梯度的序列"""
    # 计算需要多少个窗口
    num_windows = (length + window_size - 1) // window_size
    
    # 创建GC含量梯度
    gc_gradient = np.linspace(min_gc, max_gc, num_windows)
    
    # 为每个窗口生成序列
    sequence = ""
    for i in range(num_windows):
        window_length = min(window_size, length - i * window_size)
        if window_length <= 0:
            break
        window_gc = gc_gradient[i]
        window_seq = generate_sequence_with_gc(window_length, window_gc)
        sequence += window_seq
    
    return sequence

def define_chromosome_regions(length, telomere_size=30000, centromere_size=100000):
    """定义染色体的特殊区域"""
    regions = OrderedDict()
    
    # 定义端粒区域（两端）
    regions['telomere_1'] = (0, telomere_size - 1)
    regions['telomere_2'] = (length - telomere_size, length - 1)
    
    # 定义着丝粒区域（居中）
    centromere_start = (length - centromere_size) // 2
    regions['centromere'] = (centromere_start, centromere_start + centromere_size - 1)
    
    # 定义普通染色体臂区域
    regions['arm_1'] = (telomere_size, centromere_start - 1)
    regions['arm_2'] = (centromere_start + centromere_size, length - telomere_size - 1)
    
    return regions

def generate_structured_genome(length, regions, gc_content=0.5, gc_gradient=False, min_gc=0.35, max_gc=0.65):
    """生成具有结构的基因组序列"""
    # 初始化空序列
    genome = ['N'] * length
    
    # 为每个区域生成序列
    for region_name, (start, end) in regions.items():
        region_length = end - start + 1
        
        # 根据区域类型调整GC含量
        region_gc = gc_content
        if 'telomere' in region_name:
            # 端粒区域通常AT富集
            region_gc = 0.3
        elif 'centromere' in region_name:
            # 着丝粒区域通常AT富集
            region_gc = 0.35
        elif gc_gradient and 'arm' in region_name:
            # 染色体臂可以有GC梯度
            if region_name == 'arm_1':
                region_seq = generate_gc_gradient_sequence(region_length, min_gc, max_gc)
            else:  # arm_2
                region_seq = generate_gc_gradient_sequence(region_length, max_gc, min_gc)
            genome[start:end+1] = list(region_seq)
            continue
        
        # 生成区域序列
        region_seq = generate_sequence_with_gc(region_length, region_gc)
        genome[start:end+1] = list(region_seq)
    
    return ''.join(genome)

def write_fasta(sequence, output_file, chunk_size=60):
    """将序列写入FASTA格式文件"""
    with open(output_file, 'w') as f:
        f.write(">chr1 structured random genome\n")
        # 每行写入chunk_size个碱基
        for i in range(0, len(sequence), chunk_size):
            f.write(sequence[i:i+chunk_size] + '\n')

def write_regions_bed(regions, output_file):
    """将染色体区域信息写入BED格式文件"""
    with open(output_file, 'w') as f:
        f.write("# 染色体区域定义\n")
        f.write("# Region\tStart\tEnd\tLength\n")
        
        for region_name, (start, end) in regions.items():
            length = end - start + 1
            f.write(f"{region_name}\t{start}\t{end}\t{length}\n")

def calculate_gc_content(sequence, window_size=10000, step_size=5000):
    """计算序列的GC含量分布"""
    gc_content = []
    positions = []
    
    for i in range(0, len(sequence) - window_size + 1, step_size):
        window = sequence[i:i+window_size]
        gc = (window.count('G') + window.count('C')) / window_size
        gc_content.append(gc)
        positions.append(i + window_size // 2)  # 窗口中心位置
    
    return positions, gc_content

def write_gc_content(positions, gc_content, output_file):
    """将GC含量分布写入文件"""
    with open(output_file, 'w') as f:
        f.write("# 染色体GC含量分布\n")
        f.write("# Position\tGC_Content\n")
        
        for pos, gc in zip(positions, gc_content):
            f.write(f"{pos}\t{gc:.4f}\n")

def main():
    parser = argparse.ArgumentParser(description='生成具有结构特征的随机基因组序列')
    parser.add_argument('-l', '--length', type=int, default=50000000,
                        help='基因组长度 (默认: 50,000,000 bp)')
    parser.add_argument('-o', '--output', type=str, default='base_genome.fa',
                        help='输出文件名 (默认: base_genome.fa)')
    parser.add_argument('--gc', type=float, default=0.5,
                        help='基础GC含量 (默认: 0.5)')
    parser.add_argument('--gc_gradient', action='store_true',
                        help='是否添加GC含量梯度')
    parser.add_argument('--min_gc', type=float, default=0.35,
                        help='GC梯度最小值 (默认: 0.35)')
    parser.add_argument('--max_gc', type=float, default=0.65,
                        help='GC梯度最大值 (默认: 0.65)')
    parser.add_argument('--telomere_size', type=int, default=30000,
                        help='端粒区域大小 (默认: 30,000 bp)')
    parser.add_argument('--centromere_size', type=int, default=100000,
                        help='着丝粒区域大小 (默认: 100,000 bp)')
    parser.add_argument('--seed', type=int, default=None,
                        help='随机数种子，用于结果复现 (默认: 随机)')
    parser.add_argument('--num_processes', type=int, default=None,
                        help='并行处理的进程数量 (默认: CPU核心数的75%)')
    
    args = parser.parse_args()
    
    # 参数验证
    validation_result = validate_genome_generation_params(args)
    validation_result.print_messages()
    validation_result.exit_if_invalid()
    
    # 设置随机种子
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)
    
    print(f"生成长度为 {args.length:,} bp 的结构化随机基因组...")
    
    # 定义染色体区域
    regions = define_chromosome_regions(args.length, args.telomere_size, args.centromere_size)
    
    # 生成结构化序列
    sequence = generate_structured_genome(
        args.length, regions, args.gc, args.gc_gradient, args.min_gc, args.max_gc)
    
    # 创建输出目录（如果不存在）
    output_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(output_dir, exist_ok=True)
    
    # 写入FASTA文件
    write_fasta(sequence, args.output)
    
    # 基本输出文件名（不含扩展名）
    output_base = os.path.splitext(args.output)[0]
    
    # 写入区域定义
    regions_file = f"{output_base}_regions.bed"
    write_regions_bed(regions, regions_file)
    
    # 计算并写入GC含量分布
    positions, gc_content = calculate_gc_content(sequence)
    gc_file = f"{output_base}_gc_content.txt"
    write_gc_content(positions, gc_content, gc_file)
    
    # 创建JSON配置文件，用于其他脚本
    config = {
        'genome_length': args.length,
        'regions': {name: [start, end] for name, (start, end) in regions.items()},
        'gc_content': args.gc,
        'gc_gradient': args.gc_gradient,
        'min_gc': args.min_gc,
        'max_gc': args.max_gc,
        'telomere_size': args.telomere_size,
        'centromere_size': args.centromere_size
    }
    
    config_file = f"{output_base}_config.json"
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"随机基因组已生成并保存至 {args.output}")
    print(f"染色体区域定义已保存至 {regions_file}")
    print(f"GC含量分布已保存至 {gc_file}")
    print(f"配置信息已保存至 {config_file}")

if __name__ == "__main__":
    main()

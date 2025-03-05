#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
visualization_tools.py - 染色体结构与TE分布可视化工具

此脚本生成直观的可视化图表，展示染色体结构特征、TE分布和进化情况。
"""

import os
import sys
import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter
from collections import defaultdict, Counter

def parse_fasta(fasta_file):
    """解析FASTA文件，返回序列"""
    sequence = ""
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            sequence += line
    return sequence

def load_genome_config(config_file):
    """加载基因组配置文件"""
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
        return config
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"警告: 无法加载配置文件 {config_file}: {e}")
        return None

def parse_annotations(annotation_file):
    """解析BED/GFF格式的TE注释文件"""
    annotations = []
    
    with open(annotation_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue
            
            # 检测文件格式（BED或GFF）
            if len(fields) >= 9 and fields[2] == 'transposable_element':  # GFF3格式
                start = int(fields[3]) - 1  # GFF是1-based，转换为0-based
                end = int(fields[4]) - 1
                strand = fields[6]
                
                # 解析属性字段
                attributes = {}
                for attr in fields[8].split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attributes[key] = value
                
                te_id = attributes.get('ID', 'unknown')
                family = attributes.get('Name', 'unknown')
                similarity = float(attributes.get('similarity', '0'))
                evol_time = float(attributes.get('evolutionary_time', '0'))
                completeness = attributes.get('completeness', 'unknown')
                
                annotations.append({
                    'te_id': te_id,
                    'family': family,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'similarity': similarity,
                    'evolutionary_time': evol_time,
                    'completeness': completeness
                })
            
            else:  # 假设是BED格式
                te_id = fields[0]
                family = fields[1]
                start = int(fields[2])
                end = int(fields[3])
                strand = fields[4]
                
                # 提取可能的相似度信息（如果有）
                similarity = 0
                evol_time = 0
                
                annotations.append({
                    'te_id': te_id,
                    'family': family,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'similarity': similarity,
                    'evolutionary_time': evol_time
                })
    
    return annotations

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

def calculate_te_density(annotations, genome_length, window_size=100000, step_size=50000):
    """计算TE密度分布"""
    density = []
    positions = []
    
    for i in range(0, genome_length - window_size + 1, step_size):
        window_end = i + window_size
        # 计算此窗口中的TE覆盖碱基数
        covered_bases = 0
        for te in annotations:
            # 检查TE是否与窗口重叠
            if not (te['end'] < i or te['start'] > window_end):
                overlap_start = max(i, te['start'])
                overlap_end = min(window_end, te['end'])
                covered_bases += (overlap_end - overlap_start + 1)
        
        density.append(covered_bases / window_size * 100)  # 百分比
        positions.append(i + window_size // 2)  # 窗口中心位置
    
    return positions, density

def visualize_chromosome_structure(sequence, annotations, config, output_prefix):
    """生成染色体结构与TE分布可视化图"""
    genome_length = len(sequence)
    
    # 创建输出目录
    os.makedirs(os.path.dirname(os.path.abspath(output_prefix)), exist_ok=True)
    
    # 计算GC含量
    gc_positions, gc_content = calculate_gc_content(sequence)
    
    # 计算TE密度
    te_positions, te_density = calculate_te_density(annotations, genome_length)
    
    # 提取染色体区域信息（如果有）
    chromosome_regions = {}
    if config and 'regions' in config:
        chromosome_regions = config['regions']
    
    # 创建染色体结构与TE分布图
    plt.figure(figsize=(15, 10))
    
    # 1. 绘制染色体示意图
    plt.subplot(4, 1, 1)
    plt.title('染色体结构与TE分布', fontsize=14)
    
    # 绘制染色体主体
    plt.axhline(y=1, xmin=0, xmax=1, color='grey', linewidth=10, alpha=0.5)
    
    # 如果有区域信息，绘制特殊区域
    if chromosome_regions:
        region_colors = {
            'telomere': 'red',
            'centromere': 'blue',
            'arm': 'grey'
        }
        
        for region_name, (start, end) in chromosome_regions.items():
            region_type = region_name.split('_')[0]  # 提取区域类型（不含编号）
            color = region_colors.get(region_type, 'purple')
            plt.axhline(y=1, xmin=start/genome_length, xmax=(end+1)/genome_length, 
                       color=color, linewidth=10, alpha=0.8)
    
    # 绘制TE分布
    for te in annotations:
        color = 'green'
        alpha = 0.7
        
        # 如果有进化时间信息，使用颜色梯度
        if hasattr(te, 'evolutionary_time') and te['evolutionary_time'] > 0:
            # 颜色从红（古老）到绿（年轻）
            evol_time = min(te['evolutionary_time'], 10) / 10  # 规范化到0-1范围
            r = evol_time
            g = 1 - evol_time
            b = 0
            color = (r, g, b)
        
        # 如果有完整性信息，调整透明度
        if hasattr(te, 'completeness'):
            if te['completeness'] == 'fragmented':
                alpha = 0.4
        
        plt.axhline(y=1.1, xmin=te['start']/genome_length, 
                   xmax=(te['end']+1)/genome_length, 
                   color=color, linewidth=3, alpha=alpha)
    
    plt.ylim(0.5, 1.5)
    plt.xlim(0, genome_length)
    plt.ylabel('染色体结构')
    plt.gca().xaxis.set_major_formatter(
        FuncFormatter(lambda x, _: f'{x/1e6:.1f}Mb'))
    plt.gca().set_yticks([])
    
    # 创建图例
    legend_elements = []
    if chromosome_regions:
        for region_type, color in region_colors.items():
            if any(region_type in name for name in chromosome_regions.keys()):
                legend_elements.append(
                    mpatches.Patch(color=color, alpha=0.8, label=region_type.capitalize()))
    
    legend_elements.append(mpatches.Patch(color='green', alpha=0.7, label='TE'))
    plt.legend(handles=legend_elements, loc='upper right')
    
    # 2. 绘制GC含量曲线
    plt.subplot(4, 1, 2)
    plt.plot(gc_positions, gc_content, 'b-', linewidth=2)
    plt.xlim(0, genome_length)
    plt.ylabel('GC含量')
    plt.gca().xaxis.set_major_formatter(
        FuncFormatter(lambda x, _: f'{x/1e6:.1f}Mb'))
    
    # 3. 绘制TE密度曲线
    plt.subplot(4, 1, 3)
    plt.plot(te_positions, te_density, 'g-', linewidth=2)
    plt.xlim(0, genome_length)
    plt.ylabel('TE密度 (%)')
    plt.gca().xaxis.set_major_formatter(
        FuncFormatter(lambda x, _: f'{x/1e6:.1f}Mb'))
    
    # 4. 如果有进化时间信息，绘制TE进化时间分布
    evol_times = [te.get('evolutionary_time', 0) for te in annotations]
    if any(evol_times):
        plt.subplot(4, 1, 4)
        
        # 创建散点图，X轴为位置，Y轴为进化时间，颜色表示相似度
        positions = [(te['start'] + te['end']) / 2 for te in annotations]
        similarities = [te.get('similarity', 0) for te in annotations]
        
        # 创建颜色映射：从红（低相似度）到绿（高相似度）
        cmap = LinearSegmentedColormap.from_list('similarity', ['red', 'yellow', 'green'])
        normalized_similarities = [s/100 for s in similarities]  # 规范化到0-1
        
        plt.scatter(positions, evol_times, c=normalized_similarities, 
                   cmap=cmap, alpha=0.7, s=30)
        
        plt.xlim(0, genome_length)
        plt.yscale('log')
        plt.ylabel('进化时间 (百万年)')
        plt.xlabel('染色体位置')
        plt.gca().xaxis.set_major_formatter(
            FuncFormatter(lambda x, _: f'{x/1e6:.1f}Mb'))
        
        # 添加颜色条
        cbar = plt.colorbar()
        cbar.set_label('序列相似度 (%)')
        cbar.set_ticks([0, 0.5, 1])
        cbar.set_ticklabels(['0%', '50%', '100%'])
    else:
        plt.subplot(4, 1, 4)
        # 绘制TE家族分布
        family_counts = Counter([te['family'] for te in annotations])
        top_families = dict(family_counts.most_common(10))
        
        positions = np.arange(len(top_families))
        plt.bar(positions, top_families.values(), align='center', alpha=0.7)
        plt.xticks(positions, top_families.keys(), rotation=45, ha='right')
        plt.ylabel('TE数量')
        plt.xlabel('TE家族')
        plt.tight_layout()
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_chromosome_structure.png", dpi=300)
    plt.close()
    
    # 创建GC含量与TE密度关系图
    plt.figure(figsize=(10, 8))
    
    # 将GC含量和TE密度数据对齐到相同的窗口
    aligned_gc = []
    aligned_te = []
    
    for gc_pos, gc_val in zip(gc_positions, gc_content):
        for te_pos, te_val in zip(te_positions, te_density):
            if abs(gc_pos - te_pos) < 5000:  # 如果位置接近
                aligned_gc.append(gc_val)
                aligned_te.append(te_val)
                break
    
    # 散点图
    plt.scatter(aligned_gc, aligned_te, alpha=0.6)
    
    # 拟合线
    if aligned_gc:
        z = np.polyfit(aligned_gc, aligned_te, 1)
        p = np.poly1d(z)
        plt.plot(sorted(aligned_gc), p(sorted(aligned_gc)), "r--", linewidth=2)
        
        # 计算相关系数
        correlation = np.corrcoef(aligned_gc, aligned_te)[0, 1]
        plt.annotate(f"相关系数: {correlation:.3f}", 
                    xy=(0.05, 0.95), xycoords='axes fraction')
    
    plt.xlabel('GC含量')
    plt.ylabel('TE密度 (%)')
    plt.title('GC含量与TE密度关系')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(f"{output_prefix}_gc_content_vs_te_density.png", dpi=300)
    plt.close()
    
    # 如果有进化时间信息，创建TE家族与进化时间关系图
    evol_times = [te.get('evolutionary_time', 0) for te in annotations]
    if any(evol_times):
        plt.figure(figsize=(12, 8))
        
        # 按家族分组进化时间
        family_times = defaultdict(list)
        for te in annotations:
            if te.get('evolutionary_time', 0) > 0:
                family_times[te['family']].append(te['evolutionary_time'])
        
        # 筛选至少有5个实例的家族
        filtered_families = {f: t for f, t in family_times.items() if len(t) >= 5}
        
        if filtered_families:
            # 创建箱线图
            plt.boxplot([times for times in filtered_families.values()], 
                        labels=filtered_families.keys(),
                        vert=False, showfliers=False)
            
            plt.xscale('log')
            plt.xlabel('进化时间 (百万年)')
            plt.ylabel('TE家族')
            plt.title('不同TE家族的进化时间分布')
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.tight_layout()
            plt.savefig(f"{output_prefix}_family_evolution_time.png", dpi=300)
            plt.close()

def main():
    parser = argparse.ArgumentParser(description='生成染色体结构与TE分布可视化图表')
    parser.add_argument('-g', '--genome', type=str, required=True,
                        help='基因组FASTA文件')
    parser.add_argument('-a', '--annotations', type=str, required=True,
                        help='TE注释文件 (BED或GFF格式)')
    parser.add_argument('-c', '--config', type=str, default=None,
                        help='基因组配置文件 (包含染色体区域信息)')
    parser.add_argument('-o', '--output_prefix', type=str, default='visualization',
                        help='输出文件前缀 (默认: visualization)')
    parser.add_argument('-f', '--format', type=str, choices=['png', 'pdf', 'svg'], default='png',
                        help='输出图像格式 (默认: png)')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.genome):
        sys.exit(f"错误: 基因组文件 {args.genome} 不存在")
    if not os.path.exists(args.annotations):
        sys.exit(f"错误: 注释文件 {args.annotations} 不存在")
    
    # 创建输出目录
    os.makedirs(os.path.dirname(os.path.abspath(args.output_prefix)), exist_ok=True)
    
    # 解析输入文件
    print("解析基因组和TE注释...")
    sequence = parse_fasta(args.genome)
    annotations = parse_annotations(args.annotations)
    
    # 加载配置文件（如果提供）
    config = None
    if args.config and os.path.exists(args.config):
        config = load_genome_config(args.config)
    
    print(f"基因组长度: {len(sequence):,} bp")
    print(f"TE注释数量: {len(annotations)}")
    
    # 生成可视化
    print("生成可视化图表...")
    visualize_chromosome_structure(sequence, annotations, config, args.output_prefix)
    
    print(f"完成! 可视化图表已保存:")
    print(f"染色体结构图: {args.output_prefix}_chromosome_structure.{args.format}")
    print(f"GC含量与TE密度关系图: {args.output_prefix}_gc_content_vs_te_density.{args.format}")
    print(f"家族进化时间图: {args.output_prefix}_family_evolution_time.{args.format}")

if __name__ == "__main__":
    main()

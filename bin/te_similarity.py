#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
te_similarity.py - 分析TE序列相似度和估算进化时间

此脚本分析插入的TE与原始TE序列的相似度，估算进化时间，并生成详细报告。
"""

import os
import sys
import argparse
import json
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from Bio import pairwise2
from Bio.Seq import Seq
import datetime
# 导入参数验证和并行处理模块
from parameter_validation import validate_similarity_params
from parallel_processing import calculate_similarity_metrics_parallel, reverse_complement, get_optimal_processes

def parse_fasta(fasta_file):
    """解析FASTA文件，返回序列字典"""
    sequences = {}
    current_id = None
    description = ""
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # 提取序列ID和描述
                parts = line[1:].split(None, 1)
                current_id = parts[0]
                description = parts[1] if len(parts) > 1 else ""
                sequences[current_id] = {'seq': [], 'desc': description}
            elif current_id:
                # 添加序列行
                sequences[current_id]['seq'].append(line)
    
    # 合并每个序列的所有行
    for seq_id in sequences:
        sequences[seq_id]['seq'] = ''.join(sequences[seq_id]['seq'])
    
    return sequences

def parse_te_annotations(annotation_file):
    """解析增强的BED格式TE注释文件"""
    annotations = []
    
    with open(annotation_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
                
            te_id = fields[0]
            family = fields[1]
            start = int(fields[2])
            end = int(fields[3])
            strand = fields[4]
            nested_in = fields[6] if len(fields) > 6 and fields[6] != 'None' else None
            
            # 尝试解析插入时间
            insertion_time = fields[9] if len(fields) > 9 else None
            
            annotations.append({
                'te_id': te_id,
                'family': family,
                'start': start,
                'end': end,
                'strand': strand,
                'nested_in': nested_in,
                'insertion_time': insertion_time
            })
    
    return annotations

def estimate_evolutionary_time(similarity, mutation_rate=1e-8):
    """
    基于序列相似度估算进化时间
    
    参数:
    similarity -- 序列相似度百分比
    mutation_rate -- 每碱基每代的突变率

    返回:
    time -- 估计的进化时间（单位：年）
    """
    # 相似度转换为差异度（0-1范围）
    divergence = (100 - similarity) / 100
    
    # Jukes-Cantor校正（处理多次突变的同一位点）
    if divergence < 0.75:  # JC69模型在差异度接近0.75时不适用
        corrected_divergence = -0.75 * np.log(1 - (4/3) * divergence)
    else:
        corrected_divergence = 1.0
    
    # 估算时间 = 校正后的差异度 / (2 * 突变率)
    # 假设平均每年20代，转换为年
    generations = corrected_divergence / (2 * mutation_rate)
    years = generations / 20
    
    return years

def generate_similarity_report(metrics, output_file, annotations=None):
    """生成相似度分析报告"""
    with open(output_file, 'w') as f:
        f.write("# TE序列相似度和进化时间分析报告\n")
        f.write(f"生成时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## 总体统计\n")
        similarities = [m['similarity'] for m in metrics.values()]
        f.write(f"分析的TE数量: {len(metrics)}\n")
        f.write(f"平均相似度: {np.mean(similarities):.2f}%\n")
        f.write(f"相似度范围: {min(similarities):.2f}% - {max(similarities):.2f}%\n\n")
        
        # 按家族分组统计
        family_metrics = defaultdict(list)
        for te_id, metric in metrics.items():
            family = metric['original_id']
            family_metrics[family].append(metric['similarity'])
        
        f.write("## 家族统计\n")
        f.write("家族\t数量\t平均相似度(%)\t平均进化时间(百万年)\n")
        for family, sims in sorted(family_metrics.items(), key=lambda x: -np.mean(x[1])):
            avg_sim = np.mean(sims)
            avg_time = estimate_evolutionary_time(avg_sim) / 1e6  # 转换为百万年
            f.write(f"{family}\t{len(sims)}\t{avg_sim:.2f}\t{avg_time:.2f}\n")
        
        f.write("\n## 嵌套TE分析\n")
        nested_tes = [te_id for te_id, metric in metrics.items() if metric['nested_in'] != 'None']
        f.write(f"嵌套TE数量: {len(nested_tes)}\n")
        
        if nested_tes:
            nested_sims = [metrics[te_id]['similarity'] for te_id in nested_tes]
            non_nested_sims = [metrics[te_id]['similarity'] for te_id in metrics if te_id not in nested_tes]
            
            f.write(f"嵌套TE平均相似度: {np.mean(nested_sims):.2f}%\n")
            f.write(f"非嵌套TE平均相似度: {np.mean(non_nested_sims):.2f}%\n\n")
        
        f.write("## 各TE详细信息\n")
        f.write("TE_ID\t家族\t相似度(%)\t身份值(%)\t原始长度\t插入长度\t嵌套于\t进化时间(年)\n")
        
        for te_id, metric in sorted(metrics.items(), key=lambda x: -x[1]['similarity']):
            evol_time = estimate_evolutionary_time(metric['similarity'])
            f.write(f"{te_id}\t{metric['original_id']}\t{metric['similarity']:.2f}\t"
                   f"{metric['identity']:.2f}\t{metric['original_length']}\t"
                   f"{metric['inserted_length']}\t{metric['nested_in']}\t"
                   f"{evol_time:.0f}\n")

def generate_visualization(metrics, output_prefix, annotations=None):
    """生成相似度和进化时间可视化"""
    # 提取数据
    similarities = [m['similarity'] for m in metrics.values()]
    evol_times = [estimate_evolutionary_time(m['similarity']) / 1e6 for m in metrics.values()]  # 百万年
    
    # 创建目录
    os.makedirs(os.path.dirname(os.path.abspath(output_prefix)), exist_ok=True)
    
    # 图1：相似度分布直方图
    plt.figure(figsize=(10, 6))
    plt.hist(similarities, bins=20, color='skyblue', edgecolor='black')
    plt.xlabel('相似度 (%)')
    plt.ylabel('TE数量')
    plt.title('TE序列相似度分布')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(f"{output_prefix}_similarity_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 图2：进化时间直方图
    plt.figure(figsize=(10, 6))
    plt.hist(evol_times, bins=20, color='lightgreen', edgecolor='black')
    plt.xlabel('进化时间 (百万年)')
    plt.ylabel('TE数量')
    plt.title('TE进化时间分布')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(f"{output_prefix}_evolutionary_time.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 图3：家族相似度箱线图
    family_similarities = defaultdict(list)
    for te_id, metric in metrics.items():
        family = metric['original_id']
        family_similarities[family].append(metric['similarity'])
    
    # 只保留至少有3个实例的家族，按中位数相似度排序
    filtered_families = {f: s for f, s in family_similarities.items() if len(s) >= 3}
    sorted_families = sorted(filtered_families.items(), key=lambda x: np.median(x[1]), reverse=True)
    
    if sorted_families:
        plt.figure(figsize=(12, 8))
        box_data = [sims for _, sims in sorted_families]
        labels = [family for family, _ in sorted_families]
        
        plt.boxplot(box_data, labels=labels, vert=False)
        plt.xlabel('相似度 (%)')
        plt.ylabel('TE家族')
        plt.title('不同TE家族的相似度分布')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_family_similarity.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # 图4：嵌套vs非嵌套相似度比较
    nested_tes = [te_id for te_id, metric in metrics.items() if metric['nested_in'] != 'None']
    
    if nested_tes:
        nested_sims = [metrics[te_id]['similarity'] for te_id in nested_tes]
        non_nested_sims = [metrics[te_id]['similarity'] for te_id in metrics if te_id not in nested_tes]
        
        plt.figure(figsize=(8, 6))
        plt.boxplot([non_nested_sims, nested_sims], labels=['非嵌套', '嵌套'])
        plt.ylabel('相似度 (%)')
        plt.title('嵌套vs非嵌套TE相似度比较')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.savefig(f"{output_prefix}_nested_comparison.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # 生成R脚本用于高级可视化
    with open(f"{output_prefix}_visualization.R", 'w') as f:
        f.write("""#!/usr/bin/env Rscript
# TE相似度和进化时间可视化脚本

# 加载必要的库
library(ggplot2)
library(reshape2)

# 读取数据
te_data <- read.table("%s_similarity_data.tsv", header=TRUE, sep="\\t")

# 绘制相似度-进化时间散点图
pdf("%s_similarity_vs_time.pdf", width=10, height=8)
ggplot(te_data, aes(x=Similarity, y=EvolutionaryTime, color=Family)) +
  geom_point(alpha=0.7) +
  theme_minimal() +
  labs(title="TE相似度与进化时间关系", 
       x="相似度 (%%)", y="进化时间 (百万年)") +
  scale_y_log10() +
  theme(legend.position="right")
dev.off()

# 绘制家族进化时间热图
if (length(unique(te_data$Family)) > 1) {
  family_time <- aggregate(EvolutionaryTime ~ Family, data=te_data, FUN=median)
  family_time <- family_time[order(family_time$EvolutionaryTime),]
  
  pdf("%s_family_time_heatmap.pdf", width=12, height=8)
  ggplot(family_time, aes(x=reorder(Family, -EvolutionaryTime), y=1, fill=EvolutionaryTime)) +
    geom_tile() +
    scale_fill_gradient(low="yellow", high="red") +
    theme_minimal() +
    labs(title="TE家族进化时间热图", x="TE家族", y="") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  dev.off()
}

print("可视化完成，生成了PDF文件")
""" % (output_prefix, output_prefix, output_prefix))
    
    # 生成数据文件供R脚本使用
    with open(f"{output_prefix}_similarity_data.tsv", 'w') as f:
        f.write("TE_ID\tFamily\tSimilarity\tIdentity\tOriginalLength\tInsertedLength\tNestedIn\tEvolutionaryTime\n")
        
        for te_id, metric in metrics.items():
            evol_time = estimate_evolutionary_time(metric['similarity']) / 1e6  # 百万年
            f.write(f"{te_id}\t{metric['original_id']}\t{metric['similarity']:.2f}\t"
                   f"{metric['identity']:.2f}\t{metric['original_length']}\t"
                   f"{metric['inserted_length']}\t{metric['nested_in']}\t"
                   f"{evol_time:.6f}\n")

def export_gff3_with_evolutionary_data(metrics, annotations, output_file):
    """导出带有进化数据的GFF3文件"""
    with open(output_file, 'w') as f:
        f.write("##gff-version 3\n")
        f.write(f"##date {datetime.datetime.now().strftime('%Y-%m-%d')}\n")
        f.write("##source TE_simulation_with_evolutionary_analysis\n")
        
        for anno in annotations:
            te_id = anno['te_id']
            
            if te_id in metrics:
                metric = metrics[te_id]
                similarity = metric['similarity']
                evol_time = estimate_evolutionary_time(similarity)
                evol_time_mya = evol_time / 1e6  # 百万年
                
                attributes = [
                    f"ID={te_id}",
                    f"Name={anno['family']}",
                    f"similarity={similarity:.2f}",
                    f"evolutionary_time={evol_time_mya:.2f}Mya"
                ]
                
                if anno['nested_in'] and anno['nested_in'] != 'None':
                    attributes.append(f"Parent={anno['nested_in']}")
                
                attributes_str = ';'.join(attributes)
                
                f.write(f"chr1\tTEsimulation\ttransposable_element\t{anno['start']+1}\t{anno['end']+1}\t.\t"
                       f"{anno['strand']}\t.\t{attributes_str}\n")

def main():
    parser = argparse.ArgumentParser(description='分析TE序列相似度并估算进化时间')
    parser.add_argument('-t', '--te_library', type=str, required=True,
                        help='原始TE文库FASTA文件')
    parser.add_argument('-s', '--te_sequences', type=str, required=True,
                        help='插入的TE序列FASTA文件')
    parser.add_argument('-a', '--annotations', type=str, default=None,
                        help='TE注释文件 (BED格式)')
    parser.add_argument('-o', '--output_prefix', type=str, default='te_similarity',
                        help='输出文件前缀 (默认: te_similarity)')
    parser.add_argument('-r', '--mutation_rate', type=float, default=1e-8,
                        help='每碱基每代的突变率 (默认: 1e-8)')
    parser.add_argument('--num_processes', type=int, default=None,
                        help='并行处理的进程数量 (默认: CPU核心数的75%)')
    
    args = parser.parse_args()
    
    # 参数验证
    validation_result = validate_similarity_params(args)
    validation_result.print_messages()
    validation_result.exit_if_invalid()
    
    # 设置并行处理进程数
    num_processes = args.num_processes if args.num_processes else get_optimal_processes()
    print(f"使用 {num_processes} 个进程进行并行计算")
    
    # 创建输出目录
    os.makedirs(os.path.dirname(os.path.abspath(args.output_prefix)), exist_ok=True)
    
    # 解析输入文件
    print("解析TE文库和插入的TE序列...")
    library_seqs = parse_fasta(args.te_library)
    inserted_seqs = parse_fasta(args.te_sequences)
    
    annotations = None
    if args.annotations and os.path.exists(args.annotations):
        print("解析TE注释...")
        annotations = parse_te_annotations(args.annotations)
    
    print(f"TE文库中TE数量: {len(library_seqs)}")
    print(f"插入的TE数量: {len(inserted_seqs)}")
    
    # 并行计算相似度指标
    print("计算序列相似度和进化指标...")
    similarity_metrics = calculate_similarity_metrics_parallel(
        library_seqs, inserted_seqs, num_processes)
    
    # 生成输出文件
    report_file = f"{args.output_prefix}_report.txt"
    gff_file = f"{args.output_prefix}_evolutionary.gff3"
    
    print("生成报告和可视化...")
    generate_similarity_report(similarity_metrics, report_file, annotations)
    generate_visualization(similarity_metrics, args.output_prefix, annotations)
    
    if annotations:
        export_gff3_with_evolutionary_data(similarity_metrics, annotations, gff_file)
    
    print(f"完成! 分析结果已保存:")
    print(f"相似度报告: {report_file}")
    print(f"可视化文件: {args.output_prefix}_*.png")
    print(f"R可视化脚本: {args.output_prefix}_visualization.R")
    if annotations:
        print(f"带进化数据的GFF3: {gff_file}")
    
    print("\n要生成高级可视化图表，请运行: Rscript %s_visualization.R" % args.output_prefix)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_annotation.py - 生成最终注释和统计报告

此脚本分析最终基因组和TE注释文件，生成详细统计报告和可视化友好的GFF3格式注释。
用于TE扩增模拟工作流的最终步骤。
"""

import os
import sys
import argparse
import json
import re
from collections import defaultdict, Counter
import datetime
import matplotlib.pyplot as plt
import numpy as np
# 导入并行处理功能
from parallel_processing import get_optimal_processes

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
            if len(fields) < 6:
                continue
                
            te_id = fields[0]
            family = fields[1]
            start = int(fields[2])
            end = int(fields[3])
            strand = fields[4]
            length = int(fields[5])
            nested_in = fields[6] if len(fields) > 6 and fields[6] != 'None' else None
            mutations = fields[7] if len(fields) > 7 and fields[7] != 'None' else None
            
            annotations.append({
                'te_id': te_id,
                'family': family,
                'start': start,
                'end': end,
                'strand': strand,
                'length': length,
                'nested_in': nested_in,
                'mutations': mutations
            })
    
    return annotations

def analyze_mutations(mutations_str):
    """分析突变字符串，提取SNP和SV信息"""
    if not mutations_str or mutations_str == 'None':
        return {'snp_count': 0, 'sv_count': 0, 'sv_types': {}}
    
    parts = mutations_str.split(';')
    snp_count = sum(1 for p in parts if p.startswith('SNP:'))
    sv_count = sum(1 for p in parts if p.startswith('SV:'))
    
    sv_types = Counter()
    for part in parts:
        if part.startswith('SV:'):
            # 提取SV类型信息
            sv_fields = part.split(':')
            if len(sv_fields) >= 3:
                sv_type = sv_fields[2]
                # 去除可能的详细信息
                if ',' in sv_type:
                    sv_type = sv_type.split(',')[0]
                sv_types[sv_type] += 1
    
    return {
        'snp_count': snp_count,
        'sv_count': sv_count,
        'sv_types': dict(sv_types)
    }

def calculate_gc_content(sequence):
    """计算序列的GC含量，支持字符串和字节对象"""
    if not sequence:
        return 0
    
    # 检查序列类型并相应处理
    if isinstance(sequence, (bytes, bytearray)):
        g_count = sequence.count(ord('G')) + sequence.count(ord('g'))
        c_count = sequence.count(ord('C')) + sequence.count(ord('c'))
        return (g_count + c_count) / len(sequence) * 100
    else:
        # 原始字符串处理方式
        g_count = sequence.upper().count('G')
        c_count = sequence.upper().count('C')
        return (g_count + c_count) / len(sequence) * 100

# 这个函数从analyze_annotations内部移到全局作用域
def process_annotation_batch(annotations_batch, genome_bytes):
    """批量处理注释以提高效率"""
    result = {
        'te_sequences': [],
        'positions': set(),
        'total_bases': 0,
        'strand_counts': {'forward': 0, 'reverse': 0},
        'lengths': [],
        'family_counts': Counter(),
        'nested_count': 0,
        'mutation_stats': {'snp_count': 0, 'sv_count': 0, 'sv_types': Counter()},
        'te_positions': []
    }
    
    for anno in annotations_batch:
        # 处理一批注释...
        result['total_bases'] += anno['length']
        result['strand_counts']['forward' if anno['strand'] == '+' else 'reverse'] += 1
        result['lengths'].append(anno['length'])
        result['family_counts'][anno['family']] += 1
        
        if anno['nested_in']:
            result['nested_count'] += 1
        
        # 突变统计
        if anno.get('mutations'):
            mutation_info = analyze_mutations(anno['mutations'])
            result['mutation_stats']['snp_count'] += mutation_info['snp_count']
            result['mutation_stats']['sv_count'] += mutation_info['sv_count']
            for sv_type, count in mutation_info['sv_types'].items():
                result['mutation_stats']['sv_types'][sv_type] += count
        
        # 提取TE序列
        if anno['start'] >= 0 and anno['end'] < len(genome_bytes):
            # 使用切片而不是索引循环提高效率
            te_seq = genome_bytes[anno['start']:anno['end']+1]
            result['te_sequences'].append(bytes(te_seq))
            
            # 添加TE覆盖的位置
            for pos in range(anno['start'], anno['end']+1):
                result['positions'].add(pos)
        
        # 记录TE位置
        result['te_positions'].append((anno['start'], anno['end']))
    
    return result

# 改进内存使用和处理效率
def analyze_annotations(annotations, genome_seq):
    """分析TE注释，生成统计信息 - 优化版"""
    # 使用bytearray储存序列可以减少内存使用
    if isinstance(genome_seq, str):
        genome_bytes = bytearray(genome_seq.encode('ascii'))
    else:
        genome_bytes = genome_seq
        
    genome_length = len(genome_bytes)
    
    # 预分配内存跟踪TE位置
    all_te_positions = set()
    # 预计算数组大小可加快执行速度
    all_te_positions_capacity = sum(a['end'] - a['start'] + 1 for a in annotations)
    if hasattr(set, "reserve"):  # Python 3.11+
        all_te_positions.reserve(all_te_positions_capacity)
    
    # 预分配统计集合和数组
    stats = {
        'total_te_count': len(annotations),
        'total_te_bases': 0,
        'genome_length': genome_length,
        'te_percentage': 0,
        'strand_distribution': {'forward': 0, 'reverse': 0},
        'length_stats': {'min': float('inf'), 'max': 0, 'avg': 0, 'median': 0},
        'family_counts': Counter(),
        'nested_te_count': 0,
        'mutation_stats': {'snp_count': 0, 'sv_count': 0, 'sv_types': Counter()},
        'gc_content': {'genome': calculate_gc_content(genome_bytes), 'te': 0, 'non_te': 0},
        'te_positions': []
    }
    
    # 多线程处理注释详情
    import multiprocessing as mp
    from functools import partial
    
    # 划分任务进行并行处理
    num_processes = get_optimal_processes(task_type='cpu_bound')
    batch_size = max(1, len(annotations) // (num_processes * 2))
    annotation_batches = [annotations[i:i+batch_size] for i in range(0, len(annotations), batch_size)]
    
    # 使用全局函数而不是局部函数，以避免pickle错误
    process_func = partial(process_annotation_batch, genome_bytes=genome_bytes)
    
    # 并行处理注释
    with mp.Pool(processes=num_processes) as pool:
        batch_results = pool.map(process_func, annotation_batches)
    
    # 合并批处理结果
    lengths = []
    te_sequences = []
    
    for result in batch_results:
        stats['total_te_bases'] += result['total_bases']
        stats['strand_distribution']['forward'] += result['strand_counts']['forward']
        stats['strand_distribution']['reverse'] += result['strand_counts']['reverse']
        lengths.extend(result['lengths'])
        stats['family_counts'].update(result['family_counts'])
        stats['nested_te_count'] += result['nested_count']
        stats['mutation_stats']['snp_count'] += result['mutation_stats']['snp_count']
        stats['mutation_stats']['sv_count'] += result['mutation_stats']['sv_count']
        stats['mutation_stats']['sv_types'].update(result['mutation_stats']['sv_types'])
        stats['te_positions'].extend(result['te_positions'])
        te_sequences.extend(result['te_sequences'])
        all_te_positions.update(result['positions'])
    
    # 计算长度统计
    if lengths:
        stats['length_stats']['min'] = min(lengths)
        stats['length_stats']['max'] = max(lengths)
        stats['length_stats']['avg'] = sum(lengths) / len(lengths)
        stats['length_stats']['median'] = sorted(lengths)[len(lengths) // 2]
    
    # 计算TE占比
    stats['te_percentage'] = (stats['total_te_bases'] / stats['genome_length']) * 100
    
    # 计算TE和非TE区域的GC含量
    if te_sequences:
        # 连接所有TE序列计算GC含量
        all_te_seq = b''.join(te_sequences)
        stats['gc_content']['te'] = (all_te_seq.count(b'G') + all_te_seq.count(b'C')) / len(all_te_seq) * 100
    
    # 高效地计算非TE区域的GC含量
    if len(all_te_positions) < genome_length:
        non_te_positions = [i for i in range(genome_length) if i not in all_te_positions]
        if len(non_te_positions) > 10000:  # 只取样本计算GC含量
            import random
            sample_size = 10000
            non_te_sample = random.sample(non_te_positions, sample_size)
            non_te_seq = bytearray(genome_bytes[i] for i in sorted(non_te_sample))
        else:
            non_te_seq = bytearray(genome_bytes[i] for i in sorted(non_te_positions))
        
        if non_te_seq:
            stats['gc_content']['non_te'] = (non_te_seq.count(b'G') + non_te_seq.count(b'C')) / len(non_te_seq) * 100
    
    return stats

def generate_statistics_report(stats, output_file):
    """生成统计报告"""
    with open(output_file, 'w') as f:
        f.write("# TE扩增模拟统计报告\n")
        f.write(f"生成时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## 基本统计\n")
        f.write(f"基因组长度: {stats['genome_length']:,} bp\n")
        f.write(f"TE总数: {stats['total_te_count']:,}\n")
        f.write(f"TE总碱基数: {stats['total_te_bases']:,} bp\n")
        f.write(f"TE占比: {stats['te_percentage']:.2f}%\n\n")
        
        f.write("## GC含量统计\n")
        f.write(f"全基因组GC含量: {stats['gc_content']['genome']:.2f}%\n")
        f.write(f"TE区域GC含量: {stats['gc_content']['te']:.2f}%\n")
        f.write(f"非TE区域GC含量: {stats['gc_content']['non_te']:.2f}%\n\n")
        
        f.write("## TE长度统计\n")
        f.write(f"最小长度: {stats['length_stats']['min']:,} bp\n")
        f.write(f"最大长度: {stats['length_stats']['max']:,} bp\n")
        f.write(f"平均长度: {stats['length_stats']['avg']:.2f} bp\n")
        f.write(f"中位长度: {stats['length_stats']['median']:,} bp\n\n")
        
        f.write("## 链分布\n")
        f.write(f"正链: {stats['strand_distribution']['forward']} ({stats['strand_distribution']['forward']/stats['total_te_count']*100:.2f}%)\n")
        f.write(f"负链: {stats['strand_distribution']['reverse']} ({stats['strand_distribution']['reverse']/stats['total_te_count']*100:.2f}%)\n\n")
        
        f.write("## 嵌套TE统计\n")
        f.write(f"嵌套TE数量: {stats['nested_te_count']} ({stats['nested_te_count']/stats['total_te_count']*100:.2f}%)\n\n")
        
        f.write("## 突变统计\n")
        f.write(f"SNP总数: {stats['mutation_stats']['snp_count']:,}\n")
        f.write(f"SV总数: {stats['mutation_stats']['sv_count']:,}\n")
        
        if stats['mutation_stats']['sv_types']:
            f.write("SV类型分布:\n")
            for sv_type, count in sorted(stats['mutation_stats']['sv_types'].items(), key=lambda x: -x[1]):
                f.write(f"  - {sv_type}: {count} ({count/stats['mutation_stats']['sv_count']*100:.2f}%)\n")
        f.write("\n")
        
        f.write("## TE家族分布\n")
        for family, count in stats['family_counts'].most_common():
            f.write(f"{family}: {count} ({count/stats['total_te_count']*100:.2f}%)\n")

def generate_gff3(annotations, output_file):
    """生成增强的GFF3格式注释文件，包含丰富的TE信息"""
    with open(output_file, 'w') as f:
        f.write("##gff-version 3\n")
        f.write(f"##date {datetime.datetime.now().strftime('%Y-%m-%d')}\n")
        f.write("##source TE_simulation_pipeline\n")
        
        for anno in annotations:
            # 提取突变信息
            mutation_info = analyze_mutations(anno.get('mutations', 'None'))
            snp_count = mutation_info['snp_count']
            sv_count = mutation_info['sv_count']
            
            # 确定完整性状态
            completeness = "complete"
            if sv_count > 0 and any('deletion' in sv_type for sv_type in mutation_info['sv_types']):
                completeness = "fragmented"
            
            # 构建属性字段
            attributes = [
                f"ID={anno['te_id']}",
                f"Name={anno['family']}",
                f"length={anno['length']}",
                f"snp_count={snp_count}",
                f"sv_count={sv_count}",
                f"completeness={completeness}"
            ]
            
            # 添加嵌套信息
            if anno.get('nested_in') and anno['nested_in'] != 'None':
                attributes.append(f"Parent={anno['nested_in']}")
            
            attributes_str = ';'.join(attributes)
            
            # 写入GFF3行
            f.write(f"chr1\tTE_simulation\ttransposable_element\t{anno['start']+1}\t{anno['end']+1}\t.\t"
                  f"{anno['strand']}\t.\t{attributes_str}\n")

def generate_te_distribution_plot(stats, output_file):
    """生成TE分布图"""
    plt.figure(figsize=(12, 6))
    
    # 创建基因组覆盖率图
    genome_length = stats['genome_length']
    bin_size = max(1000, genome_length // 1000)  # 自适应bin大小
    bins = np.zeros(genome_length // bin_size + 1)
    
    # 统计每个bin的TE覆盖
    for start, end in stats['te_positions']:
        start_bin = start // bin_size
        end_bin = end // bin_size
        
        for bin_idx in range(start_bin, end_bin + 1):
            if bin_idx < len(bins):
                bins[bin_idx] += 1
    
    # 绘制分布图
    bin_positions = np.arange(len(bins)) * bin_size
    plt.plot(bin_positions, bins, 'b-', linewidth=2)
    
    plt.title('基因组TE分布')
    plt.xlabel('基因组位置 (bp)')
    plt.ylabel('TE数量')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # 格式化x轴标签
    plt.gca().xaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f'{x/1e6:.1f}Mb'))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def generate_te_family_chart(stats, output_file):
    """生成TE家族分布饼图"""
    plt.figure(figsize=(10, 8))
    
    # 获取前10个家族，其余归为"其他"
    top_families = stats['family_counts'].most_common(10)
    labels = [family for family, _ in top_families]
    sizes = [count for _, count in top_families]
    
    # 如果有其他家族，添加"其他"类别
    others_count = sum(count for family, count in stats['family_counts'].items() 
                     if family not in labels)
    if others_count > 0:
        labels.append('其他')
        sizes.append(others_count)
    
    # 创建饼图
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, shadow=True)
    plt.axis('equal')  # 确保饼图是圆的
    plt.title('TE家族分布')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def generate_te_length_histogram(stats, annotations, output_file):
    """生成TE长度分布直方图"""
    plt.figure(figsize=(10, 6))
    
    # 提取TE长度
    lengths = [anno['length'] for anno in annotations]
    
    # 创建直方图
    bins = np.logspace(np.log10(min(lengths)), np.log10(max(lengths)), 20)
    plt.hist(lengths, bins=bins, alpha=0.7, color='skyblue', edgecolor='black')
    
    plt.title('TE长度分布')
    plt.xlabel('长度 (bp) [对数尺度]')
    plt.ylabel('TE数量')
    plt.xscale('log')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def generate_r_visualization_script(output_prefix):
    """生成R语言可视化脚本"""
    with open(f"{output_prefix}_visualization.R", 'w') as f:
        f.write("""#!/usr/bin/env Rscript
# TE扩增模拟可视化脚本

# 加载必要的库
library(ggplot2)
library(reshape2)
library(dplyr)

# 读取数据
te_anno <- read.delim("%s_te_annotations.tsv", header=TRUE, sep="\\t")

# 创建TE家族分布条形图
pdf("%s_family_distribution.pdf", width=10, height=8)
family_counts <- table(te_anno$Family)
family_df <- data.frame(Family=names(family_counts), Count=as.numeric(family_counts))
family_df <- family_df[order(-family_df$Count),]

# 只展示前15个家族，其余归为"其他"
if(nrow(family_df) > 15) {
  top_families <- family_df[1:15,]
  other_count <- sum(family_df[16:nrow(family_df), "Count"])
  top_families <- rbind(top_families, data.frame(Family="Other", Count=other_count))
  family_df <- top_families
}

ggplot(family_df, aes(x=reorder(Family, -Count), y=Count, fill=Family)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(title="TE家族分布", x="TE家族", y="数量") +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none")
dev.off()

# 创建TE长度分布箱线图
if("Length" %%in%% colnames(te_anno)) {
  pdf("%s_length_distribution.pdf", width=10, height=8)
  # 按家族分组
  top_families <- names(sort(table(te_anno$Family), decreasing=TRUE)[1:min(10, length(unique(te_anno$Family)))])
  family_data <- te_anno[te_anno$Family %%in%% top_families,]
  
  ggplot(family_data, aes(x=reorder(Family, -Length, median), y=Length, fill=Family)) +
    geom_boxplot() +
    scale_y_log10() +
    theme_minimal() +
    labs(title="TE长度分布(按家族)", x="TE家族", y="长度(bp, 对数尺度)") +
    theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none")
  dev.off()
}

# 创建TE位置分布图
if(all(c("Start", "End") %%in%% colnames(te_anno))) {
  pdf("%s_position_distribution.pdf", width=12, height=6)
  # 计算TE密度
  te_anno$Position <- (te_anno$Start + te_anno$End) / 2
  ggplot(te_anno, aes(x=Position)) +
    geom_density(fill="skyblue", alpha=0.7) +
    theme_minimal() +
    labs(title="TE在基因组中的分布", x="基因组位置", y="密度") +
    scale_x_continuous(labels=function(x) paste0(x/1e6, "Mb"))
  dev.off()
}

print("可视化完成，生成了PDF文件")
""" % (output_prefix, output_prefix, output_prefix, output_prefix))

def export_annotations_tsv(annotations, output_file):
    """将注释导出为TSV格式，用于R脚本处理"""
    with open(output_file, 'w') as f:
        # 写入表头
        headers = ["TE_ID", "Family", "Start", "End", "Strand", "Length", "Nested_In"]
        f.write('\t'.join(headers) + '\n')
        
        # 写入数据
        for anno in annotations:
            nested = str(anno.get('nested_in', 'None'))
            if nested == 'None':
                nested = 'NA'
                
            row = [
                anno['te_id'],
                anno['family'],
                str(anno['start']),
                str(anno['end']),
                anno['strand'],
                str(anno['length']),
                nested
            ]
            f.write('\t'.join(row) + '\n')

def main():
    parser = argparse.ArgumentParser(description='生成最终注释和统计报告')
    parser.add_argument('-g', '--genome', type=str, required=True,
                        help='输入基因组FASTA文件')
    parser.add_argument('-a', '--annotations', type=str, required=True,
                        help='TE注释文件 (BED格式)')
    parser.add_argument('-o', '--output_prefix', type=str, default='te_annotation',
                        help='输出文件前缀 (默认: te_annotation)')
    parser.add_argument('--num_processes', type=int, default=None,
                        help='并行处理的进程数量 (默认: CPU核心数的75%)')
    
    args = parser.parse_args()
    
    # 设置并行处理进程数
    num_processes = args.num_processes if args.num_processes else get_optimal_processes()
    print(f"使用 {num_processes} 个进程进行并行计算")
    
    # 解析输入文件
    print("解析基因组和TE注释...")
    genome_dict = parse_fasta(args.genome)
    annotations = parse_te_annotations(args.annotations)
    
    if not genome_dict:
        sys.exit("错误: 基因组文件为空或格式不正确")
    if not annotations:
        sys.exit("错误: TE注释为空或格式不正确")
    
    # 获取第一条染色体序列
    chromosome_id = list(genome_dict.keys())[0]
    genome_sequence = genome_dict[chromosome_id]['seq']
    
    print(f"基因组长度: {len(genome_sequence):,} bp")
    print(f"TE注释数量: {len(annotations)}")
    
    # 创建输出目录
    os.makedirs(os.path.dirname(os.path.abspath(args.output_prefix)), exist_ok=True)
    
    # 分析注释
    print("分析TE注释并生成统计信息...")
    stats = analyze_annotations(annotations, genome_sequence)
    
    # 生成输出文件
    print("生成输出文件...")
    
    # 统计报告
    stats_file = f"{args.output_prefix}_te_statistics.txt"
    generate_statistics_report(stats, stats_file)
    
    # GFF3文件
    gff_file = f"{args.output_prefix}_te_annotations.gff3"
    generate_gff3(annotations, gff_file)
    
    # TSV文件
    tsv_file = f"{args.output_prefix}_te_annotations.tsv"
    export_annotations_tsv(annotations, tsv_file)
    
    # 可视化图表
    te_distribution_file = f"{args.output_prefix}_te_distribution.png"
    generate_te_distribution_plot(stats, te_distribution_file)
    
    te_family_file = f"{args.output_prefix}_te_family_distribution.png"
    generate_te_family_chart(stats, te_family_file)
    
    te_length_file = f"{args.output_prefix}_te_length_distribution.png"
    generate_te_length_histogram(stats, annotations, te_length_file)
    
    # R脚本
    generate_r_visualization_script(args.output_prefix)
    
    print(f"完成! 输出文件已保存:")
    print(f"统计报告: {stats_file}")
    print(f"GFF3注释: {gff_file}")
    print(f"TSV注释: {tsv_file}")
    print(f"TE分布图: {te_distribution_file}")
    print(f"TE家族分布图: {te_family_file}")
    print(f"TE长度分布图: {te_length_file}")
    print(f"R可视化脚本: {args.output_prefix}_visualization.R")
    
    print("\n要生成高级可视化图表，请运行: Rscript %s_visualization.R" % args.output_prefix)

if __name__ == "__main__":
    main()

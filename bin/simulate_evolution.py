#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
simulate_evolution.py - 模拟基因组进化过程

此脚本根据指定的进化水平引入SNP和结构变异(SV)，包括缺失、倒位和重复。
增强版本支持TE序列跟踪和进化时间估算。
"""

import random
import argparse
import os
import sys
import re
from collections import defaultdict
import datetime
import json
import copy
# 导入参数验证和并行处理模块
from parameter_validation import validate_evolution_params
from parallel_processing import introduce_snp_parallel, get_optimal_processes

# 定义进化水平参数
EVOLUTION_LEVELS = {
    'low': {'snp': 0.01, 'sv': 0.005},   # 低: SNP 1%, SV 0.5%
    'medium': {'snp': 0.03, 'sv': 0.01}, # 中: SNP 3%, SV 1%
    'high': {'snp': 0.05, 'sv': 0.03}    # 高: SNP 5%, SV 3%
}

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

def parse_annotations(annotation_file):
    """解析BED格式的TE注释文件"""
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
            nested_in = fields[6] if len(fields) > 6 and fields[6] != 'None' else None
            mutations = fields[7] if len(fields) > 7 and fields[7] != 'None' else None
            
            annotations.append({
                'te_id': te_id,
                'family': family,
                'start': start,
                'end': end,
                'strand': strand,
                'nested_in': nested_in,
                'mutations': mutations,
                'new_mutations': []  # 将记录新突变
            })
    
    return annotations

def introduce_sv(sequence, rate, te_regions=None, te_multiplier=1.5, base_rate=1e-8):
    """
    引入结构变异(SV)，包括缺失、倒位和重复
    
    参数:
    sequence -- 基因组序列
    rate -- 基础SV比率
    te_regions -- TE区域列表 [(start, end), ...]
    te_multiplier -- TE区域的SV比率乘数
    base_rate -- 基础突变率，用于调整最终突变概率
    
    返回:
    modified_sequence -- 引入SV后的序列
    sv_events -- SV事件列表 [(type, start, end, details), ...]
    """
    seq_list = list(sequence)
    sv_events = []
    
    # 调整突变率，使其与基础突变率关联
    adjusted_rate = rate * (base_rate * 1e8)  # 调整比例，保持与原始相似
    
    # SV类型分布
    sv_types = {
        'small_deletion': 0.4,   # 小型缺失 (1-10bp)
        'medium_deletion': 0.2,  # 中型缺失 (11-100bp)
        'large_deletion': 0.1,   # 大型缺失 (101-1000bp)
        'inversion': 0.15,       # 倒位 (100-1000bp)
        'tandem_repeat': 0.15    # 串联重复 (10-100bp)
    }
    
    # SV大小范围
    sv_sizes = {
        'small_deletion': (1, 10),
        'medium_deletion': (11, 100),
        'large_deletion': (101, 1000),
        'inversion': (100, 1000),
        'tandem_repeat': (10, 100)
    }
    
    # 如果有TE区域，创建映射以快速查找位置是否在TE内
    te_map = {}
    if te_regions:
        for start, end in te_regions:
            for pos in range(start, end + 1):
                te_map[pos] = True
    
    # 计算要引入的SV数量
    total_sv = int(len(sequence) * adjusted_rate)
    
    # 为每种SV类型分配数量
    sv_counts = {}
    remaining = total_sv
    for sv_type, proportion in sv_types.items():
        if sv_type == list(sv_types.keys())[-1]:  # 最后一种类型
            sv_counts[sv_type] = remaining
        else:
            count = int(total_sv * proportion)
            sv_counts[sv_type] = count
            remaining -= count
    
    # 引入各种SV
    offset = 0  # 跟踪序列长度变化
    
    # 收集所有将要引入的SV
    all_svs = []
    for sv_type, count in sv_counts.items():
        for _ in range(count):
            # 随机选择位置
            pos = random.randint(0, len(sequence) - 1)
            
            # 根据位置调整SV率
            current_rate = adjusted_rate
            if te_regions and pos in te_map:
                current_rate *= te_multiplier
            
            # 只有通过随机检查才引入SV
            if random.random() > current_rate:
                continue
                
            # 确定SV大小
            min_size, max_size = sv_sizes[sv_type]
            size = random.randint(min_size, max_size)
            
            # 确保不超出序列范围
            end_pos = min(pos + size - 1, len(sequence) - 1)
            
            all_svs.append((sv_type, pos, end_pos))
    
    # 按位置排序SV，从后向前应用（避免位置冲突）
    all_svs.sort(key=lambda x: x[1], reverse=True)
    
    # 应用SV
    for sv_type, start, end in all_svs:
        size = end - start + 1
        
        if sv_type.endswith('deletion'):
            # 缺失: 移除指定区域
            deleted_seq = ''.join(seq_list[start:end+1])
            del seq_list[start:end+1]
            sv_events.append(('deletion', start, end, {'size': size, 'deleted_sequence': deleted_seq}))
            offset -= size
            
        elif sv_type == 'inversion':
            # 倒位: 反转指定区域
            inverted_seq = seq_list[start:end+1]
            inverted_seq.reverse()
            seq_list[start:end+1] = inverted_seq
            sv_events.append(('inversion', start, end, {'size': size}))
            
        elif sv_type == 'tandem_repeat':
            # 串联重复: 复制指定区域并插入
            repeat_seq = seq_list[start:end+1]
            repeat_times = random.randint(1, 3)  # 重复1-3次
            for _ in range(repeat_times):
                seq_list[end+1:end+1] = repeat_seq
            sv_events.append(('tandem_repeat', start, end, 
                             {'size': size, 'repeat_times': repeat_times, 
                              'new_end': end + size * repeat_times}))
            offset += size * repeat_times
    
    return ''.join(seq_list), sv_events

def update_annotations(annotations, snp_positions, sv_events):
    """
    根据SNP和SV事件更新TE注释
    
    参数:
    annotations -- 原始TE注释列表
    snp_positions -- SNP位置列表 [(position, old_base, new_base), ...]
    sv_events -- SV事件列表 [(type, start, end, details), ...]
    
    返回:
    updated_annotations -- 更新后的TE注释
    """
    # 创建SNP位置映射，用于快速查找
    snp_map = {pos: (old, new) for pos, old, new in snp_positions}
    
    # 按照位置排序SV事件（从前到后）
    sorted_sv = sorted(sv_events, key=lambda x: x[1])
    
    # 获取当前时间作为"变异时间"
    mutation_timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # 更新每个TE的注释
    updated_annotations = []
    
    for te in annotations:
        te_start = te['start']
        te_end = te['end']
        te_new_mutations = []
        
        # 解析现有的突变信息
        if te['mutations'] and te['mutations'] != 'None':
            te['new_mutations'] = []  # 确保新突变列表存在
        
        # 检查该TE区域的SNP
        for pos in range(te_start, te_end + 1):
            if pos in snp_map:
                old_base, new_base = snp_map[pos]
                relative_pos = pos - te_start
                te_new_mutations.append(('SNP', relative_pos, old_base, new_base))
        
        # 检查该TE是否受SV影响
        for sv_type, sv_start, sv_end, details in sorted_sv:
            # 检查SV是否与TE重叠
            if not (sv_end < te_start or sv_start > te_end):
                # SV与TE有重叠
                if sv_type == 'deletion':
                    # 如果TE完全被删除
                    if sv_start <= te_start and sv_end >= te_end:
                        te_new_mutations.append(('SV', 0, sv_type, 'complete_deletion'))
                        te_start = -1  # 标记为已删除
                        te_end = -1
                        break
                    
                    # 如果TE部分被删除
                    overlap_start = max(sv_start, te_start)
                    overlap_end = min(sv_end, te_end)
                    relative_start = overlap_start - te_start
                    relative_end = overlap_end - te_start
                    te_new_mutations.append(('SV', relative_start, sv_type, 
                                        {'start': relative_start, 'end': relative_end,
                                         'size': overlap_end - overlap_start + 1}))
                    
                elif sv_type == 'inversion':
                    # TE部分或完全倒位
                    overlap_start = max(sv_start, te_start)
                    overlap_end = min(sv_end, te_end)
                    relative_start = overlap_start - te_start
                    relative_end = overlap_end - te_start
                    te_new_mutations.append(('SV', relative_start, sv_type, 
                                        {'start': relative_start, 'end': relative_end,
                                         'size': overlap_end - overlap_start + 1}))
                    
                elif sv_type == 'tandem_repeat':
                    # TE部分被重复
                    if te_start >= sv_start and te_end <= sv_end:
                        repeat_times = details['repeat_times']
                        te_new_mutations.append(('SV', 0, sv_type, 
                                           {'repeat_times': repeat_times, 'complete': True}))
                    else:
                        overlap_start = max(sv_start, te_start)
                        overlap_end = min(sv_end, te_end)
                        relative_start = overlap_start - te_start
                        relative_end = overlap_end - te_start
                        te_new_mutations.append(('SV', relative_start, sv_type, 
                                            {'start': relative_start, 'end': relative_end,
                                             'repeat_times': details['repeat_times']}))
        
        # 如果TE未被完全删除，则更新注释
        if te_start != -1:
            te['new_mutations'] = te_new_mutations
            updated_annotations.append(te)
    
    return updated_annotations

def update_te_sequences(te_sequences, annotations, sequence, snp_positions, sv_events):
    """
    根据SNP和SV事件更新TE序列
    
    参数:
    te_sequences -- TE序列字典 {id: {'seq': sequence, 'desc': description}}
    annotations -- TE注释列表
    sequence -- 更新后的基因组序列
    snp_positions -- SNP位置列表
    sv_events -- SV事件列表
    
    返回:
    updated_sequences -- 更新后的TE序列字典
    """
    updated_sequences = {}
    
    for te in annotations:
        te_id = te['te_id']
        start = te['start']
        end = te['end']
        
        # 如果TE被完全删除，则跳过
        if start == -1 or end == -1:
            continue
        
        # 从基因组中提取更新后的TE序列
        te_seq = sequence[start:end+1]
        
        # 创建或更新描述
        desc_parts = []
        if te_id in te_sequences:
            # 提取原始描述中的信息
            orig_desc = te_sequences[te_id]['desc']
            for part in orig_desc.split():
                if not any(p in part for p in ['start=', 'end=', 'strand=', 'nested_in=']):
                    desc_parts.append(part)
        
        # 添加位置和嵌套信息
        desc_parts.append(f"start={start}")
        desc_parts.append(f"end={end}")
        desc_parts.append(f"strand={te['strand']}")
        desc_parts.append(f"nested_in={te['nested_in'] if te['nested_in'] else 'None'}")
        
        # 更新序列和描述
        updated_sequences[te_id] = {
            'seq': te_seq,
            'desc': ' '.join(desc_parts)
        }
    
    return updated_sequences

def write_evolved_annotations(annotations, output_file):
    """将更新后的TE注释写入文件"""
    with open(output_file, 'w') as f:
        f.write("# TE_ID\tFamily\tStart\tEnd\tStrand\tLength\tNested_In\tMutations\n")
        
        for te in annotations:
            # 将新突变与原有突变合并
            all_mutations = []
            
            # 处理原有突变
            if te.get('mutations') and te['mutations'] != 'None':
                all_mutations.extend(te['mutations'].split(';'))
            
            # 处理新突变
            for m in te.get('new_mutations', []):
                if m[0] == 'SNP':
                    all_mutations.append(f"SNP:{m[1]}:{m[2]}>{m[3]}")
                elif m[0] == 'SV':
                    if isinstance(m[3], str):
                        all_mutations.append(f"SV:{m[1]}:{m[2]}:{m[3]}")
                    else:
                        details = ','.join([f"{k}={v}" for k, v in m[3].items()])
                        all_mutations.append(f"SV:{m[1]}:{m[2]}:{details}")
            
            mutations_str = ';'.join(all_mutations) if all_mutations else "None"
            nested_str = te['nested_in'] if te['nested_in'] else "None"
            length = te['end'] - te['start'] + 1
            
            f.write(f"{te['te_id']}\t{te['family']}\t{te['start']}\t{te['end']}\t"
                   f"{te['strand']}\t{length}\t{nested_str}\t{mutations_str}\n")

def write_evolved_gff3(annotations, output_file):
    """将更新后的TE注释写入GFF3格式文件"""
    with open(output_file, 'w') as f:
        f.write("##gff-version 3\n")
        
        for te in annotations:
            # 合并所有突变信息
            all_mutations = []
            
            # 处理原有突变
            if te.get('mutations') and te['mutations'] != 'None':
                all_mutations.extend(te['mutations'].split(';'))
            
            # 处理新突变
            for m in te.get('new_mutations', []):
                if m[0] == 'SNP':
                    all_mutations.append(f"SNP_{m[1]}_{m[2]}to{m[3]}")
                elif m[0] == 'SV':
                    if isinstance(m[3], str):
                        all_mutations.append(f"SV_{m[1]}_{m[2]}_{m[3]}")
                    else:
                        details = '_'.join([f"{k}{v}" for k, v in m[3].items()])
                        all_mutations.append(f"SV_{m[1]}_{m[2]}_{details}")
            
            mutations_str = ','.join(all_mutations) if all_mutations else "None"
            
            # 确定完整性状态
            completeness = "complete"
            if any("deletion" in m for m in all_mutations):
                completeness = "fragmented"
            
            # 计算突变计数
            snp_count = sum(1 for m in all_mutations if m.startswith("SNP"))
            sv_count = sum(1 for m in all_mutations if m.startswith("SV"))
                
            attributes = (f"ID={te['te_id']};Name={te['family']};"
                      f"completeness={completeness};"
                      f"snp_count={snp_count};sv_count={sv_count};"
                      f"nested_in={te['nested_in'] if te['nested_in'] else 'None'}")
            
            f.write(f"chr1\tTEsimulation\ttransposable_element\t{te['start']+1}\t{te['end']+1}\t.\t"
                   f"{te['strand']}\t.\t{attributes}\n")

def write_te_sequences(sequences, output_file):
    """将TE序列写入FASTA文件"""
    with open(output_file, 'w') as f:
        for te_id, data in sequences.items():
            f.write(f">{te_id} {data['desc']}\n")
            
            # 每行写入60个碱基
            te_seq = data['seq']
            for i in range(0, len(te_seq), 60):
                f.write(te_seq[i:i+60] + '\n')

def write_fasta(sequence, output_file, chunk_size=60):
    """将序列写入FASTA格式文件"""
    with open(output_file, 'w') as f:
        f.write(">chr1 evolved genome\n")
        # 每行写入chunk_size个碱基
        for i in range(0, len(sequence), chunk_size):
            f.write(sequence[i:i+chunk_size] + '\n')

def write_evolution_log(snp_positions, sv_events, output_file):
    """将进化事件写入日志文件，用于后续分析"""
    with open(output_file, 'w') as f:
        # 写入时间戳
        f.write(f"# 进化事件日志\n")
        f.write(f"# 生成时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # 写入SNP事件
        f.write(f"## SNP事件 (共 {len(snp_positions)} 个)\n")
        f.write("Position\tOld_Base\tNew_Base\n")
        
        for pos, old_base, new_base in snp_positions:
            f.write(f"{pos}\t{old_base}\t{new_base}\n")
        
        # 写入SV事件
        f.write(f"\n## SV事件 (共 {len(sv_events)} 个)\n")
        f.write("Type\tStart\tEnd\tDetails\n")
        
        for sv_type, start, end, details in sv_events:
            # 将details转换为字符串
            if isinstance(details, dict):
                details_str = ','.join([f"{k}={v}" for k, v in details.items() if k != 'deleted_sequence'])
                if 'deleted_sequence' in details:
                    details_str += f",seq_length={len(details['deleted_sequence'])}"
            else:
                details_str = str(details)
            
            f.write(f"{sv_type}\t{start}\t{end}\t{details_str}\n")

def main():
    parser = argparse.ArgumentParser(description='模拟基因组进化过程，支持TE序列追踪和进化时间估算')
    parser.add_argument('-g', '--genome', type=str, required=True,
                        help='输入基因组FASTA文件')
    parser.add_argument('-a', '--annotations', type=str, required=True,
                        help='TE注释文件 (BED格式)')
    parser.add_argument('-s', '--te_sequences', type=str, default=None,
                        help='TE序列文件 (FASTA格式，可选)')
    parser.add_argument('-l', '--level', type=str, choices=['low', 'medium', 'high'], default='medium',
                        help='进化水平: low(SNP 1%, SV 0.5%), medium(SNP 3%, SV 1%), high(SNP 5%, SV 3%) (默认: medium)')
    parser.add_argument('-m', '--te_multiplier', type=float, default=2.0,
                        help='TE区域的变异率乘数 (默认: 2.0)')
    parser.add_argument('-o', '--output_prefix', type=str, default='evolved',
                        help='输出文件前缀 (默认: evolved)')
    parser.add_argument('-r', '--mutation_rate', type=float, default=1e-8,
                        help='基础突变率 (每碱基每代) (默认: 1e-8)')
    parser.add_argument('--seed', type=int, default=None,
                        help='随机数种子，用于结果复现 (默认: 随机)')
    parser.add_argument('--num_processes', type=int, default=None,
                        help='并行处理的进程数量 (默认: CPU核心数的75%)')
    
    args = parser.parse_args()
    
    # 参数验证
    validation_result = validate_evolution_params(args)
    validation_result.print_messages()
    validation_result.exit_if_invalid()
    
    # 设置随机种子
    if args.seed is not None:
        random.seed(args.seed)
    
    # 设置并行处理进程数
    num_processes = args.num_processes if args.num_processes else get_optimal_processes()
    print(f"使用 {num_processes} 个进程进行并行计算")
    
    # 解析输入文件
    print("解析基因组和TE注释...")
    genome_dict = parse_fasta(args.genome)
    te_annotations = parse_annotations(args.annotations)
    
    te_sequences = None
    if args.te_sequences:
        print("解析TE序列...")
        te_sequences = parse_fasta(args.te_sequences)
    
    if not genome_dict:
        sys.exit("错误: 基因组文件为空或格式不正确")
    if not te_annotations:
        sys.exit("错误: TE注释为空或格式不正确")
    
    # 获取第一条染色体序列
    chromosome_id = list(genome_dict.keys())[0]
    genome_sequence = genome_dict[chromosome_id]['seq']
    
    # 获取进化水平参数
    evolution_params = EVOLUTION_LEVELS[args.level]
    snp_rate = evolution_params['snp']
    sv_rate = evolution_params['sv']
    
    print(f"基因组长度: {len(genome_sequence):,} bp")
    print(f"TE注释数量: {len(te_annotations)}")
    print(f"进化水平: {args.level} (SNP率: {snp_rate*100}%, SV率: {sv_rate*100}%)")
    print(f"TE区域变异率乘数: {args.te_multiplier}")
    print(f"基础突变率: {args.mutation_rate}")
    
    # 提取TE区域
    te_regions = [(te['start'], te['end']) for te in te_annotations]
    
    # 引入SNP (使用并行处理)
    print("引入SNP变异...")
    modified_genome, snp_positions = introduce_snp_parallel(
        genome_sequence, 
        snp_rate, 
        te_regions, 
        args.te_multiplier, 
        args.mutation_rate,
        num_processes,
        args.seed
    )
    print(f"引入了 {len(snp_positions)} 个SNP")
    
    # 引入SV
    print("引入结构变异...")
    modified_genome, sv_events = introduce_sv(
        modified_genome, sv_rate, te_regions, args.te_multiplier, args.mutation_rate)
    print(f"引入了 {len(sv_events)} 个结构变异")
    
    # 更新TE注释
    print("更新TE注释...")
    updated_annotations = update_annotations(te_annotations, snp_positions, sv_events)
    print(f"更新后的TE注释数量: {len(updated_annotations)}")
    
    # 更新TE序列（如果提供）
    updated_sequences = None
    if te_sequences:
        print("更新TE序列...")
        updated_sequences = update_te_sequences(te_sequences, updated_annotations, modified_genome, snp_positions, sv_events)
        print(f"更新后的TE序列数量: {len(updated_sequences)}")
    
    # 创建输出目录
    os.makedirs(os.path.dirname(os.path.abspath(args.output_prefix)), exist_ok=True)
    
    # 保存结果
    genome_output = f"{args.output_prefix}_genome.fa"
    bed_output = f"{args.output_prefix}_annotations.bed"
    gff_output = f"{args.output_prefix}_annotations.gff3"
    sequences_output = f"{args.output_prefix}_te_sequences.fa"
    log_output = f"{args.output_prefix}_evolution_log.txt"
    
    print("保存结果...")
    write_fasta(modified_genome, genome_output)
    write_evolved_annotations(updated_annotations, bed_output)
    write_evolved_gff3(updated_annotations, gff_output)
    write_evolution_log(snp_positions, sv_events, log_output)
    
    if updated_sequences:
        write_te_sequences(updated_sequences, sequences_output)
    
    print(f"完成! 进化后的基因组和注释已保存")
    print(f"进化后的基因组保存至: {genome_output}")
    print(f"更新的TE注释 (BED格式) 保存至: {bed_output}")
    print(f"更新的TE注释 (GFF3格式) 保存至: {gff_output}")
    if updated_sequences:
        print(f"更新的TE序列保存至: {sequences_output}")
    print(f"进化事件日志保存至: {log_output}")

if __name__ == "__main__":
    main()

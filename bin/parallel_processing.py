#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parallel_processing.py - 为TE扩增模拟工作流提供并行处理支持

此模块提供可重用的并行处理函数，用于加速大型基因组的模拟过程。
"""

import os
import multiprocessing as mp
from functools import partial
import numpy as np
import random
from Bio import pairwise2

def get_optimal_processes(scale_factor=0.9, task_type='mixed'):
    """确定最佳进程数量，提高CPU利用率"""
    cpu_count = mp.cpu_count()
    
    # 基于任务类型智能选择进程数
    if task_type == 'io_bound':
        # I/O密集型任务使用更多进程覆盖等待时间
        processes = max(1, int(cpu_count * 1.5))
    elif task_type == 'cpu_bound':
        # 计算密集型任务轻微超订CPU
        processes = max(1, int(cpu_count * 1.1))
    else:
        # 混合型任务使用略低于核心数的进程
        processes = max(1, int(cpu_count * scale_factor))
    
    # 根据系统内存情况限制最大进程数
    try:
        import psutil
        # 如果可用内存低于2GB，则减少进程数
        if psutil.virtual_memory().available < 2 * 1024 * 1024 * 1024:
            processes = min(processes, cpu_count // 2)
    except ImportError:
        pass
    
    return processes

# TE插入相关并行处理函数
def calculate_insertion_weights_chunk(args):
    """计算一个区块的TE插入权重"""
    chunk_start, chunk_end, genome, regions, te_metadata, te_id = args
    chunk_len = chunk_end - chunk_start
    weights = [1.0] * chunk_len
    
    # 获取TE的偏好
    preferences = te_metadata.get(te_id, {})
    gc_pref = preferences.get('gc_preference', 'neutral')
    region_prefs = preferences.get('region_preference', [])
    
    # 应用区域偏好
    if regions and region_prefs:
        for region_name, (start, end) in regions.items():
            if end < chunk_start or start > chunk_end:
                continue  # 区域与当前块不重叠
                
            region_type = region_name.split('_')[0]  # 提取区域类型（不含编号）
            weight_multiplier = 5.0 if region_type in region_prefs else 0.2
            
            overlap_start = max(chunk_start, start)
            overlap_end = min(chunk_end, end + 1)
            
            if overlap_start < overlap_end:
                for i in range(overlap_start - chunk_start, overlap_end - chunk_start):
                    weights[i] *= weight_multiplier
    
    # 应用GC含量偏好
    if gc_pref != 'neutral':
        window_size = 500  # 计算GC含量的窗口大小
        step_size = window_size // 2
        
        for i in range(0, chunk_len, step_size):
            end_pos = min(i + window_size, chunk_len)
            if end_pos - i < 100:  # 忽略过小的窗口
                continue
                
            # 对于跨区块边界的窗口，尝试获取完整窗口
            window_start = chunk_start + i
            window_end = window_start + window_size
            
            if window_end <= len(genome):
                window = genome[window_start:window_end]
                gc_content = (window.count('G') + window.count('C')) / len(window)
                
                gc_multiplier = 1.0
                if gc_pref == 'high' and gc_content > 0.55:
                    gc_multiplier = 3.0
                elif gc_pref == 'high' and gc_content < 0.45:
                    gc_multiplier = 0.3
                elif gc_pref == 'low' and gc_content < 0.45:
                    gc_multiplier = 3.0
                elif gc_pref == 'low' and gc_content > 0.55:
                    gc_multiplier = 0.3
                
                for j in range(i, min(i + window_size, chunk_len)):
                    weights[j] *= gc_multiplier
    
    return (chunk_start, weights)

def calculate_insertion_weights_parallel(genome, regions, te_metadata, te_id, num_processes=None):
    """优化并行计算的权重函数"""
    if num_processes is None:
        num_processes = get_optimal_processes(task_type='cpu_bound')
    
    genome_length = len(genome)
    
    # 动态计算最佳分块大小
    chunk_size = min(500000, max(10000, genome_length // (num_processes * 4)))
    chunks = []
    
    # 将工作分配为更多小块以改善负载均衡
    for i in range(0, genome_length, chunk_size):
        end = min(i + chunk_size, genome_length)
        chunks.append((i, end, genome, regions, te_metadata, te_id))
    
    # 使用进程池的imap而非map以支持动态任务分配
    with mp.Pool(processes=num_processes) as pool:
        results = list(pool.imap_unordered(calculate_insertion_weights_chunk, chunks, chunksize=1))
    
    # 合并结果
    weights = [1.0] * genome_length
    for start, chunk_weights in results:
        for i, w in enumerate(chunk_weights):
            if start + i < genome_length:
                weights[start + i] = w
    
    return weights

# 序列变异相关并行处理函数
def introduce_snp_chunk(args):
    """在序列区块中引入SNP变异"""
    chunk, start_idx, rate, te_map, te_multiplier, base_rate, random_seed = args
    
    # 设置随机种子，确保可重现性
    if random_seed is not None:
        random.seed(random_seed + start_idx)  # 为每个区块设置不同的种子
    
    seq_list = list(chunk)
    snp_positions = []
    
    # 调整突变率，使其与基础突变率关联
    adjusted_rate = rate * (base_rate * 1e8)
    
    bases = ['A', 'T', 'C', 'G']
    for i in range(len(seq_list)):
        global_pos = start_idx + i
        
        # 确定当前位置的SNP率
        current_rate = adjusted_rate
        if te_map and global_pos in te_map:
            current_rate *= te_multiplier
        
        # 随机决定是否引入SNP
        if random.random() < current_rate:
            old_base = seq_list[i]
            if old_base in bases:
                # 确保新碱基与原碱基不同
                available_bases = [b for b in bases if b != old_base]
                new_base = random.choice(available_bases)
                seq_list[i] = new_base
                
                # 记录SNP位置
                snp_positions.append((global_pos, old_base, new_base))
    
    return (''.join(seq_list), snp_positions)

def build_te_map(te_regions):
    """构建TE区域的快速查找映射"""
    if not te_regions:
        return None
        
    te_map = {}
    for start, end in te_regions:
        for pos in range(start, end + 1):
            te_map[pos] = True
    return te_map

def introduce_snp_parallel(sequence, rate, te_regions=None, te_multiplier=1.5, base_rate=1e-8, 
                         num_processes=None, random_seed=None):
    """
    并行版本的SNP引入函数
    
    参数:
    sequence -- 基因组序列
    rate -- SNP引入率
    te_regions -- TE区域列表 [(start, end), ...]
    te_multiplier -- TE区域的变异率乘数
    base_rate -- 基础突变率
    num_processes -- 并行进程数量
    random_seed -- 随机数种子
    
    返回:
    (modified_sequence, snp_positions) -- 修改后的序列和SNP位置列表
    """
    if num_processes is None:
        num_processes = get_optimal_processes()
    
    # 构建TE区域映射以加速查找
    te_map = build_te_map(te_regions)
    
    # 将序列划分为区块
    genome_length = len(sequence)
    chunk_size = max(500000, genome_length // num_processes)
    chunks = []
    
    for i in range(0, genome_length, chunk_size):
        end = min(i + chunk_size, genome_length)
        chunk_seed = None if random_seed is None else random_seed + i
        chunks.append((sequence[i:end], i, rate, te_map, te_multiplier, base_rate, chunk_seed))
    
    # 并行处理各区块
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(introduce_snp_chunk, chunks)
    
    # 合并结果
    modified_sequence = ''
    snp_positions = []
    
    for seq, positions in results:
        modified_sequence += seq
        snp_positions.extend(positions)
    
    return modified_sequence, snp_positions

# 序列相似度计算相关并行函数
def calculate_sequence_similarity_pair(args):
    """计算一对序列的相似度"""
    te_id, inserted_seq, original_seq, info = args
    
    # 使用BioPython的pairwise2模块进行全局比对
    alignments = pairwise2.align.globalms(original_seq, inserted_seq, 2, -1, -2, -0.5, one_alignment_only=True)
    
    if alignments:
        alignment = alignments[0]
        # 计算相似度：匹配碱基数 / 比对长度
        matches = sum(a == b for a, b in zip(alignment.seqA, alignment.seqB) if a != '-' and b != '-')
        total = max(len(original_seq), len(inserted_seq))
        similarity = (matches / total) * 100
        
        # 计算身份值（完全相同的碱基比例）
        identity = sum(a == b for a, b in zip(original_seq, inserted_seq) if a != '-' and b != '-') / len(original_seq) * 100
        
        return (te_id, {
            'similarity': similarity,
            'identity': identity,
            'original_length': len(original_seq),
            'inserted_length': len(inserted_seq),
            'original_id': info['original_id'],
            'nested_in': info.get('nested_in', 'None'),
        })
    else:
        return (te_id, {
            'similarity': 0.0,
            'identity': 0.0,
            'original_length': len(original_seq),
            'inserted_length': len(inserted_seq),
            'original_id': info['original_id'],
            'nested_in': info.get('nested_in', 'None'),
        })


def calculate_sequence_similarity_batch(batch_tasks):
    """批量计算序列相似度"""
    batch_results = []

    for te_id, inserted_seq, original_seq, info in batch_tasks:
        # 简化的相似度计算
        min_len = min(len(original_seq), len(inserted_seq))

        matches = 0
        for i in range(min_len):
            if original_seq[i] == inserted_seq[i]:
                matches += 1

        similarity = (matches / min_len) * 100 if min_len > 0 else 0
        identity = (matches / len(original_seq)) * 100 if len(original_seq) > 0 else 0

        batch_results.append((te_id, {
            'similarity': similarity,
            'identity': identity,
            'original_length': len(original_seq),
            'inserted_length': len(inserted_seq),
            'original_id': info['original_id'],
            'nested_in': info.get('nested_in', 'None')
        }))

    return batch_results

# 下面是需要替换的函数
def calculate_similarity_metrics_parallel(library_seqs, inserted_seqs, num_processes=None):
    """并行计算TE序列相似度指标"""
    if num_processes is None:
        num_processes = get_optimal_processes(task_type='cpu_bound')

    # 准备计算任务
    tasks = []

    # 预处理：收集所有需要比较的序列对
    for te_id, inserted_data in inserted_seqs.items():
        inserted_seq = inserted_data['seq']
        desc = inserted_data['desc']

        # 从描述中提取信息
        info = {'original_id': te_id.rsplit('_', 1)[0], 'nested_in': 'None'}
        for item in desc.split():
            if '=' in item:
                key, value = item.split('=', 1)
                info[key] = value

        # 获取原始TE序列
        original_id = info['original_id']
        if original_id in library_seqs:
            original_seq = library_seqs[original_id]['seq']

            # 如果是反向互补插入，需要调整原始序列
            if info.get('strand') == '-':
                original_seq = reverse_complement(original_seq)

            tasks.append((te_id, inserted_seq, original_seq, info))

    # 批处理任务
    total_tasks = len(tasks)
    batch_size = max(1, min(100, total_tasks // (num_processes * 2)))

    # 对任务进行分组
    batched_tasks = []
    for i in range(0, len(tasks), batch_size):
        batched_tasks.append(tasks[i:i+batch_size])

    # 使用进程池并行处理批次任务
    print(f"TE序列相似度计算：将{total_tasks}个任务分为{len(batched_tasks)}批进行处理")

    # 使用map而非imap_unordered来避免某些版本的Python中的序列化问题
    with mp.Pool(processes=num_processes) as pool:
        all_batch_results = pool.map(calculate_sequence_similarity_batch, batched_tasks)

        # 收集结果
        metrics = {}
        for batch_results in all_batch_results:
            for te_id, data in batch_results:
                metrics[te_id] = data

    return metrics

def reverse_complement(sequence):
    """返回DNA序列的反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

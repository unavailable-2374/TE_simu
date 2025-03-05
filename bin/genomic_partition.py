#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
genomic_partition.py - 修复版基因组分区并行处理模块

此脚本修复了多进程序列化问题和其他逻辑错误。
"""

import os
import sys
import argparse
import subprocess
import tempfile
import random
import multiprocessing as mp
from functools import partial
import shutil
import time

def split_genome(genome_file, num_partitions):
    """将基因组分割为多个区段"""
    # 读取基因组
    with open(genome_file, 'r') as f:
        header = f.readline().strip()
        sequence = ''
        for line in f:
            sequence += line.strip()
    
    # 计算每个分区的大小
    genome_length = len(sequence)
    partition_size = genome_length // num_partitions
    
    partitions = []
    for i in range(num_partitions):
        start = i * partition_size
        end = (i + 1) * partition_size if i < num_partitions - 1 else genome_length
        partition_seq = sequence[start:end]
        
        # 创建临时文件存储分区
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.fa', mode='w')
        temp_file.write(f"{header}_part_{i+1}\n")
        # 分块写入序列，每行60个碱基
        for j in range(0, len(partition_seq), 60):
            temp_file.write(f"{partition_seq[j:j+60]}\n")
        temp_file.close()
        
        partitions.append({
            'file': temp_file.name,
            'start': start,
            'end': end,
            'length': end - start,
            'id': i
        })
    
    return partitions

# 定义为模块级函数而非局部函数，以解决pickle问题
def process_partition(args):
    """处理单个基因组分区"""
    partition, partition_seed, te_library, percentage, use_weights, config = args
    partition_file = partition['file']
    partition_id = partition['id']
    
    # 计算比例TE百分比
    part_percentage = percentage
    
    # 创建临时输出目录
    temp_dir = tempfile.mkdtemp()
    part_output = os.path.join(temp_dir, f"part_{partition_id}")
    
    # 构建命令
    cmd_parts = [
        "python", "optimized_te_insertion.py",
        "-g", partition_file,
        "-t", te_library,
        "-p", str(part_percentage),
        "-o", part_output,
        "--num_processes", "1"  # 每个分区用1个进程，避免嵌套并行
    ]
    
    if use_weights:
        cmd_parts.append("--use_weights")
        if config:
            cmd_parts.extend(["-c", config])
    
    if partition_seed is not None:
        cmd_parts.extend(["--seed", str(partition_seed)])
    
    # 执行命令
    try:
        subprocess.run(cmd_parts, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # 收集结果
        genome_output = f"{part_output}_genome.fa"
        annotations_output = f"{part_output}_te_annotations.bed"
        
        # 调整注释坐标
        if os.path.exists(annotations_output):
            adjusted_annotations = []
            with open(annotations_output, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        adjusted_annotations.append(line)
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) >= 4:
                        # 调整起始和结束位置
                        fields[2] = str(int(fields[2]) + partition['start'])
                        fields[3] = str(int(fields[3]) + partition['start'])
                        adjusted_annotations.append('\t'.join(fields))
            
            # 保存调整后的注释
            with open(annotations_output, 'w') as f:
                for line in adjusted_annotations:
                    f.write(f"{line}\n")
        
        return {
            'genome': genome_output,
            'annotations': annotations_output if os.path.exists(annotations_output) else None,
            'partition': partition
        }
    except Exception as e:
        print(f"处理分区 {partition_id} 出错: {e}")
        return None

def parallelize_te_insertion(genome_file, te_library, percentage, partitions, output_prefix, 
                          use_weights=False, config=None, seed=None, num_processes=8):
    """并行在基因组分区上插入TE"""
    # 为每个分区分配对应的TE百分比
    total_length = sum(p['length'] for p in partitions)
    # 添加总长度信息到每个分区
    for p in partitions:
        p['total_length'] = total_length
    
    # 处理每个分区
    print(f"并行处理 {len(partitions)} 个基因组分区...")
    
    # 为每个分区设置不同的种子
    partition_seeds = None
    if seed is not None:
        random.seed(seed)
        partition_seeds = [random.randint(1, 100000) for _ in range(len(partitions))]
    
    # 准备参数
    args = []
    for i, partition in enumerate(partitions):
        part_seed = partition_seeds[i] if partition_seeds else None
        args.append((partition, part_seed, te_library, percentage, use_weights, config))
    
    # 并行处理
    pool = mp.Pool(processes=min(len(partitions), mp.cpu_count()))
    results = pool.map(process_partition, args)
    pool.close()
    pool.join()
    
    # 过滤出非None的结果
    valid_results = [r for r in results if r is not None]
    
    if not valid_results:
        raise Exception("所有分区处理都失败了")
    
    # 合并结果
    merged_genome_file = f"{output_prefix}_genome.fa"
    merged_annotations_file = f"{output_prefix}_te_annotations.bed"
    
    # 提取并合并序列
    with open(merged_genome_file, 'w') as f_out:
        f_out.write(">chr1 merged genome with inserted TEs\n")
        
        for result in sorted(valid_results, key=lambda x: x['partition']['start']):
            if os.path.exists(result['genome']):
                # 跳过头部
                with open(result['genome'], 'r') as f_in:
                    next(f_in)  # 跳过FASTA头
                    for line in f_in:
                        f_out.write(line)
    
    # 合并注释
    with open(merged_annotations_file, 'w') as f_out:
        f_out.write("# Merged TE annotations\n")
        f_out.write("# TE_ID\tFamily\tStart\tEnd\tStrand\tLength\tNested_In\tMutations\n")
        
        # 调整TE ID以避免冲突
        te_counter = {}
        combined_annotations = []
        
        for i, result in enumerate(sorted(valid_results, key=lambda x: x['partition']['start'])):
            if result['annotations'] and os.path.exists(result['annotations']):
                with open(result['annotations'], 'r') as f_in:
                    # 跳过头部
                    for line in f_in:
                        if not line.startswith('#'):
                            fields = line.strip().split('\t')
                            if len(fields) >= 2:
                                # 调整TE ID
                                te_id_parts = fields[0].rsplit('_', 1)
                                original_id = te_id_parts[0]
                                if original_id not in te_counter:
                                    te_counter[original_id] = 0
                                
                                te_counter[original_id] += 1
                                fields[0] = f"{original_id}_{te_counter[original_id]}"
                                
                                # 添加到组合列表
                                combined_annotations.append(fields)
        
        # 写入组合的注释
        for fields in combined_annotations:
            f_out.write('\t'.join(fields) + '\n')
    
    # 清理临时文件
    for partition in partitions:
        if os.path.exists(partition['file']):
            os.remove(partition['file'])
    
    for result in valid_results:
        if os.path.exists(os.path.dirname(result['genome'])):
            shutil.rmtree(os.path.dirname(result['genome']), ignore_errors=True)
    
    return merged_genome_file, merged_annotations_file

def parallelize_evolution(input_genome, annotations, evolution_level, mutation_rate, te_multiplier, 
                        output_prefix, num_processes=8, seed=None):
    """使用分区并行模拟进化"""
    # 使用现有的simulate_evolution.py脚本
    cmd_parts = [
        "python", "simulate_evolution.py",
        "-g", input_genome,
        "-a", annotations,
        "-l", evolution_level,
        "-m", str(te_multiplier),
        "-r", str(mutation_rate),
        "-o", output_prefix,
        "--num_processes", str(num_processes)
    ]
    
    if seed is not None:
        cmd_parts.extend(["--seed", str(seed)])
    
    # 执行命令
    try:
        print(f"执行命令: {' '.join(cmd_parts)}")
        subprocess.run(cmd_parts, check=True)
    except subprocess.CalledProcessError as e:
        print(f"进化模拟失败: {e}")
        return None, None, None
    
    # 返回输出文件路径
    evolved_genome = f"{output_prefix}_genome.fa"
    evolved_anno = f"{output_prefix}_annotations.bed"
    evolved_seq = f"{output_prefix}_te_sequences.fa"
    
    # 验证文件是否存在
    if not os.path.exists(evolved_genome):
        print(f"警告: 进化后的基因组文件 {evolved_genome} 不存在")
    if not os.path.exists(evolved_anno):
        print(f"警告: 进化后的注释文件 {evolved_anno} 不存在")
    
    return evolved_genome, evolved_anno, evolved_seq

def simulate_stage(stage, input_genome, te_library, annotations, te_sequences, 
                 output_dir, te_percentage, evolution_level, config_file, 
                 te_multiplier, mutation_rate, seed, num_partitions, use_weights=False):
    """模拟一个进化阶段，包括TE插入和变异引入"""
    stage_dir = os.path.join(output_dir, f"stage_{stage}")
    os.makedirs(stage_dir, exist_ok=True)
    
    stage_prefix = os.path.join(stage_dir, f"stage_{stage}")
    
    print(f"\n=== 阶段 {stage} (分区并行) ===")
    
    # 设置随机种子
    stage_seed = seed + stage if seed is not None else None
    
    # 1. 创建分区并插入TE
    if te_percentage > 0:
        print(f"将基因组分为 {num_partitions} 个分区进行TE插入...")
        partitions = split_genome(input_genome, num_partitions)
        
        # 2. TE插入
        print(f"并行插入TE (目标占比: {te_percentage}%)...")
        te_output_prefix = f"{stage_prefix}_te_inserted"
        genome_with_te, annotations_with_te = parallelize_te_insertion(
            input_genome, te_library, te_percentage, partitions, 
            te_output_prefix, use_weights, config_file, stage_seed, num_partitions
        )
        
        # 验证TE插入结果
        if not os.path.exists(genome_with_te):
            print(f"错误: TE插入后的基因组文件 {genome_with_te} 不存在")
            return None, None, None
        
        if not os.path.exists(annotations_with_te):
            print(f"错误: TE插入后的注释文件 {annotations_with_te} 不存在")
            return None, None, None
        
        print(f"TE插入完成，基因组文件: {genome_with_te}")
        print(f"TE注释文件: {annotations_with_te}")
    else:
        # 如果不需要插入TE，直接使用输入文件
        print("跳过TE插入阶段 (TE百分比为0)")
        genome_with_te = input_genome
        annotations_with_te = annotations
    
    # 如果没有注释文件，无法继续进化模拟
    if annotations_with_te is None:
        print("警告: 没有找到有效的TE注释文件，无法进行进化模拟")
        return genome_with_te, None, None
    
    # 3. 进化模拟
    print(f"引入变异 (进化水平: {evolution_level})...")
    evolved_prefix = f"{stage_prefix}_evolved"
    evolved_genome, evolved_anno, evolved_seq = parallelize_evolution(
        genome_with_te, annotations_with_te, evolution_level, mutation_rate, 
        te_multiplier, evolved_prefix, num_processes=num_partitions, seed=stage_seed
    )
    
    if evolved_genome is None:
        print("进化模拟失败，返回TE插入后的文件")
        return genome_with_te, annotations_with_te, None
    
    print(f"进化模拟完成，基因组文件: {evolved_genome}")
    print(f"注释文件: {evolved_anno}")
    
    return evolved_genome, evolved_anno, evolved_seq

def main():
    parser = argparse.ArgumentParser(description='分区并行TE扩增模拟 (修复版)')
    parser.add_argument('--genome', type=str, required=True, help='基因组FASTA文件')
    parser.add_argument('--te_library', type=str, required=True, help='TE文库FASTA文件')
    parser.add_argument('--stages', type=int, default=3, help='进化阶段数')
    parser.add_argument('--te_percentages', type=str, required=True, help='每个阶段的TE百分比(逗号分隔)')
    parser.add_argument('--evolution_levels', type=str, required=True, help='每个阶段的进化水平(逗号分隔)')
    parser.add_argument('--output_dir', type=str, required=True, help='输出目录')
    parser.add_argument('--config', type=str, help='基因组配置文件')
    parser.add_argument('--seed', type=int, help='随机数种子')
    parser.add_argument('--num_partitions', type=int, default=4, help='基因组分区数量')
    parser.add_argument('--mutation_rate', type=float, default=1e-8, help='突变率')
    parser.add_argument('--te_multiplier', type=float, default=2.0, help='TE区域突变率乘数')
    parser.add_argument('--fast_mode', type=str, default="False", help='使用快速模式')
    
    args = parser.parse_args()
    
    # 解析参数
    te_percentages = [float(p) for p in args.te_percentages.split(',')]
    evolution_levels = args.evolution_levels.split(',')
    
    # 转换fast_mode
    use_weights = True
    if args.fast_mode.lower() in ('true', 't', 'yes', 'y', '1'):
        use_weights = False
    
    # 确保参数数量一致
    if len(te_percentages) < args.stages:
        te_percentages.extend([te_percentages[-1]] * (args.stages - len(te_percentages)))
    if len(evolution_levels) < args.stages:
        evolution_levels.extend([evolution_levels[-1]] * (args.stages - len(evolution_levels)))
    
    # 截断多余的参数
    te_percentages = te_percentages[:args.stages]
    evolution_levels = evolution_levels[:args.stages]
    
    # 记录每个阶段的信息
    stages_info = []
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 初始文件
    current_genome = args.genome
    current_annotations = None
    current_te_sequences = None
    
    # 处理每个阶段
    for i in range(args.stages):
        stage_info = {
            'stage': i+1,
            'te_percentage': te_percentages[i],
            'evolution_level': evolution_levels[i]
        }
        
        # 模拟当前阶段
        start_time = time.time()
        
        # 打印当前状态
        print(f"开始阶段 {i+1}:")
        print(f"  使用基因组文件: {current_genome}")
        print(f"  使用注释文件: {current_annotations}")
        print(f"  使用TE序列文件: {current_te_sequences}")
        print(f"  TE百分比: {te_percentages[i]}%")
        print(f"  进化水平: {evolution_levels[i]}")
        
        evolved_genome, evolved_annotations, evolved_sequences = simulate_stage(
            i+1, current_genome, args.te_library, current_annotations, current_te_sequences,
            args.output_dir, te_percentages[i], evolution_levels[i], args.config,
            args.te_multiplier, args.mutation_rate, args.seed, args.num_partitions, use_weights
        )
        
        # 检查模拟是否成功
        if evolved_genome is None:
            print(f"阶段 {i+1} 模拟失败，终止模拟")
            break
        
        end_time = time.time()
        print(f"阶段 {i+1} 完成，用时: {end_time - start_time:.2f} 秒")
        
        # 更新文件路径用于下一阶段
        current_genome = evolved_genome
        current_annotations = evolved_annotations
        current_te_sequences = evolved_sequences
        
        # 记录该阶段的输出
        stage_info['genome'] = evolved_genome
        stage_info['annotations'] = evolved_annotations
        stage_info['te_sequences'] = evolved_sequences
        stages_info.append(stage_info)
    
    # 检查是否完成了所有阶段
    if len(stages_info) < args.stages:
        print(f"警告: 只完成了 {len(stages_info)}/{args.stages} 个阶段")
    
    if not stages_info:
        print("错误: 没有完成任何阶段，无法生成最终结果")
        return
    
    # 为最终结果创建符号链接
    final_genome = os.path.join(args.output_dir, "final_genome.fa")
    final_annotations = os.path.join(args.output_dir, "final_annotations.bed")
    final_gff = os.path.join(args.output_dir, "final_annotations.gff3")
    final_te_sequences = os.path.join(args.output_dir, "final_te_sequences.fa")
    
    # 删除已存在的文件
    for file_path in [final_genome, final_annotations, final_gff, final_te_sequences]:
        if os.path.exists(file_path):
            os.remove(file_path)
    
    # 复制最终文件
    shutil.copy2(current_genome, final_genome)
    
    if current_annotations and os.path.exists(current_annotations):
        shutil.copy2(current_annotations, final_annotations)
    
    if current_te_sequences and os.path.exists(current_te_sequences):
        shutil.copy2(current_te_sequences, final_te_sequences)
    
    # 复制GFF3文件
    if current_annotations:
        gff_path = current_annotations.replace(".bed", ".gff3")
        if os.path.exists(gff_path):
            shutil.copy2(gff_path, final_gff)
    
    # 生成摘要报告
    summary_file = os.path.join(args.output_dir, "evolution_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("# 模拟进化过程摘要\n")
        f.write(f"生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"基础突变率: {args.mutation_rate}\n\n")
        
        for i, stage in enumerate(stages_info):
            f.write(f"## 阶段 {i+1}\n")
            f.write(f"TE插入比例: {stage['te_percentage']}%\n")
            f.write(f"进化水平: {stage['evolution_level']}\n")
            f.write(f"基因组文件: {stage['genome']}\n")
            f.write(f"注释文件: {stage['annotations'] if stage['annotations'] else 'None'}\n")
            if 'te_sequences' in stage and stage['te_sequences']:
                f.write(f"TE序列文件: {stage['te_sequences']}\n")
            f.write("\n")
        
        f.write("## 最终结果\n")
        final_stage = stages_info[-1]
        f.write(f"最终基因组: {final_genome}\n")
        f.write(f"最终注释: {final_annotations}\n")
        if 'te_sequences' in final_stage and final_stage['te_sequences']:
            f.write(f"最终TE序列: {final_te_sequences}\n")
    
    # 检查最终文件大小差异
    try:
        base_size = os.path.getsize(args.genome)
        final_size = os.path.getsize(final_genome)
        size_diff = final_size - base_size
        
        print(f"\n基础基因组大小: {base_size:,} 字节")
        print(f"最终基因组大小: {final_size:,} 字节")
        print(f"大小差异: {size_diff:,} 字节 ({(size_diff/base_size)*100:.2f}%)")
        
        if abs(size_diff) < 0.01 * base_size:
            print("警告: 最终基因组与初始基因组大小几乎相同，可能没有正确插入TE")
    except Exception as e:
        print(f"无法比较文件大小: {e}")
    
    print("\n=== 模拟完成 ===")
    print(f"总共完成 {len(stages_info)} 个进化阶段")
    print(f"最终基因组: {final_genome}")
    print(f"最终注释文件: {final_annotations}")
    print(f"最终TE序列文件: {final_te_sequences}")
    print(f"摘要报告: {summary_file}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
time_evolution.py - 模拟时间线上的TE插入和变异过程

此脚本模拟多个进化阶段，每个阶段包括TE插入和变异引入，以更真实地模拟基因组进化过程。
增强版本考虑了染色体结构和TE插入偏好，并记录详细的序列变化以支持进化时间分析。
"""

import os
import sys
import argparse
import random
import subprocess
import datetime
import json
import shutil
# 导入参数验证模块
from parameter_validation import validate_timeline_params
from parallel_processing import get_optimal_processes

def run_command(command):
    """执行外部命令并返回输出"""
    try:
        result = subprocess.run(
            command, shell=True, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"命令执行失败: {e}")
        print(f"错误信息: {e.stderr}")
        sys.exit(1)

def simulate_stage(stage, input_genome, te_library, annotations, te_sequences, 
                 output_dir, te_percentage, evolution_level, config_file, 
                 te_multiplier, mutation_rate, seed, num_processes, use_weights=False):
    """
    模拟一个进化阶段，包括TE插入和变异引入
    
    参数:
    stage -- 当前阶段编号
    input_genome -- 输入基因组文件
    te_library -- TE文库文件
    annotations -- 当前TE注释文件 (如果有)
    te_sequences -- 当前TE序列文件 (如果有)
    output_dir -- 输出目录
    te_percentage -- 该阶段要插入的TE百分比
    evolution_level -- 进化水平 (low/medium/high)
    config_file -- 基因组配置文件
    te_multiplier -- TE区域的变异率乘数
    mutation_rate -- 突变率参数
    seed -- 随机数种子
    num_processes -- 并行处理的进程数量
    use_weights -- 是否使用权重计算 (默认: False)
    
    返回:
    (evolved_genome, evolved_annotations, evolved_sequences) -- 该阶段结束后的文件路径
    """
    stage_dir = os.path.join(output_dir, f"stage_{stage}")
    os.makedirs(stage_dir, exist_ok=True)
    
    stage_prefix = os.path.join(stage_dir, f"stage_{stage}")
    
    print(f"\n=== 阶段 {stage} ===")
    
    # 准备并行处理参数
    processes_param = f"--num_processes {num_processes}" if num_processes else ""
    
    # 设置权重参数
    fast_param = "--use_weights" if use_weights else ""
    
    # 1. TE插入
    if te_percentage > 0:
        print(f"插入TE (目标占比: {te_percentage}%)...")
        
        te_output_prefix = f"{stage_prefix}_te_inserted"
        te_output_genome = f"{te_output_prefix}_genome.fa"
        te_output_anno = f"{te_output_prefix}_te_annotations.bed"
        te_output_seq = f"{te_output_prefix}_te_sequences.fa"
        
        seed_param = f"--seed {seed}" if seed is not None else ""
        config_param = f"-c {config_file}" if config_file and os.path.exists(config_file) else ""
        
        fast_param = "--use_weights" if use_weights else ""
        cmd = (f"python optimized_te_insertion.py -g {input_genome} -t {te_library} "
            f"-p {te_percentage} -o {te_output_prefix} {config_param} {seed_param} {processes_param} {fast_param}") 
        run_command(cmd)
    else:
        # 如果不需要插入TE，直接使用输入文件
        te_output_genome = input_genome
        te_output_anno = annotations
        te_output_seq = te_sequences if te_sequences else None
    
    # 2. 引入变异
    print(f"引入变异 (进化水平: {evolution_level}, 突变率基数: {mutation_rate})...")
    
    evolved_prefix = f"{stage_prefix}_evolved"
    evolved_genome = f"{evolved_prefix}_genome.fa"
    evolved_anno = f"{evolved_prefix}_annotations.bed"
    evolved_seq = f"{evolved_prefix}_te_sequences.fa"
    
    cmd = (f"python simulate_evolution.py -g {te_output_genome} -a {te_output_anno} "
          f"-l {evolution_level} -m {te_multiplier} -r {mutation_rate} -o {evolved_prefix} {processes_param}")
    
    if te_output_seq and os.path.exists(te_output_seq):
        cmd += f" -s {te_output_seq}"
    
    if seed is not None:
        cmd += f" --seed {seed + stage}"  # 每个阶段使用不同的种子
    
    run_command(cmd)
    
    return evolved_genome, evolved_anno, evolved_seq

def generate_summary(stages_info, output_file, mutation_rate):
    """生成进化过程摘要"""
    with open(output_file, 'w') as f:
        f.write("# 模拟进化过程摘要\n")
        f.write(f"生成时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"基础突变率: {mutation_rate}\n\n")
        
        for i, stage in enumerate(stages_info):
            f.write(f"## 阶段 {i+1}\n")
            f.write(f"TE插入比例: {stage['te_percentage']}%\n")
            f.write(f"进化水平: {stage['evolution_level']}\n")
            f.write(f"基因组文件: {stage['genome']}\n")
            f.write(f"注释文件: {stage['annotations']}\n")
            if 'te_sequences' in stage and stage['te_sequences']:
                f.write(f"TE序列文件: {stage['te_sequences']}\n")
            f.write("\n")
        
        f.write("## 最终结果\n")
        final_stage = stages_info[-1]
        f.write(f"最终基因组: {final_stage['genome']}\n")
        f.write(f"最终注释: {final_stage['annotations']}\n")
        if 'te_sequences' in final_stage and final_stage['te_sequences']:
            f.write(f"最终TE序列: {final_stage['te_sequences']}\n")

def main():
    parser = argparse.ArgumentParser(description='模拟时间线上的TE插入和变异过程，支持染色体结构和进化时间分析')
    parser.add_argument('-g', '--base_genome', type=str, required=True,
                        help='初始基因组FASTA文件')
    parser.add_argument('-t', '--te_library', type=str, required=True,
                        help='TE文库FASTA文件')
    parser.add_argument('-s', '--stages', type=int, default=3,
                        help='进化阶段数量 (默认: 3)')
    parser.add_argument('-p', '--te_percentages', type=str, default="10,10,10",
                        help='每个阶段的TE插入比例，逗号分隔 (默认: 10,10,10)')
    parser.add_argument('-l', '--evolution_levels', type=str, default="low,medium,high",
                        help='每个阶段的进化水平，逗号分隔 (默认: low,medium,high)')
    parser.add_argument('-m', '--te_multiplier', type=float, default=2.0,
                        help='TE区域的变异率乘数 (默认: 2.0)')
    parser.add_argument('-o', '--output_dir', type=str, default='evolution_output',
                        help='输出目录 (默认: evolution_output)')
    parser.add_argument('-c', '--config', type=str, default=None,
                        help='基因组配置文件 (包含染色体结构信息)')
    parser.add_argument('-r', '--mutation_rate', type=float, default=1e-8,
                        help='基础突变率 (每碱基每代) (默认: 1e-8)')
    parser.add_argument('--seed', type=int, default=None,
                        help='随机数种子，用于结果复现 (默认: 随机)')
    parser.add_argument('--num_processes', type=int, default=None,
                        help='并行处理的进程数量 (默认: CPU核心数的75%)')
    parser.add_argument('--fast_mode', action='store_true',
                        help='使用快速TE插入模式 (不考虑区域偏好和GC含量)')
    
    args = parser.parse_args()
    
    # 参数验证
    validation_result = validate_timeline_params(args)
    validation_result.print_messages()
    validation_result.exit_if_invalid()
    
    # 设置随机种子
    if args.seed is not None:
        random.seed(args.seed)
    
    # 设置并行处理进程数
    num_processes = args.num_processes if args.num_processes else get_optimal_processes()
    print(f"使用 {num_processes} 个进程进行并行计算")
    
    # 解析每个阶段的参数
    te_percentages = [float(p) for p in args.te_percentages.split(',')]
    evolution_levels = args.evolution_levels.split(',')
    
    # 确保参数数量一致
    if len(te_percentages) < args.stages:
        te_percentages.extend([te_percentages[-1]] * (args.stages - len(te_percentages)))
    if len(evolution_levels) < args.stages:
        evolution_levels.extend([evolution_levels[-1]] * (args.stages - len(evolution_levels)))
    
    # 截断多余的参数
    te_percentages = te_percentages[:args.stages]
    evolution_levels = evolution_levels[:args.stages]
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f"模拟 {args.stages} 个进化阶段...")
    print(f"TE插入比例: {', '.join(map(str, te_percentages))}%")
    print(f"进化水平: {', '.join(evolution_levels)}")
    print(f"基础突变率: {args.mutation_rate}")
    print(f"快速模式: {'启用' if args.fast_mode else '禁用'}")
    
    # 初始文件
    current_genome = args.base_genome
    current_annotations = None  # 初始没有注释
    current_te_sequences = None  # 初始没有TE序列
    
    # 记录每个阶段的信息
    stages_info = []
    
    # 模拟每个阶段
    for i in range(args.stages):
        stage_info = {
            'stage': i+1,
            'te_percentage': te_percentages[i],
            'evolution_level': evolution_levels[i]
        }
        
        # 运行当前阶段
        evolved_genome, evolved_annotations, evolved_sequences = simulate_stage(
            i+1, current_genome, args.te_library, current_annotations, current_te_sequences,
            args.output_dir, te_percentages[i], evolution_levels[i], args.config,
            args.te_multiplier, args.mutation_rate, args.seed, num_processes, not args.fast_mode
        )
        
        # 更新文件路径用于下一阶段
        current_genome = evolved_genome
        current_annotations = evolved_annotations
        current_te_sequences = evolved_sequences
        
        # 记录该阶段的输出
        stage_info['genome'] = evolved_genome
        stage_info['annotations'] = evolved_annotations
        stage_info['te_sequences'] = evolved_sequences
        stages_info.append(stage_info)
    
    # 为最终结果创建符号链接
    final_genome = os.path.join(args.output_dir, "final_genome.fa")
    final_annotations = os.path.join(args.output_dir, "final_annotations.bed")
    final_gff = os.path.join(args.output_dir, "final_annotations.gff3")
    final_te_sequences = os.path.join(args.output_dir, "final_te_sequences.fa")
    
    # 删除已存在的文件（防止覆盖错误）
    for file_path in [final_genome, final_annotations, final_gff, final_te_sequences]:
        if os.path.exists(file_path):
            os.remove(file_path)
    
    # 复制最终文件（而不是符号链接，以防文件被移动）
    shutil.copy2(current_genome, final_genome)
    shutil.copy2(current_annotations, final_annotations)
    
    if current_te_sequences and os.path.exists(current_te_sequences):
        shutil.copy2(current_te_sequences, final_te_sequences)
    
    # 复制GFF3文件
    gff_path = current_annotations.replace(".bed", ".gff3")
    if os.path.exists(gff_path):
        shutil.copy2(gff_path, final_gff)
    
    # 生成摘要报告
    summary_file = os.path.join(args.output_dir, "evolution_summary.txt")
    generate_summary(stages_info, summary_file, args.mutation_rate)
    
    print("\n=== 模拟完成 ===")
    print(f"总共完成 {args.stages} 个进化阶段")
    print(f"最终基因组: {final_genome}")
    print(f"最终注释文件: {final_annotations}")
    print(f"最终TE序列文件: {final_te_sequences}")
    print(f"摘要报告: {summary_file}")

if __name__ == "__main__":
    main()

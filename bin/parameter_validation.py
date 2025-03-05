#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parameter_validation.py - 为TE扩增模拟工作流提供参数验证功能

此模块提供参数验证函数，确保输入参数的合理性，提高工作流的健壮性。
"""

import os
import sys
import re
import math
from collections import defaultdict

class ValidationResult:
    """参数验证结果类"""
    def __init__(self):
        self.errors = []
        self.warnings = []
    
    def add_error(self, message):
        """添加错误消息"""
        self.errors.append(message)
    
    def add_warning(self, message):
        """添加警告消息"""
        self.warnings.append(message)
    
    def is_valid(self):
        """检查是否通过验证（无错误）"""
        return len(self.errors) == 0
    
    def print_messages(self):
        """打印所有错误和警告消息"""
        if self.errors:
            print("\n[错误] 参数验证失败:")
            for error in self.errors:
                print(f"  - {error}")
        
        if self.warnings:
            print("\n[警告] 参数验证产生警告:")
            for warning in self.warnings:
                print(f"  - {warning}")
        
        if self.is_valid() and not self.warnings:
            print("\n[信息] 所有参数验证通过")
        elif self.is_valid():
            print("\n[信息] 参数验证通过，但请注意上述警告")
    
    def exit_if_invalid(self):
        """如果验证失败则退出程序"""
        if not self.is_valid():
            self.print_messages()
            sys.exit(1)

# 基因组生成参数验证
def validate_genome_generation_params(args):
    """
    验证基因组生成参数
    
    参数:
    args -- 参数命名空间，必须包含:
            length, gc, gc_gradient, min_gc, max_gc, 
            telomere_size, centromere_size
    
    返回:
    result -- ValidationResult对象
    """
    result = ValidationResult()
    
    # 验证基因组长度
    if args.length <= 0:
        result.add_error(f"基因组长度必须为正数，当前值: {args.length}")
    elif args.length < 10000:
        result.add_warning(f"基因组长度过小 ({args.length} bp)，可能不足以模拟真实染色体结构")
    elif args.length > 1e9:
        result.add_warning(f"基因组长度过大 ({args.length:,} bp)，可能导致内存不足或运行时间过长")
    
    # 验证GC含量
    if args.gc < 0 or args.gc > 1:
        result.add_error(f"GC含量必须在0-1范围内，当前值: {args.gc}")
    
    # 验证GC梯度参数
    if args.gc_gradient:
        if args.min_gc < 0 or args.min_gc > 1:
            result.add_error(f"最小GC含量必须在0-1范围内，当前值: {args.min_gc}")
        if args.max_gc < 0 or args.max_gc > 1:
            result.add_error(f"最大GC含量必须在0-1范围内，当前值: {args.max_gc}")
        if args.min_gc >= args.max_gc:
            result.add_error(f"最小GC含量 ({args.min_gc}) 必须小于最大GC含量 ({args.max_gc})")
    
    # 验证染色体区域大小
    if args.telomere_size < 0:
        result.add_error(f"端粒区域大小不能为负数，当前值: {args.telomere_size}")
    if args.centromere_size < 0:
        result.add_error(f"着丝粒区域大小不能为负数，当前值: {args.centromere_size}")
    
    # 验证区域大小与基因组长度的关系
    total_special_regions = 2 * args.telomere_size + args.centromere_size
    if total_special_regions >= args.length:
        result.add_error(f"特殊区域总大小 ({total_special_regions:,} bp) 不能大于等于基因组长度 ({args.length:,} bp)")
    elif total_special_regions > args.length * 0.5:
        result.add_warning(f"特殊区域占比过高 ({total_special_regions/args.length:.2%})，可能导致普通染色体臂区域过小")
    
    return result

# TE插入参数验证
def validate_te_insertion_params(args):
    """
    验证TE插入参数
    
    参数:
    args -- 参数命名空间，必须包含:
            genome, te_library, percentage, output_prefix, config(可选)
    
    返回:
    result -- ValidationResult对象
    """
    result = ValidationResult()
    
    # 验证输入文件存在
    if not os.path.exists(args.genome):
        result.add_error(f"基因组文件不存在: {args.genome}")
    if not os.path.exists(args.te_library):
        result.add_error(f"TE文库文件不存在: {args.te_library}")
    if hasattr(args, 'config') and args.config and not os.path.exists(args.config):
        result.add_warning(f"配置文件不存在: {args.config}，将使用默认参数")
    
    # 验证TE占比
    if args.percentage <= 0:
        result.add_error(f"TE占比必须为正数，当前值: {args.percentage}%")
    elif args.percentage > 90:
        result.add_warning(f"TE占比过高 ({args.percentage}%)，自然基因组中TE占比通常不超过70%")
    
    # 验证输出目录
    output_dir = os.path.dirname(os.path.abspath(args.output_prefix))
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except PermissionError:
            result.add_error(f"无权限创建输出目录: {output_dir}")
    
    # 检查输入文件格式
    if os.path.exists(args.genome):
        try:
            with open(args.genome, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    result.add_error(f"基因组文件不是有效的FASTA格式: {args.genome}")
        except UnicodeDecodeError:
            result.add_error(f"基因组文件不是有效的文本文件: {args.genome}")
        except Exception as e:
            result.add_error(f"读取基因组文件时出错: {e}")
    
    if os.path.exists(args.te_library):
        try:
            with open(args.te_library, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    result.add_error(f"TE文库文件不是有效的FASTA格式: {args.te_library}")
        except UnicodeDecodeError:
            result.add_error(f"TE文库文件不是有效的文本文件: {args.te_library}")
        except Exception as e:
            result.add_error(f"读取TE文库文件时出错: {e}")
    
    return result

# 演化模拟参数验证
def validate_evolution_params(args):
    """
    验证演化模拟参数
    
    参数:
    args -- 参数命名空间，必须包含:
            genome, annotations, level, te_multiplier, 
            mutation_rate, output_prefix
    
    返回:
    result -- ValidationResult对象
    """
    result = ValidationResult()
    
    # 验证输入文件存在
    if not os.path.exists(args.genome):
        result.add_error(f"基因组文件不存在: {args.genome}")
    if not os.path.exists(args.annotations):
        result.add_error(f"注释文件不存在: {args.annotations}")
    if hasattr(args, 'te_sequences') and args.te_sequences and not os.path.exists(args.te_sequences):
        result.add_warning(f"TE序列文件不存在: {args.te_sequences}，将跳过TE序列追踪")
    
    # 验证进化水平
    valid_levels = ['low', 'medium', 'high']
    if args.level not in valid_levels:
        result.add_error(f"无效的进化水平: {args.level}，有效选项为 {', '.join(valid_levels)}")
    
    # 验证TE变异率乘数
    if args.te_multiplier <= 0:
        result.add_error(f"TE区域变异率乘数必须为正数，当前值: {args.te_multiplier}")
    elif args.te_multiplier > 10:
        result.add_warning(f"TE区域变异率乘数过高 ({args.te_multiplier})，自然情况下通常为1.5-3")
    
    # 验证突变率
    if args.mutation_rate <= 0:
        result.add_error(f"突变率必须为正数，当前值: {args.mutation_rate}")
    elif args.mutation_rate > 1e-5:
        result.add_warning(f"突变率过高 ({args.mutation_rate})，自然基因组中每碱基每代突变率通常在1e-8至1e-9范围内")
    elif args.mutation_rate < 1e-10:
        result.add_warning(f"突变率过低 ({args.mutation_rate})，可能导致几乎没有变异")
    
    # 检查注释文件格式
    if os.path.exists(args.annotations):
        try:
            with open(args.annotations, 'r') as f:
                # 跳过注释行
                line = f.readline()
                while line.startswith('#'):
                    line = f.readline()
                    if not line:
                        break
                
                if line:
                    fields = line.strip().split('\t')
                    if len(fields) < 5:
                        result.add_warning(f"注释文件格式可能不正确，每行应至少包含5列")
        except Exception as e:
            result.add_warning(f"读取注释文件时出错: {e}")
    
    return result

# 时间线进化参数验证
def validate_timeline_params(args):
    """
    验证时间线进化参数
    
    参数:
    args -- 参数命名空间，必须包含:
            base_genome, te_library, stages, te_percentages,
            evolution_levels, output_dir, mutation_rate
    
    返回:
    result -- ValidationResult对象
    """
    result = ValidationResult()
    
    # 验证输入文件存在
    if not os.path.exists(args.base_genome):
        result.add_error(f"基因组文件不存在: {args.base_genome}")
    if not os.path.exists(args.te_library):
        result.add_error(f"TE文库文件不存在: {args.te_library}")
    if hasattr(args, 'config') and args.config and not os.path.exists(args.config):
        result.add_warning(f"配置文件不存在: {args.config}，将使用默认参数")
    
    # 验证阶段数
    if args.stages <= 0:
        result.add_error(f"进化阶段数必须为正数，当前值: {args.stages}")
    elif args.stages > 10:
        result.add_warning(f"进化阶段数较多 ({args.stages})，可能导致运行时间过长")
    
    # 解析TE百分比和进化水平
    try:
        te_percentages = [float(p) for p in args.te_percentages.split(',')]
    except ValueError:
        result.add_error(f"无效的TE百分比列表: {args.te_percentages}，应为逗号分隔的数字")
        te_percentages = []
    
    evolution_levels = args.evolution_levels.split(',')
    
    # 验证TE百分比
    if len(te_percentages) < args.stages:
        result.add_warning(f"TE百分比数量 ({len(te_percentages)}) 少于阶段数 ({args.stages})，将使用最后一个值填充")
    
    for percentage in te_percentages:
        if percentage < 0:
            result.add_error(f"TE百分比不能为负数: {percentage}%")
    
    if te_percentages and sum(te_percentages[:min(len(te_percentages), args.stages)]) > 90:
        result.add_warning(f"总TE占比过高 ({sum(te_percentages[:min(len(te_percentages), args.stages)])}%)，自然基因组中TE占比通常不超过70%")
    
    # 验证进化水平
    valid_levels = ['low', 'medium', 'high']
    if len(evolution_levels) < args.stages:
        result.add_warning(f"进化水平数量 ({len(evolution_levels)}) 少于阶段数 ({args.stages})，将使用最后一个值填充")
    
    for level in evolution_levels:
        if level not in valid_levels:
            result.add_error(f"无效的进化水平: {level}，有效选项为 {', '.join(valid_levels)}")
    
    # 验证TE区域变异率乘数
    if hasattr(args, 'te_multiplier'):
        if args.te_multiplier <= 0:
            result.add_error(f"TE区域变异率乘数必须为正数，当前值: {args.te_multiplier}")
        elif args.te_multiplier > 10:
            result.add_warning(f"TE区域变异率乘数过高 ({args.te_multiplier})，自然情况下通常为1.5-3")
    
    # 验证突变率
    if args.mutation_rate <= 0:
        result.add_error(f"突变率必须为正数，当前值: {args.mutation_rate}")
    elif args.mutation_rate > 1e-5:
        result.add_warning(f"突变率过高 ({args.mutation_rate})，自然基因组中每碱基每代突变率通常在1e-8至1e-9范围内")
    elif args.mutation_rate < 1e-10:
        result.add_warning(f"突变率过低 ({args.mutation_rate})，可能导致几乎没有变异")
    
    return result

# TE相似度分析参数验证
def validate_similarity_params(args):
    """
    验证TE相似度分析参数
    
    参数:
    args -- 参数命名空间，必须包含:
            te_library, te_sequences, annotations(可选),
            output_prefix, mutation_rate
    
    返回:
    result -- ValidationResult对象
    """
    result = ValidationResult()
    
    # 验证输入文件存在
    if not os.path.exists(args.te_library):
        result.add_error(f"TE文库文件不存在: {args.te_library}")
    if not os.path.exists(args.te_sequences):
        result.add_error(f"TE序列文件不存在: {args.te_sequences}")
    if hasattr(args, 'annotations') and args.annotations and not os.path.exists(args.annotations):
        result.add_warning(f"注释文件不存在: {args.annotations}，将跳过部分分析")
    
    # 验证突变率
    if args.mutation_rate <= 0:
        result.add_error(f"突变率必须为正数，当前值: {args.mutation_rate}")
    elif args.mutation_rate > 1e-5:
        result.add_warning(f"突变率过高 ({args.mutation_rate})，自然基因组中每碱基每代突变率通常在1e-8至1e-9范围内")
    elif args.mutation_rate < 1e-10:
        result.add_warning(f"突变率过低 ({args.mutation_rate})，可能导致进化时间估算不准确")
    
    # 检查FASTA文件格式
    if os.path.exists(args.te_library):
        try:
            with open(args.te_library, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    result.add_error(f"TE文库文件不是有效的FASTA格式: {args.te_library}")
        except Exception as e:
            result.add_error(f"读取TE文库文件时出错: {e}")
    
    if os.path.exists(args.te_sequences):
        try:
            with open(args.te_sequences, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    result.add_error(f"TE序列文件不是有效的FASTA格式: {args.te_sequences}")
        except Exception as e:
            result.add_error(f"读取TE序列文件时出错: {e}")
    
    # 创建输出目录
    output_dir = os.path.dirname(os.path.abspath(args.output_prefix))
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except PermissionError:
            result.add_error(f"无权限创建输出目录: {output_dir}")
    
    return result

# 可视化参数验证
def validate_visualization_params(args):
    """
    验证可视化参数
    
    参数:
    args -- 参数命名空间，必须包含:
            genome, annotations, config(可选),
            output_prefix, format(可选)
    
    返回:
    result -- ValidationResult对象
    """
    result = ValidationResult()
    
    # 验证输入文件存在
    if not os.path.exists(args.genome):
        result.add_error(f"基因组文件不存在: {args.genome}")
    if not os.path.exists(args.annotations):
        result.add_error(f"注释文件不存在: {args.annotations}")
    if hasattr(args, 'config') and args.config and not os.path.exists(args.config):
        result.add_warning(f"配置文件不存在: {args.config}，将使用默认配置")
    
    # 验证输出格式
    if hasattr(args, 'format'):
        valid_formats = ['png', 'pdf', 'svg']
        if args.format not in valid_formats:
            result.add_error(f"无效的输出格式: {args.format}，有效选项为 {', '.join(valid_formats)}")
    
    # 创建输出目录
    output_dir = os.path.dirname(os.path.abspath(args.output_prefix))
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except PermissionError:
            result.add_error(f"无权限创建输出目录: {output_dir}")
    
    # 检查文件格式
    if os.path.exists(args.genome):
        try:
            with open(args.genome, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    result.add_error(f"基因组文件不是有效的FASTA格式: {args.genome}")
        except Exception as e:
            result.add_error(f"读取基因组文件时出错: {e}")
    
    # 检查注释文件格式（支持BED和GFF）
    if os.path.exists(args.annotations):
        try:
            with open(args.annotations, 'r') as f:
                # 跳过注释行
                line = f.readline()
                while line.startswith('#'):
                    line = f.readline()
                    if not line:
                        break
                
                if line:
                    fields = line.strip().split('\t')
                    # BED格式至少有5列，GFF格式应有9列
                    if len(fields) < 5 and len(fields) != 9:
                        result.add_warning(f"注释文件格式可能不正确，应为BED格式(至少5列)或GFF格式(9列)")
        except Exception as e:
            result.add_warning(f"读取注释文件时出错: {e}")
    
    return result

# 工作流总参数验证
def validate_workflow_params(args):
    """
    验证完整工作流参数
    
    参数:
    args -- 参数命名空间，必须包含多个子模块所需的参数
    
    返回:
    result -- ValidationResult对象
    """
    result = ValidationResult()
    
    # 验证TE文库
    if not os.path.exists(args.te_library):
        result.add_error(f"TE文库文件不存在: {args.te_library}")
    else:
        try:
            with open(args.te_library, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    result.add_error(f"TE文库文件不是有效的FASTA格式: {args.te_library}")
        except Exception as e:
            result.add_error(f"读取TE文库文件时出错: {e}")
    
    # 验证基因组长度
    if args.length <= 0:
        result.add_error(f"基因组长度必须为正数，当前值: {args.length}")
    elif args.length < 10000:
        result.add_warning(f"基因组长度过小 ({args.length} bp)，可能不足以模拟真实染色体结构")
    elif args.length > 1e9:
        result.add_warning(f"基因组长度过大 ({args.length:,} bp)，可能导致内存不足或运行时间过长")
    
    # 验证TE占比
    if args.percentage <= 0:
        result.add_error(f"TE占比必须为正数，当前值: {args.percentage}%")
    elif args.percentage > 90:
        result.add_warning(f"TE占比过高 ({args.percentage}%)，自然基因组中TE占比通常不超过70%")
    
    # 验证进化水平
    valid_levels = ['low', 'medium', 'high']
    if args.evolution not in valid_levels:
        result.add_error(f"无效的进化水平: {args.evolution}，有效选项为 {', '.join(valid_levels)}")
    
    # 验证阶段数
    if args.stages <= 0:
        result.add_error(f"进化阶段数必须为正数，当前值: {args.stages}")
    elif args.stages > 10:
        result.add_warning(f"进化阶段数较多 ({args.stages})，可能导致运行时间过长")
    
    # 验证突变率
    if args.mutation_rate <= 0:
        result.add_error(f"突变率必须为正数，当前值: {args.mutation_rate}")
    elif args.mutation_rate > 1e-5:
        result.add_warning(f"突变率过高 ({args.mutation_rate})，自然基因组中每碱基每代突变率通常在1e-8至1e-9范围内")
    elif args.mutation_rate < 1e-10:
        result.add_warning(f"突变率过低 ({args.mutation_rate})，可能导致几乎没有变异")
    
    # 验证GC含量参数
    if args.gc < 0 or args.gc > 1:
        result.add_error(f"GC含量必须在0-1范围内，当前值: {args.gc}")
    
    # 验证GC梯度参数
    if args.gc_gradient:
        if args.min_gc < 0 or args.min_gc > 1:
            result.add_error(f"最小GC含量必须在0-1范围内，当前值: {args.min_gc}")
        if args.max_gc < 0 or args.max_gc > 1:
            result.add_error(f"最大GC含量必须在0-1范围内，当前值: {args.max_gc}")
        if args.min_gc >= args.max_gc:
            result.add_error(f"最小GC含量 ({args.min_gc}) 必须小于最大GC含量 ({args.max_gc})")
    
    # 验证染色体区域大小
    if args.telomere_size < 0:
        result.add_error(f"端粒区域大小不能为负数，当前值: {args.telomere_size}")
    if args.centromere_size < 0:
        result.add_error(f"着丝粒区域大小不能为负数，当前值: {args.centromere_size}")
    
    # 验证区域大小与基因组长度的关系
    total_special_regions = 2 * args.telomere_size + args.centromere_size
    if total_special_regions >= args.length:
        result.add_error(f"特殊区域总大小 ({total_special_regions:,} bp) 不能大于等于基因组长度 ({args.length:,} bp)")
    elif total_special_regions > args.length * 0.5:
        result.add_warning(f"特殊区域占比过高 ({total_special_regions/args.length:.2%})，可能导致普通染色体臂区域过小")
    
    # 验证输出目录
    if not os.path.exists(args.output):
        try:
            os.makedirs(args.output, exist_ok=True)
        except PermissionError:
            result.add_error(f"无权限创建输出目录: {args.output}")
    
    return result

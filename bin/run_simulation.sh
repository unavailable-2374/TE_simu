#!/bin/bash
# run_simulation.sh - 执行增强版TE扩增模拟工作流

# 系统优化设置
# --------------------------------------------------
# 使用taskset绑定进程到特定CPU内核，提高缓存效率
if [ -x "$(command -v taskset)" ]; then
    TASKSET_CMD="taskset -c 0-47"  # 假设有48个核心可用
else
    TASKSET_CMD=""
fi

# 设置资源限制 - 增加文件打开数量限制
ulimit -n 65535 2>/dev/null || true

# 优化Python性能的环境变量
export PYTHONUNBUFFERED=1
export PYTHONHASHSEED=0
export PYTHONFAULTHANDLER=1
export PYTHONIOENCODING=utf-8

# 关闭NumPy/SciPy内部的多线程，防止与我们的多进程冲突
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export PYTHONWARNINGS="ignore"

# 获取脚本所在目录的绝对路径
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# --------------------------------------------------
# 默认参数
GENOME_LENGTH=50000000
TE_PERCENTAGE=30
EVOLUTION_LEVEL="medium"
STAGES=3
TE_LIBRARY=""
OUTPUT_DIR="te_simulation_output"
SEED=$(date +%s)  # 使用当前时间戳作为默认种子
GC_CONTENT=0.5
GC_GRADIENT=false
MIN_GC=0.35
MAX_GC=0.65
TELOMERE_SIZE=30000
CENTROMERE_SIZE=100000
MUTATION_RATE=1e-8
TE_DISTRIBUTION="60,30,10"
# 添加并行处理参数
NUM_PROCESSES=0  # 0表示自动确定（使用CPU核心数的75%）
# 添加快速模式参数
FAST_MODE=false  # 默认不使用快速模式

# 显示用法
usage() {
    echo "用法: $0 [选项]"
    echo "选项:"
    echo "  -l, --length GENOME_LENGTH    初始基因组长度 (默认: 50000000)"
    echo "  -p, --percentage TE_PERCENTAGE    目标TE占比 (默认: 30)"
    echo "  -e, --evolution EVOLUTION_LEVEL    进化水平 (low/medium/high, 默认: medium)"
    echo "  -s, --stages STAGES    进化阶段数量 (默认: 3)"
    echo "  -t, --te_library TE_LIBRARY    TE文库文件 (必需)"
    echo "  -o, --output OUTPUT_DIR    输出目录 (默认: te_simulation_output)"
    echo "  -d, --distribution TE_DISTRIBUTION    TE插入阶段分布 (逗号分隔的百分比, 默认: 60,30,10)"
    echo "  --seed SEED    随机数种子 (默认: 当前时间戳)"
    echo "  --gc GC_CONTENT    基础GC含量 (默认: 0.5)"
    echo "  --gc_gradient    启用GC含量梯度"
    echo "  --min_gc MIN_GC    GC梯度最小值 (默认: 0.35)"
    echo "  --max_gc MAX_GC    GC梯度最大值 (默认: 0.65)"
    echo "  --telomere_size TELOMERE_SIZE    端粒区域大小 (默认: 30000 bp)"
    echo "  --centromere_size CENTROMERE_SIZE    着丝粒区域大小 (默认: 100000 bp)"
    echo "  --mutation_rate MUTATION_RATE    每碱基每代的突变率 (默认: 1e-8)"
    echo "  --num_processes NUM_PROCESSES    并行处理的进程数量 (默认: 自动)"
    echo "  --fast                           使用快速TE插入算法 (不考虑区域偏好和GC含量)"
    echo "  -h, --help    显示此帮助信息"
    exit 1
}

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -l|--length)
            GENOME_LENGTH="$2"
            shift 2
            ;;
        -p|--percentage)
            TE_PERCENTAGE="$2"
            shift 2
            ;;
        -e|--evolution)
            EVOLUTION_LEVEL="$2"
            shift 2
            ;;
        -s|--stages)
            STAGES="$2"
            shift 2
            ;;
        -t|--te_library)
            TE_LIBRARY="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --seed)
            SEED="$2"
            shift 2
            ;;
        --gc)
            GC_CONTENT="$2"
            shift 2
            ;;
        --gc_gradient)
            GC_GRADIENT=true
            shift
            ;;
        --min_gc)
            MIN_GC="$2"
            shift 2
            ;;
        --max_gc)
            MAX_GC="$2"
            shift 2
            ;;
        --telomere_size)
            TELOMERE_SIZE="$2"
            shift 2
            ;;
        --centromere_size)
            CENTROMERE_SIZE="$2"
            shift 2
            ;;
        --mutation_rate)
            MUTATION_RATE="$2"
            shift 2
            ;;
        --num_processes)
            NUM_PROCESSES="$2"
            shift 2
            ;;
        --fast)
            FAST_MODE=true
            shift
            ;;
        -d|--distribution)
            TE_DISTRIBUTION="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "未知选项: $1"
            usage
            ;;
    esac
done

# 检查必需参数
if [ -z "$TE_LIBRARY" ]; then
    echo "错误: 必须指定TE文库文件 (-t, --te_library)"
    usage
fi

if [ ! -f "$TE_LIBRARY" ]; then
    echo "错误: TE文库文件 '$TE_LIBRARY' 不存在"
    exit 1
fi

# 参数验证（简化版）
if [ "$GENOME_LENGTH" -le 0 ]; then
    echo "错误: 基因组长度必须为正数，当前值: $GENOME_LENGTH"
    exit 1
fi

if [[ ! "$EVOLUTION_LEVEL" =~ ^(low|medium|high)$ ]]; then
    echo "错误: 无效的进化水平: $EVOLUTION_LEVEL，有效选项为 low, medium, high"
    exit 1
fi

if (( $(echo "$GC_CONTENT < 0 || $GC_CONTENT > 1" | bc -l) )); then
    echo "错误: GC含量必须在0-1范围内，当前值: $GC_CONTENT"
    exit 1
fi

IFS=',' read -r -a DIST_ARRAY <<< "$TE_DISTRIBUTION"

if [ ${#DIST_ARRAY[@]} -ne $STAGES ]; then
    echo "警告: TE分布比例数量 (${#DIST_ARRAY[@]}) 与阶段数 ($STAGES) 不匹配"
    
    # 如果提供的分布比例数量少于阶段数，使用最后一个值填充
    if [ ${#DIST_ARRAY[@]} -lt $STAGES ]; then
        LAST_VALUE=${DIST_ARRAY[${#DIST_ARRAY[@]}-1]}
        for ((i=${#DIST_ARRAY[@]}; i<$STAGES; i++)); do
            DIST_ARRAY[$i]=$LAST_VALUE
        done
    else
        # 如果提供的分布比例数量多于阶段数，截断多余的值
        DIST_ARRAY=("${DIST_ARRAY[@]:0:$STAGES}")
    fi
    
    # 重新构建TE_DISTRIBUTION字符串
    TE_DISTRIBUTION=$(IFS=,; echo "${DIST_ARRAY[*]}")
    echo "调整后的TE分布比例: $TE_DISTRIBUTION"
fi

TOTAL_DIST=0
for val in "${DIST_ARRAY[@]}"; do
    TOTAL_DIST=$(echo "$TOTAL_DIST + $val" | bc)
done

# 构建TE百分比序列，根据分布比例分配总百分比
TE_PERCENTAGES=""
for val in "${DIST_ARRAY[@]}"; do
    # 计算该阶段的实际TE百分比
    STAGE_PERCENTAGE=$(echo "scale=4; $TE_PERCENTAGE * $val / $TOTAL_DIST" | bc)

    if [ -z "$TE_PERCENTAGES" ]; then
        TE_PERCENTAGES="$STAGE_PERCENTAGE"
    else
        TE_PERCENTAGES="$TE_PERCENTAGES,$STAGE_PERCENTAGE"
    fi
done

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 自动设置进程数和分区数
if [ "$NUM_PROCESSES" -eq 0 ]; then
    NUM_PROCESSES=$(nproc)
    echo "自动检测到 $NUM_PROCESSES 个CPU核心"
fi

# 设置分区数 - 根据核心数动态调整
NUM_PARTITIONS=$((NUM_PROCESSES / 2))
if [ "$NUM_PARTITIONS" -lt 2 ]; then
    NUM_PARTITIONS=2
elif [ "$NUM_PARTITIONS" -gt 32 ]; then
    NUM_PARTITIONS=32  # 限制最大分区数
fi

# 准备并行处理参数
PROCESSES_PARAM="--num_processes $NUM_PROCESSES"

# 准备快速模式参数
FAST_MODE_PARAM=""
if [ "$FAST_MODE" = true ]; then
    FAST_MODE_PARAM="--fast_mode"
else
    FAST_MODE_PARAM="--use_weights"
fi

# 显示参数
echo "=== 增强版TE扩增模拟参数 ==="
echo "初始基因组长度: $GENOME_LENGTH bp"
echo "目标TE占比: $TE_PERCENTAGE%"
echo "进化水平: $EVOLUTION_LEVEL"
echo "进化阶段数: $STAGES"
echo "TE分布比例: $TE_DISTRIBUTION"
echo "每个阶段的TE插入百分比: $TE_PERCENTAGES"
echo "TE文库: $TE_LIBRARY"
echo "输出目录: $OUTPUT_DIR"
echo "随机数种子: $SEED"
echo "基础GC含量: $GC_CONTENT"
echo "GC含量梯度: $GC_GRADIENT"
if [ "$GC_GRADIENT" = true ]; then
    echo "GC梯度范围: $MIN_GC - $MAX_GC"
fi
echo "端粒区域大小: $TELOMERE_SIZE bp"
echo "着丝粒区域大小: $CENTROMERE_SIZE bp"
echo "突变率: $MUTATION_RATE"
echo "并行处理进程数: $NUM_PROCESSES"
echo "基因组分区数: $NUM_PARTITIONS"
echo "使用快速模式: $FAST_MODE"
echo

# 保存参数到文件
PARAMS_FILE="$OUTPUT_DIR/simulation_parameters.txt"
echo "# 增强版TE扩增模拟参数" > "$PARAMS_FILE"
echo "初始基因组长度: $GENOME_LENGTH bp" >> "$PARAMS_FILE"
echo "目标TE占比: $TE_PERCENTAGE%" >> "$PARAMS_FILE"
echo "进化水平: $EVOLUTION_LEVEL" >> "$PARAMS_FILE"
echo "进化阶段数: $STAGES" >> "$PARAMS_FILE"
echo "TE分布比例: $TE_DISTRIBUTION"
echo "每个阶段的TE插入百分比: $TE_PERCENTAGES"
echo "TE文库: $TE_LIBRARY" >> "$PARAMS_FILE"
echo "输出目录: $OUTPUT_DIR" >> "$PARAMS_FILE"
echo "随机数种子: $SEED" >> "$PARAMS_FILE"
echo "基础GC含量: $GC_CONTENT" >> "$PARAMS_FILE"
echo "GC含量梯度: $GC_GRADIENT" >> "$PARAMS_FILE"
if [ "$GC_GRADIENT" = true ]; then
    echo "GC梯度范围: $MIN_GC - $MAX_GC" >> "$PARAMS_FILE"
fi
echo "端粒区域大小: $TELOMERE_SIZE bp" >> "$PARAMS_FILE"
echo "着丝粒区域大小: $CENTROMERE_SIZE bp" >> "$PARAMS_FILE"
echo "突变率: $MUTATION_RATE" >> "$PARAMS_FILE"
echo "并行处理进程数: $NUM_PROCESSES" >> "$PARAMS_FILE"
echo "基因组分区数: $NUM_PARTITIONS" >> "$PARAMS_FILE"
echo "快速模式: $FAST_MODE" >> "$PARAMS_FILE"
echo "运行时间: $(date)" >> "$PARAMS_FILE"

# 步骤1: 生成结构化初始随机基因组
echo "=== 步骤1: 生成结构化初始随机基因组 ==="
BASE_GENOME="$OUTPUT_DIR/base_genome.fa"

# 构建基因组生成命令
GENOME_CMD="$TASKSET_CMD python $SCRIPT_DIR/generate_base_genome.py -l $GENOME_LENGTH -o $BASE_GENOME --gc $GC_CONTENT --seed $SEED"
GENOME_CMD="$GENOME_CMD --telomere_size $TELOMERE_SIZE --centromere_size $CENTROMERE_SIZE $PROCESSES_PARAM"

if [ "$GC_GRADIENT" = true ]; then
    GENOME_CMD="$GENOME_CMD --gc_gradient --min_gc $MIN_GC --max_gc $MAX_GC"
fi

# 执行命令
eval $GENOME_CMD
echo

# 根据阶段数调整TE百分比
STAGE_PERCENTAGE=$(echo "scale=2; $TE_PERCENTAGE / $STAGES" | bc)
echo "每个阶段TE插入百分比: $STAGE_PERCENTAGE%"

# 构建进化水平序列 - 简化版本
EVOLUTION_LEVELS=""
if [ "$EVOLUTION_LEVEL" = "low" ]; then
    # 所有阶段都使用低水平
    PATTERN=("low")
elif [ "$EVOLUTION_LEVEL" = "medium" ]; then
    # 逐步增加: 低 -> 中 -> 高 -> 循环
    PATTERN=("low" "medium" "high")
elif [ "$EVOLUTION_LEVEL" = "high" ]; then
    # 开始为中等，其余为高水平
    PATTERN=("medium" "high")
fi

# 构建进化水平字符串
for ((i=1; i<=$STAGES; i++)); do
    # 确定模式数组中的索引（循环使用）
    idx=$(( (i-1) % ${#PATTERN[@]} ))
    level=${PATTERN[$idx]}
    
    if [ -z "$EVOLUTION_LEVELS" ]; then
        EVOLUTION_LEVELS="$level"
    else
        EVOLUTION_LEVELS="$EVOLUTION_LEVELS,$level"
    fi
done

echo "进化水平序列: $EVOLUTION_LEVELS"

# 构建TE百分比序列
TE_PERCENTAGES=""
for ((i=1; i<=$STAGES; i++)); do
    if [ -z "$TE_PERCENTAGES" ]; then
        TE_PERCENTAGES="$STAGE_PERCENTAGE"
    else
        TE_PERCENTAGES="$TE_PERCENTAGES,$STAGE_PERCENTAGE"
    fi
done

# 获取基因组配置文件
BASE_CONFIG="${BASE_GENOME%.*}_config.json"

# 步骤2: 运行增强版时间线演化模拟
echo "=== 步骤2: 运行增强版时间线演化模拟 ==="

$TASKSET_CMD python "$SCRIPT_DIR/genomic_partition.py" \
    --genome "$BASE_GENOME" \
    --te_library "$TE_LIBRARY" \
    --stages "$STAGES" \
    --te_percentages "$TE_PERCENTAGES" \
    --evolution_levels "$EVOLUTION_LEVELS" \
    --output_dir "$OUTPUT_DIR/timeline" \
    --config "$BASE_CONFIG" \
    --mutation_rate "$MUTATION_RATE" \
    --seed "$SEED" \
    --num_partitions "$NUM_PARTITIONS" \
    --te_multiplier 2.0 \
    --fast_mode "$FAST_MODE"

# 步骤3: 生成最终注释和统计
echo "=== 步骤3: 生成最终注释和统计 ==="
FINAL_GENOME="$OUTPUT_DIR/timeline/final_genome.fa"
FINAL_ANNOTATIONS="$OUTPUT_DIR/timeline/final_annotations.bed"
FINAL_TE_SEQUENCES="$OUTPUT_DIR/timeline/final_te_sequences.fa"
FINAL_OUTPUT="$OUTPUT_DIR/final"

$TASKSET_CMD python "$SCRIPT_DIR/generate_annotation.py" \
    -g "$FINAL_GENOME" \
    -a "$FINAL_ANNOTATIONS" \
    -o "$FINAL_OUTPUT" \
    $PROCESSES_PARAM

# 步骤4: 分析TE相似度和进化时间
echo "=== 步骤4: 分析TE相似度和进化时间 ==="
$TASKSET_CMD python "$SCRIPT_DIR/te_similarity.py" \
    -t "$TE_LIBRARY" \
    -s "$FINAL_TE_SEQUENCES" \
    -a "$FINAL_ANNOTATIONS" \
    -o "$FINAL_OUTPUT/te_similarity" \
    -r "$MUTATION_RATE" \
    $PROCESSES_PARAM

echo
echo "=== 模拟完成 ==="
echo "最终文件位于: $OUTPUT_DIR"
echo "统计报告: $FINAL_OUTPUT""_te_statistics.txt"
echo "相似度报告: $FINAL_OUTPUT/te_similarity_report.txt"
echo "可视化脚本: $FINAL_OUTPUT""_visualization.R"
echo
echo "要生成可视化图表，请运行:"
echo "Rscript $FINAL_OUTPUT""_visualization.R"
echo
echo "要生成TE相似度和进化时间图表，请运行:"
echo "Rscript $FINAL_OUTPUT/te_similarity_visualization.R"
echo
echo "完整的模拟参数已保存至: $PARAMS_FILE"

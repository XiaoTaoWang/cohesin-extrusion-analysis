#!/bin/bash
#SBATCH -A wangxiaotao
#SBATCH -p cu,fat,gpu
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --job-name=quant
# 这里的 %A 是主作业 ID，%a 是阵列索引
#SBATCH --array=0-8
#SBATCH --output=quant.%A.out
#SBATCH --error=quant.%A.err

# 1. 定义需要处理的任务目录列表
# 您可以将目录名放入一个数组中
dirs=(
    "loading_50pct"
    "loading_25pct"
    "loading_10pct"
    "loading_0pct"
    "ctcf_50pct"
    "ctcf_25pct"
    "ctcf_10pct"
    "ctcf_0pct"
    "wild"
)

# 获取当前任务对应的目录名
target_dir=${dirs[$SLURM_ARRAY_TASK_ID]}

echo "$target_dir"

# 4. 运行 Python 脚本
#python quant.py "$target_dir" "$target_dir"
python quant3.py "$target_dir" "$target_dir"

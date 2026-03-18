#!/bin/bash
#SBATCH -A wangxiaotao
#SBATCH -p cu,gpu,fat
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=3
#SBATCH --mem=3G
#SBATCH --job-name=agg
# 这里的 %A 是主作业 ID，%a 是阵列索引
#SBATCH --array=0-8
#SBATCH --output=logs/agg.%A_%a.out
#SBATCH --error=logs/agg.%A_%a.err

# 1. 定义需要处理的任务目录列表
# 您可以将目录名放入一个数组中
dirs=(
    "loading_0pct"
    "loading_10pct"
    "loading_25pct"
    "loading_50pct"
    "ctcf_0pct"
    "ctcf_10pct"
    "ctcf_25pct"
    "ctcf_50pct"
    "wild"
)

# 获取当前任务对应的目录名
target_dir=${dirs[$SLURM_ARRAY_TASK_ID]}
input_dir="${target_dir}/cool_outputs"

python pileup.py "$input_dir" "$target_dir"
python plot_agg4.py "$target_dir"

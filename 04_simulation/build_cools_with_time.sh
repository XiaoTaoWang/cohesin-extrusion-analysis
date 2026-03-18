#!/bin/bash
#SBATCH -A wangxiaotao
#SBATCH -p cu,fat,gpu
#SBATCH -t 24:00:00        # 缩短时间限制，800小时通常会被排队系统降权
#SBATCH -N 1               # 保持单节点，因为该任务不是跨节点分布式
#SBATCH -n 1               # 单个任务
#SBATCH --cpus-per-task=12  # 增加 CPU 核心数，用于加速构建和潜在的 joblib 并行
#SBATCH --mem=32G          # 增加内存。读取大型 HDF5 和构建接触图（1000x1000）需要一定内存空间
#SBATCH --job-name=build_cools
#SBATCH --output=job.%j.out
#SBATCH --error=job.%j.err

#python build_cools.py loading_0pct loading_0pct/cool_outputs 
#python build_cools.py loading_50pct loading_50pct/cool_outputs 
#python build_cools_with_time.py wild wild/cool_with_time
#python build_cools_with_time.py speed1 speed1/cool_with_time
#python build_cools_with_time.py speed0.5 speed0.5/cool_with_time
#python build_cools_with_time.py control control/cool_with_time
#python build_cools_with_time.py control control/cool_with_time_u199
#python build_cools_with_time.py wild wild/cool_with_time_u199
#python build_cools_with_time.py control control/cool_with_time_u39
#python build_cools_with_time.py wild wild/cool_with_time_u39
#python build_cools_with_time.py wild_speed1/ wild_speed1/cool_with_time_u39 39
#python build_cools_with_time.py wild_speed1/ wild_speed1/cool_with_time_u19 19
#python build_cools_with_time.py wild_speed1/ wild_speed1/cool_with_time_u199 199
#python build_cools_with_time.py wild_speed1/ wild_speed1/cool_with_time_u999 999
#python build_cools_with_time.py wild_stability0.01/ wild_stability0.01/cool_with_time_u299 299
#python build_cools_with_time.py control_speed0.5 control_speed0.5/cool_with_time_u299 299
#python build_cools_with_time.py wild_stability0.01/ wild_stability0.01/cool_with_time_u899 899
#python build_cools_with_time.py control_speed0.5 control_speed0.5/cool_with_time_u899 899

#python build_cools_with_time.py wild_speed1 wild_speed1/cool_with_time_u899 899
#for up in 19 39 199 299; 
for up in 999; 
do 
    python build_cools_with_time.py wild_speed1 wild_speed1/cool_with_time_u$up $up
done

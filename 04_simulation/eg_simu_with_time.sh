#!/bin/bash
#SBATCH -A wangxiaotao
#SBATCH -p cu,fat,gpu
#SBATCH -t 24:00:00        # 缩短时间限制，800小时通常会被排队系统降权
#SBATCH -N 1               # 保持单节点，因为该任务不是跨节点分布式
#SBATCH -n 1               # 单个任务
#SBATCH --cpus-per-task=1  # 增加 CPU 核心数，用于加速构建和潜在的 joblib 并行
#SBATCH --mem=3G          # 增加内存。读取大型 HDF5 和构建接触图（1000x1000）需要一定内存空间
#SBATCH --job-name=eg_simu_with_time
#SBATCH --output=job.%j.out
#SBATCH --error=job.%j.err

#python eg_simu_with_time.py control/cool_with_time_u199
#python eg_simu_with_time.zoom.py control/cool_with_time_u199
#python eg_simu_with_time.py wild/cool_with_time_u199
#python eg_simu_with_time.zoom.py wild/cool_with_time_u199
#python eg_simu_with_time.py wild/cool_with_time_u39
#python eg_simu_with_time.zoom.py wild/cool_with_time_u39
#python eg_simu_with_time.py control/cool_with_time_u39
#python eg_simu_with_time.zoom.py control/cool_with_time_u39
#python eg_simu_with_time.zoom.py wild_speed1/cool_with_time_u19
#python eg_simu_with_time.zoom.py wild_speed1/cool_with_time_u39
#python eg_simu_with_time.zoom.py wild_speed1/cool_with_time_u199
#python eg_simu_with_time.zoom.py wild_speed1/cool_with_time_u999

python eg_simu_with_time.py wild_speed1/cool_with_time_u19
#python eg_simu_with_time.py wild_speed1/cool_with_time_u39
#python eg_simu_with_time.py wild_speed1/cool_with_time_u199
#python eg_simu_with_time.py wild_speed1/cool_with_time_u999

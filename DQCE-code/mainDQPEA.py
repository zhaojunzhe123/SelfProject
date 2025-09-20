# coding:utf-8
import copy
import os
from pathlib import Path  # 引入pathlib处理路径更方便

from DataRead import DataReadDHFJSP
import numpy as np
from Initial import GHInitial
from inital import initial
from fitFJSP import CalfitDHFJFP
from Tselection import *
from EA import evolution
from Tool import *
from FastNDSort import FastNDS
from EnergySave import EnergysavingDHFJSP
from LocalSearch import *
from DQN_model import DQN
import torch

# 问题实例组合
Combination=[[10,2],[20,2],[30,2],[40,2],
             [20,3],[30,3],[40,3],[50,3],
             [40,4],[50,4],[100,4],
             [50,5],[100,5],[150,5],
             [100,6],[150,6],[200,6],
             [100,7],[150,7],[200,7]]

# 数据路径设置
datapath = '/home/zhaojunzhe/Project/code/selfProject/DQCE-code/DATASET/'
FileName = []
ResultPath = []

for i in range(20):
    J = Combination[i][0]
    F1 = Combination[i][1]
    O = 5
    temp = f"{datapath}{J}J{F1}F.txt"  # 使用f-string更简洁
    temp2 = f"{J}J{F1}F"
    FileName.append(temp)
    ResultPath.append(temp2)

TF = 20
# FileName = np.array(FileName).reshape(TF, 1)
# ResultPath = np.array(ResultPath).reshape(TF, 1)

# 读取算法参数
param_path = "/home/zhaojunzhe/Project/code/selfProject/DQCE-code/parameter.txt"
with open(param_path, "r", encoding='utf-8') as f:  # 使用with语句更安全
    params = f.read().split(' ')
    ps = int(params[0])
    Pc = float(params[1])
    Pm = float(params[2])
    lr = float(params[3])
    batch_size = int(params[4])
    EPSILON = float(params[5])
    GAMMA = float(params[6])
    MEMORY_CAPACITY = int(params[7])

IndependentRun = 10

# DQN参数（覆盖从文件读取的值）
lr = 0.001
batch_size = 16
EPSILON = 0.9               # 贪婪策略
GAMMA = 0.9                 # 奖励折扣
MEMORY_CAPACITY = 512
TARGET_REPLACE_ITER = 7     # 目标网络更新频率
N_ACTIONS = 9               # 可选算子数量
EPOCH = 1


print("CUDA可用状态:", torch.cuda.is_available())

# 执行算法
CCF = 19
for file in range(CCF, CCF + 1):
    # 读取数据
    N, F, TM, H, SH, NM, M, time, ProF = DataReadDHFJSP(FileName[file])
    MaxNFEs = 20 * SH
    
    # 路径设置（使用pathlib处理路径）
    current_dir = Path.cwd()  # 获取当前工作目录
    result_dir = current_dir / "DQNV9+ES" / ResultPath[file][0]  # 构建结果目录路径
    
    # 创建结果目录（parents=True表示创建父目录，exist_ok=True避免已存在时报错）
    result_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"{ResultPath[file]} 正在优化中...\n")
    
    # 独立运行
    for rround in range(1):
        # 初始化种群
        p_chrom, m_chrom, f_chrom = initial(N, H, SH, NM, M, ps, F)
        fitness = np.zeros(shape=(ps, 3))
        NFEs = 0  # 函数评估次数
        
        # 计算初始适应度
        for i in range(ps):
            fitness[i, 0], fitness[i, 1], fitness[i, 2] = CalfitDHFJFP(
                p_chrom[i, :], m_chrom[i, :], f_chrom[i, :], 
                N, H, SH, F, TM, time
            )
        
        # 精英存档
        AP = []
        AM = []
        AF = []
        AFit = []
        iter_count = 1  # 迭代计数器
        max_iterations = 20
        # 构建DQN模型
        N_STATES = 2 * SH + N
        CountOpers = np.zeros(N_ACTIONS)
        PopCountOpers = []
        dq_net = DQN(
            N_STATES, N_ACTIONS, 
            BATCH_SIZE=batch_size, 
            LR=lr, 
            EPSILON=EPSILON, 
            GAMMA=GAMMA, 
            MEMORY_CAPACITY=MEMORY_CAPACITY, 
            TARGET_REPLACE_ITER=TARGET_REPLACE_ITER
        )
        Loss = []
        
        # 主循环
        while NFEs < MaxNFEs and iter_count < max_iterations:
            print(f"{FileName[file]} 轮次 {rround + 1} 迭代 {iter_count}")
            iter_count += 1
            
            # 初始化子代
            ChildP = np.zeros(shape=(2 * ps, SH), dtype=int)
            ChildM = np.zeros(shape=(2 * ps, SH), dtype=int)
            ChildF = np.zeros(shape=(2 * ps, N), dtype=int)
            ChildFit = np.zeros(shape=(2 * ps, 3))
            
            # 选择操作
            P_pool, M_pool, F_pool = tournamentSelection(
                p_chrom, m_chrom, f_chrom, fitness, ps, SH, N
            )
            
            # 生成子代
            for j in range(ps):
                Fit1 = np.zeros(3)
                Fit2 = np.zeros(3)
                P1, M1, F1, P2, M2, F2 = evolution(
                    P_pool, M_pool, F_pool, j, Pc, Pm, ps, SH, N, H, NM, M
                )
                
                # 计算适应度
                Fit1[0], Fit1[1], Fit1[2] = CalfitDHFJFP(
                    P1, M1, F1, N, H, SH, F, TM, time
                )
                Fit2[0], Fit2[1], Fit2[2] = CalfitDHFJFP(
                    P2, M2, F2, N, H, SH, F, TM, time
                )
                NFEs += 2
                
                # 保存子代
                t1 = j * 2
                t2 = j * 2 + 1
                ChildP[t1, :] = copy.copy(P1)
                ChildM[t1, :] = copy.copy(M1)
                ChildF[t1, :] = copy.copy(F1)
                ChildFit[t1, :] = Fit1
                
                ChildP[t2, :] = copy.copy(P2)
                ChildM[t2, :] = copy.copy(M2)
                ChildF[t2, :] = copy.copy(F2)
                ChildFit[t2, :] = Fit2
            
            # 合并父代和子代
            QP = np.vstack((p_chrom, ChildP))
            QM = np.vstack((m_chrom, ChildM))
            QF = np.vstack((f_chrom, ChildF))
            QFit = np.vstack((fitness, ChildFit))
            
            # 去重
            QP, QM, QF, QFit = DeleteReapt(QP, QM, QF, QFit, ps)
            
            # 快速非支配排序
            RQFit = QFit[:, 0:2]
            TopRank = FastNDS(RQFit, ps)
            p_chrom = QP[TopRank, :]
            m_chrom = QM[TopRank, :]
            f_chrom = QF[TopRank, :]
            fitness = QFit[TopRank, :]
            
            # 精英策略
            PF = pareto(fitness)
            if len(AFit) == 0:
                AP = copy.copy(p_chrom[PF, :])
                AM = copy.copy(m_chrom[PF, :])
                AF = copy.copy(f_chrom[PF, :])
                AFit = copy.copy(fitness[PF, :])
            else:
                AP = np.vstack((AP, p_chrom[PF, :]))
                AM = np.vstack((AM, m_chrom[PF, :]))
                AF = np.vstack((AF, f_chrom[PF, :]))
                AFit = np.vstack((AFit, fitness[PF, :]))
            
            # 对精英存档进行非支配排序和去重
            PF = pareto(AFit)
            AP = AP[PF, :]
            AM = AM[PF, :]
            AF = AF[PF, :]
            AFit = AFit[PF, :]
            AP, AM, AF, AFit = DeleteReaptE(AP, AM, AF, AFit)
            
            # 存档中的局部搜索
            L = len(AFit)
            current_state = np.zeros(N_STATES, dtype=int)
            next_state = np.zeros(N_STATES, dtype=int)
            
            for l in range(L):
                # 构建当前状态
                current_state[0:SH] = copy.copy(AP[l, :])
                current_state[SH:SH*2] = copy.copy(AM[l, :])
                current_state[SH*2:N_STATES] = copy.copy(AF[l, :])
                
                # 选择动作
                action = dq_net.choose_action(current_state)
                k = int(action)
                
                # 根据动作选择不同的局部搜索算子
                if k == 0:
                    P1, M1, F1 = N6(
                        AP[l, :], AM[l, :], AF[l, :], AFit[l, :], 
                        N, H, SH, time, TM, NM, M, F
                    )
                elif k == 1:
                    P1, M1, F1 = SwapOF(
                        AP[l, :], AM[l, :], AF[l, :], AFit[l, :], 
                        N, H, SH, time
                    )
                elif k == 2:
                    P1, M1, F1 = RandFA(
                        AP[l, :], AM[l, :], AF[l, :], AFit[l, :], 
                        N, H, SH, time, TM, NM, M, F
                    )
                elif k == 3:
                    P1, M1, F1 = RandMS(
                        AP[l, :], AM[l, :], AF[l, :], AFit[l, :], 
                        N, H, SH, time, TM, NM, M, F
                    )
                elif k == 4:
                    P1, M1, F1 = InsertOF(
                        AP[l, :], AM[l, :], AF[l, :], AFit[l, :], 
                        N, H, SH, time
                    )
                elif k == 5:
                    P1, M1, F1 = InsertIF(
                        AP[l, :], AM[l, :], AF[l, :], AFit[l, :], 
                        N, H, SH, time, F
                    )
                elif k == 6:
                    P1, M1, F1 = SwapIF(
                        AP[l, :], AM[l, :], AF[l, :], AFit[l, :], 
                        N, H, SH, time, F
                    )
                elif k == 7:
                    P1, M1, F1 = RankFA(
                        AP[l, :], AM[l, :], AF[l, :], AFit[l, :], 
                        N, H, SH, time, TM, NM, M, F, ProF
                    )
                elif k == 8:
                    P1, M1, F1 = RankMS(
                        AP[l, :], AM[l, :], AF[l, :], AFit[l, :], 
                        N, H, SH, time, TM, NM, M, F
                    )
                
                # 计算新解的适应度
                Fit1[0], Fit1[1], Fit1[2] = CalfitDHFJFP(
                    P1, M1, F1, N, H, SH, F, TM, time
                )
                NFEs += 1
                
                # 计算支配关系并给予奖励
                dom = NDS(Fit1, AFit[l, :])
                if dom == 1:
                    # 新解支配旧解
                    AP[l, :] = copy.copy(P1)
                    AM[l, :] = copy.copy(M1)
                    AF[l, :] = copy.copy(F1)
                    AFit[l, :] = copy.copy(Fit1)
                    AP = np.vstack((AP, P1))
                    AM = np.vstack((AM, M1))
                    AF = np.vstack((AF, F1))
                    AFit = np.vstack((AFit, Fit1))
                    reward = 5
                elif dom == 0 and AFit[l][0] != Fit1[0] and AFit[l][1] != Fit1[1]:
                    # 新解与旧解非支配
                    AP = np.vstack((AP, P1))
                    AM = np.vstack((AM, M1))
                    AF = np.vstack((AF, F1))
                    AFit = np.vstack((AFit, Fit1))
                    reward = 10
                else:
                    # 新解被支配
                    reward = 0
                
                # 构建下一状态并存储经验
                next_state[0:SH] = copy.copy(P1)
                next_state[SH:SH*2] = copy.copy(M1)
                next_state[SH*2:N_STATES] = copy.copy(F1)
                dq_net.store_transition(current_state, action, reward, next_state)
                
                # 学习
                if dq_net.memory_counter > 50:
                    for epoch in range(EPOCH):
                        loss = dq_net.learn()
                        Loss.append(loss)
            
            # 节能优化
            L = len(AFit)
            for j in range(L):
                P1, M1, F1 = EnergysavingDHFJSP(
                    AP[j, :], AM[j, :], AF[j, :], AFit[j, :], 
                    N, H, TM, time, SH, F
                )
                Fit1[0], Fit1[1], Fit1[2] = CalfitDHFJFP(
                    P1, M1, F1, N, H, SH, F, TM, time
                )
                NFEs += 1
                
                if NDS(Fit1, AFit[j, :]) == 1:
                    AP[j, :] = copy.copy(P1)
                    AM[j, :] = copy.copy(M1)
                    AF[j, :] = copy.copy(F1)
                    AFit[j, :] = copy.copy(Fit1)
                    AP = np.vstack((AP, P1))
                    AM = np.vstack((AM, M1))
                    AF = np.vstack((AF, F1))
                    AFit = np.vstack((AFit, Fit1))
                elif NDS(Fit1, AFit[j, :]) == 0:
                    AP = np.vstack((AP, P1))
                    AM = np.vstack((AM, M1))
                    AF = np.vstack((AF, F1))
                    AFit = np.vstack((AFit, Fit1))
        
        # 保存精英解
        PF = pareto(AFit)
        AP = AP[PF, :]
        AM = AM[PF, :]
        AF = AF[PF, :]
        AFit = AFit[PF, :]
        PF = pareto(AFit)
        l = len(PF)
        
        # 提取目标值并去重
        obj = AFit[:, 0:2]
        newobj = []
        for i in range(l):
            newobj.append(obj[PF[i], :])
        newobj = np.unique(newobj, axis=0)
        
        # 保存结果文件（使用pathlib处理路径）
        # 保存结果文件（使用pathlib处理路径）
        res_file = result_dir / f"res{rround + 1}.txt"
        with open(res_file, "w", encoding='utf-8') as f:
            for item in newobj:
                f.write(f"{item[0]:5.2f} {item[1]:6.2f}\n")
        print(f"结果已保存到: {res_file}")  # 显示结果文件路径

        # 保存损失文件
        loss_file = current_dir / "loss.txt"
        with open(loss_file, "w", encoding='utf-8') as f:
            for loss in Loss:
                f.write(f"{loss:.6f}\n")
        print(f"损失已保存到: {loss_file}")  # 显示损失文件路径
        
    print(f"完成 {FileName[file]}")
print("所有运行完成")

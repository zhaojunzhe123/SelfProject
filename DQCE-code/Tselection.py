# coding:utf-8
import numpy as np
import math
import random
from Tool import NDS
def tournamentSelection(p_chrom,m_chrom,f_chrom,fitness,ps,SH,N):
    #mating selection pool
    pool_size=ps
    P_pool = np.zeros(shape=(ps, SH), dtype=int)
    M_pool = np.zeros(shape=(ps, SH), dtype=int)
    F_pool = np.zeros(shape=(ps, N), dtype=int) #fitness of pool solutions
    # compeitor number
    tour=2

    # tournament selection
    for i in range(pool_size):
        index1=int(math.floor(random.random()*ps))
        index2 = int(math.floor(random.random() * ps))
        while index1==index2:
            index2 = int(math.floor(random.random() * ps))
        f1=fitness[index1,0:2]
        f2=fitness[index2,0:2]
        if (NDS(f1, f2) == 1):
            P_pool[i,:]=p_chrom[index1,:]
            M_pool[i,:]=m_chrom[index1,:]
            F_pool[i,:]=f_chrom[index1,:]
        elif(NDS(f1, f2) == 2):
            P_pool[i, :] = p_chrom[index2, :]
            M_pool[i, :] = m_chrom[index2, :]
            F_pool[i, :] = f_chrom[index2, :]
        else:
            if random.random() <= 0.5:
                P_pool[i, :] = p_chrom[index1, :]
                M_pool[i, :] = m_chrom[index1, :]
                F_pool[i, :] = f_chrom[index1, :]
            else:
                P_pool[i, :] = p_chrom[index2, :]
                M_pool[i, :] = m_chrom[index2, :]
                F_pool[i, :] = f_chrom[index2, :]
    return P_pool,M_pool,F_pool
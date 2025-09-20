# coding:utf-8
import copy
import random
import Tool
import numpy as np
import math
from Tool import*

def crossover(P1,M1,F1,P2,M2,F2,N,SH):
    #inital offerspring
    NP1=P1;NM1=M1
    NP2=P2;NM2=M2
    NF1=F1;NF2=F2
    #index of each operation in P1 and P2
    ci1=np.zeros(SH,dtype=int)
    ci2 = np.zeros(SH, dtype=int)
    # store some jobs in J1
    temp=[random.random() for _ in range(N) ]
    temp=Tool.mylistRound(temp)
    J1=find_all_index(temp,1)#find the index where value equal to 1
    for j in range(SH):
        if Ismemeber(P1[j], J1)==1: #if is in job set J
            ci1[j] = P1[j]+1

        if Ismemeber(P2[j], J1)==0: #if is not in job set J
            ci2[j] = P2[j]+1
    index_1_1 = find_all_index(ci1,0) # find the empty positions in ci1
    index_1_2 = find_all_index_not(ci2,0) # find the positions in ci2 which is not zero

    index_2_1 = find_all_index(ci2,0)
    index_2_2 = find_all_index_not(ci1,0)
    l1=len(index_1_1);l2=len(index_2_1)
    for j in range(l1):
        ci1[index_1_1[j]] = NP2[index_1_2[j]]
    for j in range(l2):
        ci2[index_2_1[j]] = NP1[index_2_2[j]]
    l1 = len(index_2_2);l2 = len(index_1_2)
    for j in range(l1):
        ci1[index_2_2[j]] = ci1[index_2_2[j]]-1
    for j in range(l2):
        ci2[index_1_2[j]] = ci2[index_1_2[j]] - 1
    NP1=ci1
    NP2 = ci2

    # applied uninversal crossover for machine selection
    s = [random.random() for _ in range(SH)]
    s = Tool.mylistRound(s)
    for i in range(0,SH):
        if (s[i] == 1):
            t = NM1[i]
            NM1[i] = NM2[i];
            NM2[i] = t;
    # applied uninversal crossover for machine selection
    s = [random.random() for _ in range(N)]
    s = Tool.mylistRound(s)
    for i in range(0, N):
        if (s[i] == 1):
            t = NF1[i]
            NF1[i] = NF2[i];
            NF2[i] = t;
    return NP1,NM1,NF1,NP2,NM2,NF2

def mutation(p_chrom,m_chrom,SH,N,H,NM,M):
    #swap for operation sequence as mutation operator
    p1=math.floor(random.random()*SH)
    p2 = math.floor(random.random() * SH)
    while p1==p2:
        p2 = math.floor(random.random() * SH)
    t = p_chrom[p1]
    p_chrom[p1] = p_chrom[p2]
    p_chrom[p2] = t;

    #change a machine for machine selection as mutation operator
    s1 = p_chrom
    s2 = np.zeros(SH, dtype=int)
    p = np.zeros(N, dtype=int)
    fitness = np.zeros(2)
    for i in range(SH):
        p[s1[i]] = p[s1[i]] + 1
        s2[i] = p[s1[i]]
    s3=m_chrom
    p1 = math.floor(random.random() * SH)
    p2 = math.floor(random.random() * SH)
    while p1 == p2:
        p2 = math.floor(random.random() * SH)
    #print('Job= ',s1[p1],'Op= ',s2[p1]-1,'MaxOp=',H[s1[p1]])
    n = NM[s1[p1]][s2[p1]-1]
    m = math.floor(random.random() * n)
    x = M[s1[p1]][s2[p1]-1][m]-1
    if n>1:
        while s3[p1]==x:
            m = math.floor(random.random() * n)
            x = M[s1[p1]][s2[p1]-1][m]-1
    k1 = 0
    for t2 in range(s1[p1]):#sum from 0 to s1[p1]
        k1 = k1 + H[t2]
    t1 = int(k1 + s2[p1]-1)
    m_chrom[t1] = x

    n = NM[s1[p2]][s2[p2] - 1]
    m = math.floor(random.random() * n)
    x = M[s1[p2]][s2[p2] - 1][m]-1
    if n > 1:
        while s3[p2] == x:
            m = math.floor(random.random() * n)
            x = M[s1[p2]][s2[p2] - 1][m]-1
    k1 = 0
    for t2 in range(s1[p2]):#sum from 0 to s1[p2]
        k1 = k1 + H[t2]
    t1 = int(k1 + s2[p2] - 1)
    m_chrom[t1] = x
    return p_chrom,m_chrom

def evolution(p_chrom,m_chrom,f_chrom,index,Pc,Pm,ps,SH,N,H,NM,M):
    R = math.floor(random.random() * ps)
    P1=copy.copy(p_chrom[index,:])
    P2 = copy.copy(p_chrom[R,:])
    M1=copy.copy(m_chrom[index,:])
    M2 = copy.copy(m_chrom[R,:])
    F1=copy.copy(f_chrom[index,:])
    F2 = copy.copy(f_chrom[R,:])
    while R == index:
        R = math.floor(random.random() * ps)
    if random.random()<Pc:
        P1,M1,F1,P2,M2,F2=crossover(p_chrom[index,:],m_chrom[index,:],f_chrom[index,:],p_chrom[R,:],m_chrom[R,:],f_chrom[R,:],N,SH)
    if random.random()<Pm:
        P1,M1=mutation(P1,M1,SH,N,H,NM,M)
        P2, M2 = mutation(P2, M2, SH, N, H, NM, M)
    return P1,M1,F1,P2,M2,F2

def evolution2(p_chrom,m_chrom,f_chrom,index,T,neighbour,Pc,Pm,ps,SH,N,H,NM,M):
    nei=neighbour[index,:]

    R1 = math.floor(random.random() * T)
    R1=nei[R1]
    R2 = math.floor(random.random() * T)
    R1 = nei[R2]

    P1=copy.copy(p_chrom[R1,:])
    P2 = copy.copy(p_chrom[R2,:])
    M1=copy.copy(m_chrom[R1,:])
    M2 = copy.copy(m_chrom[R2,:])
    F1=copy.copy(f_chrom[R1,:])
    F2 = copy.copy(f_chrom[R2,:])
    while R1 == R2:
        R2 = math.floor(random.random() * T)
        R2= nei[R2]
    if random.random()<Pc:
        P1,M1,F1,P2,M2,F2=crossover(p_chrom[R1,:],m_chrom[R1,:],f_chrom[R1,:],p_chrom[R2,:],m_chrom[R2,:],f_chrom[R2,:],N,SH)
    if random.random()<Pm:
        P1,M1=mutation(P1,M1,SH,N,H,NM,M)
        P2, M2 = mutation(P2, M2, SH, N, H, NM, M)
    return P1,M1,F1,P2,M2,F2
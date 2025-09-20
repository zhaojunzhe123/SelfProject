# coding:utf-8
#problem specific local search
import copy
import random

import numpy as np
from CriticalPath import FindCriticalPathDHFJSP

#Swap out of factory
def SwapOF(p_chrom,m_chrom,f_chrom,fitness,N,H,SH,time):
    Index1=np.floor(random.random()*SH)
    Index2 = np.floor(random.random() * SH)
    while Index1==Index2:
        Index2 = np.floor(random.random() * SH)
    newp=p_chrom
    newm=m_chrom
    newf=f_chrom
    Index1=int(Index1);Index2=int(Index1)
    tmp = copy.copy(newp[Index1]);
    newp[Index1] = copy.copy(newp[Index2])# 交换两点的染色体值
    newp[Index2] = copy.copy(tmp);
    return newp,newm,newf

#Swap in critical factory
def SwapIF(p_chrom,m_chrom,f_chrom,fitness,N,H,SH,time,F):
    s1 = p_chrom
    s2 = np.zeros(SH, dtype=int)
    p = np.zeros(N, dtype=int)
    for i in range(SH):
        p[s1[i]] = p[s1[i]] + 1
        s2[i] = p[s1[i]]
    P0 = [];
    P = [];
    IP = []
    FJ = []
    for f in range(F):
        P.append([])
        IP.append([])
        FJ.append([])

    for i in range(SH):
        t1 = s1[i]
        t2 = s2[i]
        t3 = f_chrom[t1]
        P[t3].append(p_chrom[i])
        IP[t3].append(i)
    for i in range(N):
        t3 = f_chrom[i]
        FJ[t3].append(i)

    cf = int(fitness[2])
    L=len(IP[cf])
    Index1=np.floor(random.random()*L)
    Index2 = np.floor(random.random() * L)
    while Index1==Index2:
        Index2 = np.floor(random.random() * L)
    newp=p_chrom
    newm=m_chrom
    newf=f_chrom
    Index1=int(Index1);Index2=int(Index1)
    Index1=IP[cf][Index1];Index2=IP[cf][Index2]
    tmp = copy.copy(newp[Index1]);
    newp[Index1] = copy.copy(newp[Index2])# 交换两点的染色体值
    newp[Index2] = copy.copy(tmp);
    return newp,newm,newf

#Insert out of factory
def InsertOF(p_chrom,m_chrom,f_chrom,fitness,N,H,SH,time):
    Index1 = np.floor(random.random() * SH)
    Index2 = np.floor(random.random() * SH)
    while Index1 == Index2:
        Index2 = np.floor(random.random() * SH)
    newp = p_chrom
    newm = m_chrom
    newf = f_chrom
    Index1 = int(Index1);
    Index2 = int(Index1)
    low=min(Index1,Index2)
    up=max(Index1,Index2)
    tmp = newp[up];
    for i in range(up,low,-1):
        newp[i]=copy.copy(newp[i-1])
    newp[low]=copy.copy(tmp)
    return newp, newm, newf

#Insert out of factory
def InsertIF(p_chrom,m_chrom,f_chrom,fitness,N,H,SH,time,F):
    s1 = p_chrom
    s2 = np.zeros(SH, dtype=int)
    p = np.zeros(N, dtype=int)
    for i in range(SH):
        p[s1[i]] = p[s1[i]] + 1
        s2[i] = p[s1[i]]
    P0 = [];
    P = [];
    IP = []
    FJ = []
    for f in range(F):
        P.append([])
        IP.append([])
        FJ.append([])

    for i in range(SH):
        t1 = s1[i]
        t2 = s2[i]
        t3 = f_chrom[t1]
        P[t3].append(p_chrom[i])
        IP[t3].append(i)
    for i in range(N):
        t3 = f_chrom[i]
        FJ[t3].append(i)

    cf = int(fitness[2])
    L = len(IP[cf])
    Index1 = np.floor(random.random() * L)
    Index2 = np.floor(random.random() * L)
    while Index1 == Index2:
        Index2 = np.floor(random.random() * L)

    newp = p_chrom
    newm = m_chrom
    newf = f_chrom
    Index1 = int(Index1);Index2 = int(Index1)
    Index1 = IP[cf][Index1];Index2 = IP[cf][Index2]
    low=min(Index1,Index2)
    up=max(Index1,Index2)
    tmp = newp[up];
    for i in range(up,low,-1):
        newp[i]=copy.copy(newp[i-1])
    newp[low]=copy.copy(tmp)
    return newp, newm, newf

#find a critical operation and randomly change a job to another factory
def RandFA(p_chrom,m_chrom,f_chrom,fitness,N,H,SH,time,TM,NM,M,F):
    s1 = p_chrom
    s2 = np.zeros(SH, dtype=int)
    p = np.zeros(N, dtype=int)
    for i in range(SH):
        p[s1[i]] = p[s1[i]] + 1
        s2[i] = p[s1[i]]
    P0 = [];
    P = [];
    FJ = []
    for f in range(F):
        P.append([])
        FJ.append([])

    for i in range(SH):
        t1 = s1[i]
        t2 = s2[i]
        t3 = f_chrom[t1]
        P[t3].append(p_chrom[i])
    for i in range(N):
        t3 = f_chrom[i]
        FJ[t3].append(i)

    cf=int(fitness[2])
    CP,_,_=FindCriticalPathDHFJSP(P[cf],m_chrom,FJ[cf],cf,N,H,time,TM)
    L=len(CP)
    IndexO=CP[int(np.floor(random.random()*L))]
    s3=copy.copy(P[cf])
    I2=s3[IndexO]
    of=int(np.floor(random.random()*F))
    while of==cf:
        of = int(np.floor(random.random() * F))
    newf=copy.copy(f_chrom)
    newf[I2]=of
    return p_chrom, m_chrom, newf

#find a critical operation and randomly change a job to another factory according to its total processing time / machine
#the lower total processing time the bigger probability of the factory will be selected
#select method is roulette
def RankFA(p_chrom,m_chrom,f_chrom,fitness,N,H,SH,time,TM,NM,M,F,ProF):
    s1 = p_chrom
    s2 = np.zeros(SH, dtype=int)
    p = np.zeros(N, dtype=int)
    for i in range(SH):
        p[s1[i]] = p[s1[i]] + 1
        s2[i] = p[s1[i]]
    P0 = [];
    P = [];
    FJ = []
    for f in range(F):
        P.append([])
        FJ.append([])

    for i in range(SH):
        t1 = s1[i]
        t2 = s2[i]
        t3 = f_chrom[t1]
        P[t3].append(p_chrom[i])
    for i in range(N):
        t3 = f_chrom[i]
        FJ[t3].append(i)

    cf=int(fitness[2])
    CP,_,_=FindCriticalPathDHFJSP(P[cf],m_chrom,FJ[cf],cf,N,H,time,TM)
    L=len(CP)
    IndexO=CP[int(np.floor(random.random()*L))]
    s3=copy.copy(P[cf])
    I2=s3[IndexO]

    pro=ProF[I2,:]
    pro_index=np.argsort(pro)
    for f in range(1,F):
        pro[pro_index[f]]=pro[pro_index[f]]+pro[pro_index[f-1]]
    x=random.random()
    for f in range(F):
        if x<pro[pro_index[f]]:
            of=pro_index[f]
            break

    of=int(np.floor(random.random()*F))
    while of==cf:
        x = random.random()
        for f in range(F):
            if x < pro[pro_index[f]]:
                of = pro_index[f]
                break
    newf=copy.copy(f_chrom)
    newf[I2]=int(of)
    return p_chrom, m_chrom, newf

#Randomly select a operation and assiged to another machine
def RandMS(p_chrom,m_chrom,f_chrom,fitness,N,H,SH,time,TM,NM,M,F):
    s1 = p_chrom
    s2 = np.zeros(SH, dtype=int)
    p = np.zeros(N, dtype=int)
    for i in range(SH):
        p[s1[i]] = p[s1[i]] + 1
        s2[i] = p[s1[i]]
    P0 = [];
    P = [];
    FJ = []
    for f in range(F):
        P.append([])
        FJ.append([])

    for i in range(SH):
        t1 = s1[i]
        t2 = s2[i]
        t3 = f_chrom[t1]
        P[t3].append(p_chrom[i])
    for i in range(N):
        t3 = f_chrom[i]
        FJ[t3].append(i)

    cf = int(fitness[2])
    CP, _, _ = FindCriticalPathDHFJSP(P[cf], m_chrom, FJ[cf], cf, N, H, time, TM)
    L = len(CP)
    IndexO = CP[int(np.floor(random.random() * L))]
    s3 = copy.copy(P[cf])
    L = len(s3)
    s4=np.zeros(L)
    p2=np.zeros(N)
    for i in range(L):
        p2[s3[i]] = p2[s3[i]] + 1
        s4[i] = p2[s3[i]]
    I2=s3[IndexO]
    J2=s4[IndexO]
    for i in range(SH):
        if s1[i]==I2 and s2[i]==J2:
            IndexO=i
            break
    I=s1[IndexO];J=s2[IndexO]
    newm=m_chrom
    t4 = 0;t1=I;t2=J;
    for k in range(t1):  # sum from 0 to t1-1
        t4 = t4 + H[k]
    tmp = m_chrom[t4 + t2 - 1]
    n=NM[I][J-1]
    cm=np.zeros(n)
    for kk in range(n):
        cm[kk]=M[I][J-1][kk]-1
    Index2=cm[int(np.floor(random.random()*n))]
    newm[t4 + t2 - 1]=int(Index2)
    newp=p_chrom
    return p_chrom, m_chrom, f_chrom

#Randomly select a operation and assiged to another machine
def RankMS(p_chrom,m_chrom,f_chrom,fitness,N,H,SH,time,TM,NM,M,F):
    s1 = p_chrom
    s2 = np.zeros(SH, dtype=int)
    p = np.zeros(N, dtype=int)
    for i in range(SH):
        p[s1[i]] = p[s1[i]] + 1
        s2[i] = p[s1[i]]
    P0 = [];
    P = [];
    FJ = []
    for f in range(F):
        P.append([])
        FJ.append([])

    for i in range(SH):
        t1 = s1[i]
        t2 = s2[i]
        t3 = f_chrom[t1]
        P[t3].append(p_chrom[i])
    for i in range(N):
        t3 = f_chrom[i]
        FJ[t3].append(i)

    cf = int(fitness[2])
    CP, _, _ = FindCriticalPathDHFJSP(P[cf], m_chrom, FJ[cf], cf, N, H, time, TM)
    L = len(CP)
    IndexO = CP[int(np.floor(random.random() * L))]
    s3 = copy.copy(P[cf])
    L = len(s3)
    s4=np.zeros(L)
    p2=np.zeros(N)
    for i in range(L):
        p2[s3[i]] = p2[s3[i]] + 1
        s4[i] = p2[s3[i]]
    I2=s3[IndexO]
    J2=s4[IndexO]
    for i in range(SH):
        if s1[i]==I2 and s2[i]==J2:
            IndexO=i
            break
    I=s1[IndexO];J=s2[IndexO]
    newm=m_chrom
    t4 = 0;t1=I;t2=J;
    for k in range(t1):  # sum from 0 to t1-1
        t4 = t4 + H[k]
    tmp = m_chrom[t4 + t2 - 1]
    n=NM[I][J-1]
    cm=np.zeros(n,dtype=int)
    cmt=np.zeros(n)
    tot=0
    for kk in range(n):
        cm[kk]=M[I][J-1][kk]-1
        cmt[kk]=time[cf][I][J-1][cm[kk]]
        tot=cmt[kk]+tot
    for kk in range(n):
        cmt[kk]=cmt[kk]/tot
    pro_index=np.argsort(cmt)
    pro = copy.copy(cmt[pro_index])
    for f in range(1,n):
        pro[f]=pro[f]+pro[f-1]
    x=random.random()
    for f in range(n):
        if x<pro[f]:
            Index2=cm[pro_index[f]]
            break

    while tmp==Index2:
        x = random.random()
        for f in range(n):
            if x < pro[f]:
                Index2 = cm[pro_index[f]]
                break
    newm[t4 + t2 - 1]=int(Index2)
    newp=p_chrom
    return p_chrom, m_chrom, f_chrom

def N6(p_chrom,m_chrom,f_chrom,fitness,N,H,SH,time,TM,NM,M,F):
    s1 = p_chrom
    s2 = np.zeros(SH, dtype=int)
    p = np.zeros(N, dtype=int)
    for i in range(SH):
        p[s1[i]] = p[s1[i]] + 1
        s2[i] = p[s1[i]]
    P0 = [];
    P = [];IP=[]
    FJ = []
    for f in range(F):
        P.append([])
        IP.append([])
        FJ.append([])

    for i in range(SH):
        t1 = s1[i]
        t2 = s2[i]
        t3 = f_chrom[t1]
        P[t3].append(p_chrom[i])
        IP[t3].append(i)
    for i in range(N):
        t3 = f_chrom[i]
        FJ[t3].append(i)

    cf = int(fitness[2])
    CP, CB, block = FindCriticalPathDHFJSP(P[cf], m_chrom, FJ[cf], cf, N, H, time, TM)
    for i in range(block):
        BL=len(CB[i])
        if BL>1:
            if i==0:
                Index1=int(np.floor(random.random()*(BL-1)))
                Index2=BL-1
                Index1=CB[i][Index1];Index2=CB[i][Index2]
                tmp=P[cf][Index1]
                for j in range(Index1,Index2):
                    P[cf][j]=P[cf][j+1]
                P[cf][Index2]=tmp
            if i==block-1:
                Index1=0
                Index2=int(np.floor(random.random()*(BL-1))+1)
                Index1 = CB[i][Index1];Index2 = CB[i][Index2]
                tmp = P[cf][Index2]
                for j in range(Index2, Index1,-1):
                    P[cf][j] = P[cf][j-1]
                P[cf][Index1] = tmp
            if i>0 and i<block-1 and BL>2:
                Index1 = int(np.floor(random.random() * (BL - 2)) + 1)
                Index2=BL-1
                Index1 = CB[i][Index1];Index2 = CB[i][Index2]
                tmp = P[cf][Index1]
                for j in range(Index1, Index2):
                    P[cf][j] = P[cf][j + 1]
                P[cf][Index2] = tmp
                Index1 = 0
                Index2 = int(np.floor(random.random() * (BL - 2)) + 1)
                Index1 = CB[i][Index1];Index2 = CB[i][Index2]
                tmp = P[cf][Index2]
                for j in range(Index2, Index1, -1):
                    P[cf][j] = P[cf][j - 1]
                P[cf][Index1] = tmp
    newm=m_chrom;newf=f_chrom
    newp=np.zeros(SH,dtype=int)
    for f in range(F):
        L=len(IP[f])
        for i in range(L):
            newp[IP[f][i]]=P[f][i]
    return newp,newm,newf

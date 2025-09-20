# coding:utf-8
import numpy as np

def CalfitFJFP(p_chrom,m_chrom,FJ,f_index,N,H,TM,time):
    #processing power and idle power
    SH=len(p_chrom)
    Ep=4;Es=1
    opmax=int(max(H))
    #initial finish time matrix to record the finish time of each operation
    finish=np.zeros(shape=(N,opmax))
    #machine finish time
    mt=np.zeros(TM)
    #sign all operation
    s1=p_chrom
    s2=np.zeros(SH,dtype=int)
    p=np.zeros(N,dtype=int)
    fitness=np.zeros(2)
    for i in range(SH):
        p[s1[i]]=p[s1[i]]+1
        s2[i]=p[s1[i]]
    #assign the machine to each operation from machine selection vector
    mm=np.zeros(SH,dtype=int)
    for i in range(SH):
        t1=s1[i]
        t2=s2[i]
        t4=0
        for k in range(t1):#sum from 0 to t1-1
            t4=t4+H[k]
        mm[i]=m_chrom[t4+t2-1]
    #initial total workload and total idle time
    TWL=0;TIT=0;
    #start decoding
    for i in range(SH):
        # because array index starts with 0 in python, number in M start from 1,
        # number in time start from 0, thus number of m_chrom need to miner 1
        t1=s1[i];t2=s2[i]-1;t3=mm[i]
        if s2[i]==0:
            mt[t3]=mt[t3]+time[f_index][t1][t2][t3]
            if time[f_index][t1][t2][t3]==0:
                print('processing time is zero')
            finish[t1][t2]=mt[t3]
            TWL=TWL+time[f_index][t1][t2][t3]
        else:
            if mt[t3]<finish[t1][t2-1]:
                TIT=TIT+finish[t1][t2-1]-mt[t3]
                mt[t3]=finish[t1][t2-1]+time[f_index][t1][t2][t3]
                if time[f_index][t1][t2][t3] == 0:
                    print('processing time is zero')
                finish[t1][t2]=mt[t3]
                TWL = TWL + time[f_index][t1][t2][t3]
            else:
                mt[t3] = mt[t3] + time[f_index][t1][t2][t3]
                finish[t1][t2] = mt[t3]
                TWL = TWL + time[f_index][t1][t2][t3]
    fitness[0]=mt[0]
    for i in range(1,TM):
        if mt[i]>fitness[0]:
            fitness[0]=mt[i]
    fitness[1]=TWL*Ep+TIT*Es
    #print(fitness[0],fitness[1],TWL,TIT)
    return fitness[0],fitness[1]

def CalfitDHFJFP(p_chrom,m_chrom,f_chrom,N,H,SH,F,TM,time):
    s1 = p_chrom
    s2 = np.zeros(SH, dtype=int)
    p = np.zeros(N, dtype=int)
    fitness = np.zeros(2)
    for i in range(SH):
        p[s1[i]] = p[s1[i]] + 1
        s2[i] = p[s1[i]]
    P0=[];P=[];FJ=[]
    for f in range(F):
        P.append([])
        FJ.append([])

    for i in range(SH):
        t1=s1[i]
        t2=s2[i]
        t3=f_chrom[t1]
        P[t3].append(p_chrom[i])
    for i in range(N):
        t3=f_chrom[i]
        FJ[t3].append(i)
    sub_f_fit=np.zeros(shape=(F,2))

    for f in range(F):
        sub_f_fit[f][0],sub_f_fit[f][1]=CalfitFJFP(P[f],m_chrom,FJ[f],f,N,H,TM,time)

    fit1=sub_f_fit[0][0]
    fit3=1;fit2=0
    for f in range(F):
        fit2=sub_f_fit[f][1]+fit2
        if fit1<sub_f_fit[f][0]:
            fit1=sub_f_fit[f][0]
            fit3=f
    return fit1,fit2,fit3
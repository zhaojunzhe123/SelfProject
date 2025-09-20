# coding:utf-8
import copy

import numpy as np
import random
import math

def initial(N,H,SH,NM,M,ps,F):
    #create operation sequence and machine selection vectors
    p_chrom=np.zeros(shape=(ps,SH),dtype=int)
    m_chrom = np.zeros(shape=(ps, SH),dtype=int)
    f_chrom = np.zeros(shape=(ps, N), dtype=int)

    chrom=np.zeros(SH,dtype=int)
    FC = np.zeros(N, dtype=int)
    #generate operation sequence randomly
    for i in range(N):
        for j in range(int(H[i])):
            k=0
            for t in range(i-1):
                k=k+H[t]
            k=int(k+j)
            chrom[k]=i
        FC[i]=i%F

    tmp=chrom;tmp2=FC
    random.shuffle(tmp)
    random.shuffle(tmp2)
    p_chrom[0,:]=tmp
    f_chrom[0,:]=tmp2

    for i in range(1,ps):
        tmp=p_chrom[i-1,:]
        random.shuffle(tmp)
        p_chrom[i,:]=tmp
        tmp2 = f_chrom[i - 1, :]
        random.shuffle(tmp2)
        f_chrom[i, :] = tmp2
    #finish generate operation sequencing sizeing ps

    #start generating machine selection vector
    for k in range(ps):
        for i in range(N):
            for j in range(int(H[i])):
                #randomly select a machine from each operation's candidate set
                t=int(math.floor(random.random()*NM[i][j]))
                k1 = 0
                for t2 in range(i):
                    k1 = k1 + H[t2]
                t1=int(k1+j)
                m_chrom[k,t1]=M[i][j][t]-1#adjust the index from 0 to TM-1 to suit time
    return p_chrom,m_chrom,f_chrom

def RFA(subps,N,F):
    f_chrom = np.zeros(shape=(subps, N), dtype=int)
    FC = np.zeros(N, dtype=int)
    # generate operation sequence randomly
    for i in range(N):
        FC[i] = i % F
    tmp2 = FC
    random.shuffle(tmp2)
    f_chrom[0, :] = tmp2
    for i in range(1,subps):
        tmp2 = f_chrom[i - 1, :]
        random.shuffle(tmp2)
        f_chrom[i, :] = tmp2
    return f_chrom

def RMS(subps,N,H,SH,NM,M):
    m_chrom = np.zeros(shape=(subps, SH), dtype=int)
    # start generating machine selection vector
    for k in range(subps):
        for i in range(N):
            for j in range(int(H[i])):
                # randomly select a machine from each operation's candidate set
                t = int(math.floor(random.random() * NM[i][j]))
                k1 = 0
                for t2 in range(i):
                    k1 = k1 + H[t2]
                t1 = int(k1 + j)
                m_chrom[k, t1] = M[i][j][t] - 1  # adjust the index from 0 to TM-1 to suit time
    return m_chrom

def MinPTF(subps,N,F,ProF):
    f_chrom = np.zeros(shape=(subps, N), dtype=int)
    for k in range(subps):
        countf = np.zeros(F, dtype=int)
        Fchrom = np.zeros(N, dtype=int)
        for i in range(N):
            pro = ProF[i,:];
            pro_index= np.argsort(pro);
            jpf = np.mean(countf);
            for f in range(F):
                if countf[pro_index[f]] <= jpf:
                    of = f;
                    break;
            Fchrom[i]= of;
            countf[of] = countf[of] + 1;
        f_chrom[k,:]=copy.copy(Fchrom);
    return f_chrom

def MinPTM(f_chrom,subps,N,H,SH,NM,time,TM,M):
    m_chrom=np.zeros((subps,SH),dtype=int)
    for k in range(subps):
        for i in range(N):
            curf=f_chrom[k][i]
            for j in range(int(H[i])):
                t1 = i
                t2 = j
                t4 = 0
                for kk in range(t1):  # sum from 0 to t1-1
                    t4 = t4 + H[kk]
                mp = t4 + t2 - 1
                index=M[i][j][0]
                f = time[curf][i][j][index]
                NM1 = int(NM[i][j])
                for t in range(NM1):
                    d = M[i][j][t]-1
                    try:
                        time[curf][i][j][d]
                    except:
                        print(i,j,d)
                    if f > time[curf][i][j][d]:
                        f = time[curf][i][j][d]
                        index = d
                m_chrom[k][mp]=index
    return m_chrom

def MinFTM(p_chrom,f_chrom,subps,N,H,SH,NM,time,TM,M,F):
    opmax = int(max(H))
    m_chrom = np.zeros((subps, SH), dtype=int)
    for k in range(subps):
        s1 = copy.copy(p_chrom[k,:])
        s2 = np.zeros(SH, dtype=int)
        p = np.zeros(N, dtype=int)
        fitness = np.zeros(2)
        for i in range(SH):
            p[s1[i]] = p[s1[i]] + 1
            s2[i] = p[s1[i]]
        # assign the machine to each operation from machine selection vector
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
            t3 = f_chrom[k][t1]
            P[t3].append(p_chrom[k][i])

        for f in range(F):
            mt=np.zeros(TM)
            SH1=len(P[f])
            mm=np.zeros(SH,dtype=int)
            s3=copy.copy(P[f])
            s4 = np.zeros(SH1, dtype=int)
            p = np.zeros(N, dtype=int)
            finish = np.zeros(shape=(N, opmax))
            for i in range(SH1):
                p[s3[i]] = p[s3[i]] + 1
                s4[i] = p[s3[i]]
            for i in range(SH1):
                t1=s3[i]
                t2=s4[i]-1
                if s4[i]==0:
                    MachineIndex = FindMinFinishTimeMachine(t1, t2, mt,TM,N,M);
                    MinPT=FindMinProcessTimeMachine(t1, t2, MachineIndex, f, time, M);
                    try:
                        mm[i] =MinPT
                    except:
                        mm[i] =MinPT[0]

                    mt[mm[i]] = mt[mm[i]]+time[f][t1][t2][mm[i]]
                    finish[t1][t2] =mt[mm[i]]
                else:
                    MachineIndex = FindMinFinishTimeMachine(t1, t2, mt, TM, N, M,NM);
                    MinPT = FindMinProcessTimeMachine(t1, t2, MachineIndex, f, time, M);
                    try:
                        mm[i] = MinPT
                    except:
                        mm[i] =MinPT[0]
                    if mt[mm[i]]<finish[t1][t2-1]:
                        mt[mm[i]]=finish[t1][t2-1]+time[f][t1][t2][mm[i]]
                    else:
                        mt[mm[i]] = mt[mm[i]] + time[f][t1][t2][mm[i]]
                        finish[t1][t2]=mt[mm[i]]
                t4=0
                for kk in range(t1):  # sum from 0 to t1-1
                    t4 = t4 + H[kk]
                mp = t4 + t2
                m_chrom[k][mp]=mm[i]
    return m_chrom

def MinWLM(p_chrom,f_chrom,subps,N,H,SH,NM,time,TM,M,F):
    m_chrom = np.zeros((subps, SH), dtype=int)
    for k in range(subps):
        s1 = copy.copy(p_chrom[k, :])
        s2 = np.zeros(SH, dtype=int)
        p = np.zeros(N, dtype=int)
        fitness = np.zeros(2)
        for i in range(SH):
            p[s1[i]] = p[s1[i]] + 1
            s2[i] = p[s1[i]]
        # assign the machine to each operation from machine selection vector
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
            t3 = f_chrom[k][t1]
            P[t3].append(p_chrom[k][i])

        for f in range(F):
            mt = np.zeros(TM)
            SH1 = len(P[f])
            mm = np.zeros(SH, dtype=int)
            s3 = copy.copy(P[f])
            s4 = np.zeros(SH1, dtype=int)
            p = np.zeros(N, dtype=int)
            for i in range(SH1):
                p[s3[i]] = p[s3[i]] + 1
                s4[i] = p[s3[i]]
            for i in range(SH1):
                t1 = s3[i]
                t2 = s4[i]-1
                n=NM[t1][t2]
                if n == 1:
                    mm[i] = M[t1][t2][0]-1;
                    mt[mm[i]] = mt[mm[i]] + time[f][t1][t2][mm[i]]
                    continue;
                else:
                    avalible_m = np.zeros(n,dtype=int);
                    avalible_m_load = np.zeros(n,dtype=int);
                    for j in range(n):
                        avalible_m[j] =M[t1][t2][j]-1;
                        avalible_m_load[j]= mt[avalible_m[j]]
                    candidateM = selectMachine(avalible_m_load);
                    sizeM = len(candidateM);
                    if (sizeM == 1):
                        mm[i] = avalible_m[candidateM];
                        mt[mm[i]] = mt[mm[i]]+time[f][t1][t2][mm[i]]
                    else:
                        t = time[f][t1][t2][avalible_m[candidateM[0]]]
                        mm[i] = avalible_m[candidateM[0]]
                        for kk in range(1,sizeM):
                            tmp = time[f][t1][t2][avalible_m[candidateM[kk]]]
                            if (t > tmp):
                                mm[i] = avalible_m[candidateM[kk]]
                        mt[mm[i]] = mt[mm[i]] + time[f][t1][t2][mm[i]]
            for i in range(SH1):
                t1 = s3[i];
                t2 = s4[i];
                t4=0
                for kk in range(t1):  # sum from 0 to t1-1
                    t4 = t4 + H[kk]
                mp = t4 + t2-1
                m_chrom[k][mp]=mm[i]
    return m_chrom

def selectMachine(mt):#输入机器负载向量 返回最小的机器负载 索引 可能是大于等于1
    L=len(mt);
    candidateM=[];
    f=mt[0];
    index=1;
    for i in range(1,L):
        if (f>mt[i]):
            f=mt[i];
            index=i;
    candidateM.append(index)
    f=mt[index];
    for i in range(L):
        if(i==index):
            continue;
        elif(f==mt[i]):
            candidateM.append(i)
    return candidateM

def FindMinFinishTimeMachine(JobIndex, OperationIndex, mt,TM,N,M,NM):
    L=NM[JobIndex][OperationIndex]
    CandidateM = np.zeros(L,dtype=int);
    for i in range(L):
        CandidateM[i] = M[JobIndex][OperationIndex][i]-1;
    if L<2:
        ProtenialMachine = CandidateM
        return ProtenialMachine
    ProtenialMachine=[]
    MinFinishMIndex = CandidateM[0];
    ProtenialMachine.append(MinFinishMIndex)
    for j in range(1,L):
        if mt[MinFinishMIndex] == mt[CandidateM[j]]:
            ProtenialMachine.append(CandidateM[j])
        else:
            if mt[MinFinishMIndex] > mt[CandidateM[j]]:
                ProtenialMachine = [];
                MinFinishMIndex = CandidateM[j];
                ProtenialMachine.append(MinFinishMIndex)
    return ProtenialMachine

def FindMinProcessTimeMachine(JobIndex,OperationIndex,MachineIndex,f_index,time,M):
    z = 0;
    L = len(MachineIndex);
    if L < 2:
        MinProessTimeMachine = MachineIndex;
        return MinProessTimeMachine
    MinProessTimeMachine = MachineIndex[0];
    MinTime = time[f_index][JobIndex][OperationIndex][MinProessTimeMachine];
    ML=np.size(MinTime);
    if ML==0:
        MinTime = z;
    for i in range(1,L):
        temp = time[f_index][JobIndex][OperationIndex][MachineIndex[i]]
        ML = np.size(temp);
        if ML == 0:
            temp = z;
        if ((MinTime > temp) and (temp != z)):
            MinProessTimeMachine = copy.copy(MachineIndex[i])
            MinTime = time[f_index][JobIndex][OperationIndex][MinProessTimeMachine]

    return MinProessTimeMachine

def GHInitial(N,H,SH,NM,M,TM,time,F,ProF,ps):
    p_chrom = np.zeros(shape=(ps, SH), dtype=int)
    m_chrom = np.zeros(shape=(ps, SH), dtype=int)
    f_chrom = np.zeros(shape=(ps, N), dtype=int)

    chrom = np.zeros(SH, dtype=int)
    # generate operation sequence randomly
    for i in range(N):
        for j in range(int(H[i])):
            k = 0
            for t in range(i - 1):
                k = k + H[t]
            k = int(k + j)
            chrom[k] = i

    tmp = chrom;
    random.shuffle(tmp)
    p_chrom[0, :] = copy.copy(tmp)

    for i in range(1, ps):
        tmp = p_chrom[i - 1, :]
        random.shuffle(tmp)
        p_chrom[i, :] = copy.copy(tmp)

    subps=math.floor(ps/5)
    for i in range(5):
        low = (i) * subps;
        up = low + subps;
        subpchrom = p_chrom[low:up,:];
        print(i)
        if i==0:
            subfchrom = MinPTF(subps,N,F,ProF);
            submchrom = RMS(subps,N,H,SH,NM,M);
        elif i==1:
            subfchrom = RFA(subps,N,F);
            submchrom = MinPTM(subfchrom,subps,N,H,SH,NM,time,TM,M);
        elif i==2:
            subfchrom = RFA(subps,N,F);
            submchrom = MinWLM(subpchrom,subfchrom,subps,N,H,SH,NM,time,TM,M,F);
        elif i==3:
            submchrom = MinFTM(subpchrom,subfchrom,subps,N,H,SH,NM,time,TM,M,F);
            subfchrom = RFA(subps,N,F);
        elif i==4:
            submchrom = RMS(subps,N,H,SH,NM,M);
            subfchrom = RFA(subps,N,F);
        m_chrom[low: up,:]=copy.copy(submchrom);
        f_chrom[low: up,:]=copy.copy(subfchrom);
    return p_chrom,m_chrom,f_chrom
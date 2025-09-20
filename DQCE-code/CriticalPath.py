# coding:utf-8
#Find Critical Path and Block
import numpy as np

def CalculateLastOperation(s1,s2,mm,digraph,last,VELTime,Vflag,time,f_index):
    Index=last
    e=0
    t1=e;t2=e
    if digraph[Index][2]!=0:
        zlast=digraph[Index][2]
        if Vflag[zlast][0]==0:
            VELTime[zlast][0],VELTime,Vflag=CalculateLastOperation(s1,s2,mm,digraph,zlast,VELTime,Vflag,time,f_index)
            Vflag[zlast][0]=1
        t1=time[f_index][s1[zlast]][s2[zlast]-1][mm[zlast]]+VELTime[zlast][0]

    if digraph[Index][3]!=0:
        zlast=digraph[Index][3]
        if Vflag[zlast][0]==0:
            VELTime[zlast][0],VELTime,Vflag=CalculateLastOperation(s1,s2,mm,digraph,zlast,VELTime,Vflag,time,f_index)
            Vflag[zlast][0]=1
        t2=time[f_index][s1[zlast]][s2[zlast]-1][mm[zlast]]+VELTime[zlast][0]
    if t1>t2:
        VELTime[Index][0]=t1
    else:
        VELTime[Index][0] = t2
    Vflag[Index][0]=1
    value=VELTime[Index][0]
    return value,VELTime,Vflag

def CalculateNextOperation(s1,s2,mm,digraph,next,VELTime,Vflag,LastNode,time,f_index):
    Index=next
    e=0
    t1=e;t2=e
    if digraph[Index][0]==0:
        t1=LastNode-time[f_index][s1[Index]][s2[Index]-1][mm[Index]]
    else:
        znext=digraph[Index][0]
        if Vflag[znext][1]==0:
            VELTime[znext][1],VELTime,Vflag=CalculateNextOperation(s1,s2,mm,digraph,znext,VELTime,Vflag,LastNode,time,f_index)
            Vflag[znext][1]=1
        t1=VELTime[znext][1]-time[f_index][s1[Index]][s2[Index]-1][mm[Index]]

    if digraph[Index][1]==0:
        t2=LastNode-time[f_index][s1[Index]][s2[Index]-1][mm[Index]]
    else:
        znext = digraph[Index][1]
        if Vflag[znext][1]==0:
            VELTime[znext][1],VELTime,Vflag=CalculateNextOperation(s1,s2,mm,digraph,znext,VELTime,Vflag,LastNode,time,f_index)
            Vflag[znext][1]=1
        t2=VELTime[znext][1]-time[f_index][s1[Index]][s2[Index]-1][mm[Index]]

    if t1 > t2:
        VELTime[Index][1] = t2
    else:
        VELTime[Index][1] = t1
    Vflag[Index][1] = 1
    value = VELTime[Index][1]
    return value, VELTime, Vflag


def FindCriticalPathDHFJSP(p_chrom,m_chrom,FJ,f_index,N,H,time,TM):
    SH=len(p_chrom)

    digraph=np.zeros(shape=(SH,4),dtype=int) #disjunctive graph
    dflag = np.zeros(shape=(SH, 4), dtype=int) #0-1 flag to note whther this node is added into graph
    s1 = p_chrom
    s2 = np.zeros(SH, dtype=int)
    p = np.zeros(N, dtype=int)
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
    #build disjuntive graph
    for i in range(SH):
        to=s1[i];tm=mm[i]
        #find next operation of current operation
        for j in range(i+1,SH):
            if to==s1[j]and dflag[i][0]==0:
                digraph[i][0]=j
                dflag[i][0]=1#next operation has been found
                digraph[j][2]=i
                dflag[j][2]=1#former operation has been found
            if tm==mm[j]and dflag[i][1]==0:
                digraph[i][1]=j
                dflag[i][1]=1#next operation has been found
                digraph[j][3]=i
                dflag[j][3]=1#former operation has been found
            if (dflag[i][0]==1 or s2[i]==H[s1[i]]) and dflag[i][1]==1:
                break
    VELTime=np.zeros(shape=(SH,2))
    Vflag = np.zeros(shape=(SH, 2))
    e=0

    #width first travse for graph
    level=np.zeros(SH,dtype=int)
    lenn=0;L=1
    while lenn<SH:
        for i in range(SH):
            if s2[i]==L:
                level[lenn]=i
                lenn=lenn+1
        L=L+1
    #calculate the earlyest start time
    LastNode=e
    for i in range(SH):
        Index=level[i]
        VELTime[Index][0]=e;Vflag[Index][0]=1
        t1=e;t2=e
        if digraph[Index][2]!=0:
            last=digraph[Index][2]
            t1=time[f_index][s1[last]][s2[last]-1][mm[last]]+VELTime[last][0]

        if digraph[Index][3]!=0:
            last=digraph[Index][3]
            if Vflag[last][0]==0:
                VELTime[last][0],VELTime,Vflag=CalculateLastOperation(s1,s2,mm,digraph,last,VELTime,Vflag,time,f_index)
                Vflag[last][0]=1
            t2=time[f_index][s1[last]][s2[last]-1][mm[last]]+VELTime[last][0]
        if t1>t2:
            VELTime[Index][0]=t1
            Vflag[Index][0] = 1
        else:
            VELTime[Index][0] = t2
            Vflag[Index][0] = 1
        if digraph[Index][0]==0:
            t1=time[f_index][s1[Index]][s2[Index]-1][mm[Index]]+VELTime[Index][0]
            if t1>LastNode:
                LastNode=t1

    # calculate the latest start time
    for i in range(SH-1,-1,-1):
        Index=level[i]
        if digraph[Index][0]==0:
            if digraph[Index][1]==0:
                VELTime[Index][1]=LastNode-time[f_index][s1[Index]][s2[Index]-1][mm[Index]]
                Vflag[Index][1] = 1
            else:
                t1=LastNode-time[f_index][s1[Index]][s2[Index]-1][mm[Index]]
                next=digraph[Index][1]
                if Vflag[next][1]==0:
                    VELTime[next][1],VELTime,Vflag =CalculateNextOperation(s1,s2,mm,digraph,next,VELTime,Vflag,LastNode,time,f_index);
                    Vflag[next][1] = 1
                t2=VELTime[next][1]-time[f_index][s1[Index]][s2[Index]-1][mm[Index]]
                if t1>t2:
                    VELTime[Index][1]=t2;Vflag[Index][1] = 1
                else:
                    VELTime[Index][1] = t1;Vflag[Index][1] = 1
            continue
        t1=e;t2=e
        if digraph[Index][0]!=0:
            next=digraph[Index][0]
            t1=VELTime[next][1]-time[f_index][s1[Index]][s2[Index]-1][mm[Index]]
            if digraph[Index][1]==0:
                VELTime[Index][1]=t1;Vflag[Index][1] = 1
                continue
        if digraph[Index][1]!=0:
            next=digraph[Index][1]
            if Vflag[next][1]==0:
                VELTime[next][1],VELTime,Vflag=CalculateNextOperation(s1,s2,mm,digraph,next,VELTime,Vflag,LastNode,time,f_index);
                Vflag[next][1] = 1
            t2=VELTime[next][1]-time[f_index][s1[Index]][s2[Index]-1][mm[Index]]
        if t1 > t2:
            VELTime[Index][1] = t2;Vflag[Index][1] = 1
        else:
            VELTime[Index][1] = t1;Vflag[Index][1] = 1

    #the latest start time miner the earlyest start time, the nodes with minimum idletime is on the critical path
    Idletime=np.zeros(SH)
    for i in range(SH):
        Idletime[i]=VELTime[i][1]-VELTime[i][0]
    #find minimum idletime operation
    MinIndex=0;

    MinIdleT=Idletime[0]

    for i in range(1,SH):
        if MinIdleT>Idletime[i]:
            MinIndex=i
            MinIdleT=Idletime[i]
    CritcalPath=[]
    for i in range(SH):
        if MinIdleT==Idletime[i]:
            CritcalPath.append(i)
    L=len(CritcalPath)
    block=0
    CritcalBlock=[]
    CritcalBlock.append([])
    CritcalBlock[block].append(CritcalPath[0])
    m0=mm[CritcalPath[0]]
    for i in range(1,L):
        m1=mm[CritcalPath[i]]
        if m0==m1:
            CritcalBlock[block].append(CritcalPath[i])
        else:
            block=block+1
            CritcalBlock.append([])
            CritcalBlock[block].append(CritcalPath[i])
        m0=m1
    return CritcalPath,CritcalBlock,block




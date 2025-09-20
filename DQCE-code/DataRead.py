# coding:utf-8
#f= open("../DATASET/J10M5O5.txt", "r",encoding='utf-8')
import numpy as np
import sys

def DataReadDHFJSP(Filepath):
    data = []
    enter='\n'
    # for example "../DATASET/J10M5O5.txt"
    with open(Filepath, "r",encoding='utf-8') as f1:
        for line in f1:
            temp=line.split(' ');
            l=len(temp)
            for i in range(l):
                if temp[i] != enter:
                    data.append(int(temp[i]))
    N=data[0]
    F=data[1]
    TM=data[2]
    H=np.zeros(N,dtype=int)
    NM=[]

#read number of selectable machine for each operation
    p=5 #index in date array
    for f in range(F):
        for j in range(N):
            H[j] = data[p]
            p=p+2
            temp = []
            for o in range(int(H[j])):
                temp.append(data[p])
                NM1 = data[p]
                p=p+1
                for k in range(int(NM1)):
                    p = p + 1
                    p = p + 1
                p=p+1
            p=p+1
            NM.append(temp)

    SH=int(np.sum(H))
    opmax=int(max(H))
    M=np.zeros(shape=(N,opmax,TM),dtype=int)
    time = np.zeros(shape=(F,N, opmax, TM))

    p = 5  # index in date array
    for f in range(F):
        for j in range(N):
            H[j] = data[p]
            p = p + 2
            for o in range(int(H[j])):
                NM1=data[p]
                p = p + 1
                for k in range(int(NM1)):
                    M[j][o][k] = data[p]
                    p = p + 1
                    t = int(M[j][o][k])
                    time[f][j][o][t - 1] = data[p]
                    p = p + 1
                p = p + 1
            p = p + 1

    f1.close()
    # convert benchmark to variable like number of job N, number of machine M, number of operation of each job H
    # number of seletable machine of each operation NM, total operation number SH,
    # the processing time of each operation on each selectable machine time
    ProF=np.zeros((N,F))
    for f in range(F):
        for i in range(N):
            toTime=0
            for j in range(int(H[i])):
                averT=0
                NM1=int(NM[i][j])
                for k in range(NM1):
                    mc=M[i][j][k]-1
                    averT=averT+time[f][i][j][mc]
                averT=averT/NM1
                toTime=toTime+averT
            ProF[i][f]=toTime
    for i in range(N):
        tot=0
        for f in range(F):
            tot=tot+ProF[i][f]
        for f in range(F):
            ProF[i][f]=ProF[i][f]/tot
    return N,F,TM,H,SH,NM,M,time,ProF
'''
data = f.read()
N=int(data[1])
TM=int(data[3])
H=np.zeros(1,N)
print(N)
print(TM)'''

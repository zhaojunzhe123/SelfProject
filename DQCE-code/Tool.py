# coding:utf-8
import numpy as np

def mylistRound(arr):
    l=len(arr)
    for i in range(l):
        arr[i]=round(arr[i])
    return arr

def find_all_index(arr, item):
    return [i for i, a in enumerate(arr) if a == item]

def find_all_index_not(arr, item):
    l=len(arr)
    flag=np.zeros(l)
    index=find_all_index(arr,item)
    flag[index]=1
    not_index=find_all_index(flag,0)
    return not_index

def NDS(fit1,fit2):
    v=0
    dom_less = 0;
    dom_equal = 0;
    dom_more = 0;
    for k in range(2):
        if fit1[k] > fit2[k]:
            dom_more = dom_more + 1;
        elif fit1[k] == fit2[k]:
            dom_equal = dom_equal + 1;
        else:
            dom_less = dom_less + 1;

    if dom_less == 0 and dom_equal != 2:
        v = 2
    if dom_more == 0 and dom_equal != 2:
        v = 1
    return v

def Ismemeber(item,list):
    l=len(list)
    flag=0
    for i in range(l):
        if list[i]==item:
            flag=1
            break
    return flag

def DeleteReapt(QP,QM,QF,QFit,ps):
    row=np.size(QFit,0)
    i=0
    while i<row:
        if i>=row:
            #print('break 1')
            break

        F=QFit[i,:]
        j=i+1
        while j<row:
            if QFit[j][0]==F[0] and QFit[j][1]==F[1]:
                QP = np.delete(QP, j, axis=0)
                QM = np.delete(QM, j, axis=0)
                QF = np.delete(QF, j, axis=0)
                QFit = np.delete(QFit, j, axis=0)
                j=j-1
                row=row-1
                if row<2*ps+1:
                    break
            j=j+1
        i=i+1
        if row < 2 * ps + 1:
            #print('break 2')
            break
    return QP,QM,QF,QFit

def DeleteReaptE(QP,QM,QF,QFit): #for elite strategy
    row=np.size(QFit,0)
    i=0
    while i<row:
        if i>=row:
            #print('break 1')
            break

        F=QFit[i,:]
        j=i+1
        while j<row:
            if QFit[j][0]==F[0] and QFit[j][1]==F[1]:
                QP = np.delete(QP, j, axis=0)
                QM = np.delete(QM, j, axis=0)
                QF = np.delete(QF, j, axis=0)
                QFit = np.delete(QFit, j, axis=0)
                j=j-1
                row=row-1
            j=j+1
        i=i+1

    return QP,QM,QF,QFit

def pareto(fitness):
    PF=[]
    L=np.size(fitness,axis=0)
    pn=np.zeros(L,dtype=int)
    for i in range(L):
        for j in range(L):
            dom_less = 0;
            dom_equal = 0;
            dom_more = 0;
            for k in range(2):#number of objectives
                if (fitness[i][k] > fitness[j][k]):
                    dom_more = dom_more + 1
                elif(fitness[i][k] == fitness[j][k]):
                    dom_equal = dom_equal + 1
                else:
                    dom_less = dom_less + 1

            if dom_less == 0 and dom_equal != 2: # i is dominated by j
                pn[i] = pn[i] + 1;
        if pn[i] == 0: # add i into pareto front
            PF.append(i)
    return PF
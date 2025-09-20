# coding:utf-8
import numpy as np

# individual in population
class individual:
    def __init__(self):
        self.n = 0  # the number of being dominated of individual i
        self.p = []  # the list of index which dominate individual i

# each pareto front
class front:
    def __init__(self):
        self.f = []  # the list of ith front

def FastNDS(fitness,ps):
    # store the indexs of top rank ps solutions
    TopRank=[]
    L=np.size(fitness,0)
    # create a object array
    S=[] # solutions
    for i in range(L):
        S.append(individual())
    # create pareto front
    F=[]
    rank=0
    F.append(front())
    #find first front
    for i in range(L):
        for j in range(L):
            dom_less = 0;
            dom_equal = 0;
            dom_more = 0;
            for k in range(2):
                if (fitness[i][k] > fitness[j][k]):
                    dom_more = dom_more + 1
                elif(fitness[i][k] == fitness[j][k]):
                    dom_equal = dom_equal + 1
                else:
                    dom_less = dom_less + 1

            if dom_less == 0 and dom_equal != 2: # i is dominated by j
                S[i].n = S[i].n + 1;
            elif dom_more == 0 and dom_equal!= 2: # i dominate j
                S[i].p.append(j)
        if S[i].n == 0: # add i into pareto front
            F[rank].f.append(i)
    # find the other front
    while len(F[rank].f)>0:
        Q=[]
        fL=len(F[rank].f)
        for i in range(fL):
            x=F[rank].f[i]
            if len(S[x].p)>0:
                pL=len(S[x].p)
                for j in range(pL):
                    k=S[x].p[j]
                    S[k].n=S[k].n-1
                    if S[k].n==0:
                        Q.append(k)
        rank=rank+1
        F.append(front())#create a new front
        F[rank].f=Q

    # crowding distance strategy
    fL=len(F);fL=fL-1#the last front is empty
    obj=fitness
    currentindex=0
    for Frt in range(fL):
        ffL=len(F[Frt].f)
        y=np.zeros(shape=(ffL,2))
        for i in range(ffL):
            t=F[Frt].f[i]
            y[i,:]=obj[t,:]
        crowd=[]
        crowd=np.zeros(shape=(ffL,2))
        for i in range(2):
            y0=y[:,i]
            sort_index=np.argsort(y0)# ascending sort
            fmin=y[sort_index[0]][i]
            fmax=y[sort_index[ffL-1]][i]
            if ffL==1 or fmax==fmin:
                for j in range(ffL):
                    crowd[j][i]=1
            else:
                sort_obj=[]
                for j in range(ffL):
                    sort_obj.append(y[sort_index[j]][i])
                gap=fmax-fmin
                for j in range(ffL):
                    if j==0:
                        crowd[sort_index[j]][i]=(sort_obj[1]-sort_obj[0])/gap
                    elif j==ffL-1:
                        crowd[sort_index[j]][i] = (sort_obj[j] - sort_obj[j-1]) / gap
                    else:
                        crowd[sort_index[j]][i] = (sort_obj[j+1] - sort_obj[j - 1]) / gap
        crowd[:,0]=crowd[:,0]+crowd[:,1]
        tmp=[];sort_index=[]
        sort_index=np.argsort(-crowd[:,0])# add '-' means sort in descending
        #tmp=F[Frt].f #light copy
        l = len(sort_index)
        for i in range(l):#deep copy
            tmp.append(F[Frt].f[i])
        for i in range(l):
            F[Frt].f[i]=tmp[sort_index[i]]

    count=0
    for Frt in range(fL):
        ffL=len(F[Frt].f)
        for i in range(ffL):
            TopRank.append(F[Frt].f[i])
            count=count+1
            if count==ps:
                break
        if count==ps:
            break
    return TopRank
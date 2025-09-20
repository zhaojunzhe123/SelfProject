# coding:utf-8
import copy

import numpy as np
import random
import math

def weight_init(ps,M,T):
    M=2
    Popsize=ps-1
    weights = np.zeros((Popsize, M));
    count = 0;
    for i in range(Popsize):# % 用循环遍历组合保证了权重向量互不相等
        weights[count][0] = i / ps;
        weights[count][1] = 1 - i / ps;
        count = count + 1;
    distance = np.zeros((Popsize, Popsize));
    neighbour = np.zeros((Popsize, T),dtype=int);
    for i in range(Popsize):
        for j in range(i+1,Popsize):
            A = copy.copy(weights[i,:]);
            B = copy.copy(weights[j,:]);
            x=np.mat((A - B))
            y=x.T
            z=x*y
            z=z[0]

            distance[i, j] = copy.copy(z)
            distance[j, i] = copy.copy(distance[i, j]);

        sindex= np.argsort(distance[i,:]);
        neighbour[i,:]=copy.copy(sindex[0: T]);
    return weights,neighbour

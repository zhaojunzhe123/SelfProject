# coding:utf-8
import copy

import numpy as np
import random
import math

def Update(weight,pre_p,pre_m,pre_f,new_p,new_m,new_f,zzold,zz,refpoint):
    flag=0;
    part = abs(zz - refpoint);
    newobj = max(weight* part);
    part = abs(zzold - refpoint);
    oldobj = max(weight* part);
    if newobj < oldobj:
        pre_p = new_p;
        pre_m = new_m;
        pre_f = new_f;
        flag = 1;
    return pre_p,pre_m,pre_f,flag

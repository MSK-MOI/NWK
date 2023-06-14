# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 16:04:50 2021

@author: jienzhu
"""

import numpy as np
import math
from scipy.optimize import linprog
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import networkx as nx

class UnbalancedOMT:
    def UnbalancedDistance(self):
        A=nx.shortest_path_length(self.G)
        A=dict(A)
        D=[A[i][j] for i in range(self.n) for j in range(self.n)]
        D=np.reshape(D,[self.n,self.n])
        DD=self.gamma*np.ones([self.n+1,self.n+1])
        DD[0:self.n,0:self.n]=D
        DD[self.n,self.n]=0
        self.DD=DD
        DD=DD.flatten()
            
        b=np.concatenate((self.rho1,self.rho0))
        res = linprog(DD, A_eq=self.AA, b_eq=b, bounds=[(0,None)]*(self.n+1)*(self.n+1))
        d=res.fun
        C=res.x
        C=np.reshape(C,[self.n+1,self.n+1])
        self.C=C
        self.d=d
        return d
    
    def interpolate(self,t):
        try:
            self.C
        except:
            self.UnbalancedDistance()

        threshold=1e-2
        rhot=np.zeros(self.n)
        if t==0:
            rhot=self.rho0
        elif t==1:
            rhot=self.rho1
        else:
            for i in range(self.n):
                for j in range(self.n):
                    if self.C[i,j]>threshold:
                        p=nx.shortest_path(self.G, source=i, target=j)
                        l=len(p)-1
                        loc=t*l
                        low=math.floor(loc)
                        up=math.ceil(loc)
                        ratio_low=loc-low
                        rhot[p[low]]=rhot[p[low]]+self.C[i,j]*(1-ratio_low)
                        rhot[p[up]]=rhot[p[up]]+self.C[i,j]*ratio_low
            for i in range(self.n):
                if self.C[i,self.n]>threshold:
                    rhot[i]=rhot[i]+self.C[i,self.n]*(1-t)
                if self.C[self.n,i]>threshold:
                    rhot[i]=rhot[i]+self.C[self.n,i]*t
        return rhot
    
    def draw(self,rhot):
        try:
            self.pp
        except:
            self.pp=nx.spring_layout(self.G)
        nx.draw(self.G,pos=self.pp,node_color=rhot,cmap=plt.cm.Blues,with_labels=True)
        
    def update(self,t):    
        self.draw(self.interpolate(t))
        
    def interpolate_path(self,N=10,output_name="tmp.gif"):
        fig=plt.figure()
        ani = animation.FuncAnimation(fig, self.update, frames=np.divide(list(range(N+1)),N))
        ani.save(output_name,writer='pillow')
    
    def __init__(self,rho0,rho1,G:nx.classes.graph.Graph,gamma):
        self.rho0=rho0
        self.rho1=rho1
        self.n=len(rho1)
        self.G=G
        self.gamma=gamma
        
        List1=list(range(0,(self.n+1)*(self.n+1)-1,self.n+1))
        List2=list(range(self.n+1))
        AA=None
        for k in range(self.n):
            row_temp=np.zeros((self.n+1)*(self.n+1))
            row_temp[np.add(List1,k)]=1
            if AA is None:
                AA=row_temp
            else:
                AA=np.vstack((AA,row_temp))
    
        for k in range(self.n):
            row_temp=np.zeros((self.n+1)*(self.n+1))
            row_temp[np.add(List2,k*(self.n+1))]=1
            AA=np.vstack((AA,row_temp))
        self.AA=AA
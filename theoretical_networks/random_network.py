# -*- coding: utf-8 -*-
"""
Random Network

PHY3030 - Projet de fin d'études
Network Science

Patrice Béchard
janvier 2017
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import networkx as nx

def main():
    
    N=100                        #number of nodes
    p=0.05                       #prob of two nodes being connected
    L=0                         #number of links set to zero at beginning
    
    position=np.zeros([N,2])    #position of the nodes
    adjacency=np.zeros([N,N])   #adjacency matrix
    kNodes=np.zeros(N)          #degree for each node
    kDistribution=np.zeros(N)   #number of nodes having
    
    
    for i in range(N):
        position[i,0]=random.random()
        position[i,1]=random.random()
    
    for i in range(N):
        for j in range(i+1,N):
            if random.random()<p:
                adjacency[i,j]=1
                adjacency[j,i]=1
                L+=1
                
                
    for i in range(N):
        for j in range(N):
            if adjacency[i,j]==1:
                plt.plot([position[i,0],position[j,0]],[position[i,1],position[j,1]],'b',linewidth=0.3)
                kNodes[i]+=1
                kNodes[j]+=1
    for i in range(N):
        kDistribution[kNodes[i]]+=1
        
    plt.plot(position[:,0],position[:,1],'r.',markersize=10)
    plt.axis('equal')
    plt.show()
    
    plt.plot(np.linspace(0,N-1,N),kDistribution,'r+')
    plt.show()
    
    
    kmoy=2*L/N                  #degré moyen du réseau
    
    print('Nombre de noeuds N : %d \nProbabilité de connection p : %2.2f \nNombre de liens L : %d \nDegré moyen <k> : %3.2f'%(N,p,L,kmoy))
    print('\nNombre de liens espéré : %5.0f \nDéviation Standard : %3.0f\nDegré moyen espéré : %3.2f\n'%(p*N*(N-1)/2,(p*(1-p)*N)**0.5,p*(N-1)))
    
    
        
main()
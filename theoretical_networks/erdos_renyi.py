# -*- coding: utf-8 -*-
# @Author: Patrice Béchard 20019173
# @Date:   2017-02-22 00:51:31
# @Last Modified time: 2017-03-30 13:29:41
#
# PHY3030 - Projet de fin d'etudes
# Erdos-Renyi Network

import matplotlib.pyplot as plt
import networkx as nx
from itertools import groupby
import numpy as np

plt.style.use('patrice')

N = 250
prob = 0.05

G=nx.erdos_renyi_graph(N,prob)
degree_sequence=sorted(nx.degree(G).values(),reverse=True) # degree sequence
d = {x:degree_sequence.count(x) for x in set(degree_sequence)}

degree = np.array(list(d.keys()))

prob_degree = np.array(list(d.values()))

meandegree = sum([degree[i]*prob_degree[i] for i in range(len(degree))])/N

prob_degree = prob_degree / N
diameter = nx.diameter(G)

fact = np.zeros(len(degree))
for i in range(len(degree)):
	fact[i] = np.math.factorial(degree[i])


degree_poisson = np.exp(-meandegree) * meandegree**degree / fact

plt.loglog(degree,prob_degree,'o')
plt.ylabel(r"distribution de degré $p_k$")
plt.xlabel(r"degré $k+1$")
plt.axvline(x=meandegree, ls='--', c='k')
plt.plot(degree,degree_poisson,'k:')
plt.text(min(degree),max(prob_degree)*0.99, r'$N = $%d , $p = $%.2f'%(N,prob) )
plt.text(min(degree),max(prob_degree)*0.85, r'$d_{max} = $%d '%diameter, fontsize=10)
plt.text(min(degree),max(prob_degree)*0.70, r'$\langle k \rangle = $%.2f '%meandegree, fontsize = 10)

plt.savefig("erdos_renyi_%d_%d_degree_dist.png"%(N,prob*100))

nx.write_gexf(G,'erdos_renyi_%d_%d.gexf'%(N,prob*100))

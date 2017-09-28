# @Author: Patrice Bechard <Patrice>
# @Date:   2017-02-22T00:24:40-05:00
# @Email:  bechardpatrice@gmail.com
# @Last modified by:   Patrice
# @Last modified time: 2017-04-14T15:02:37-04:00

# PHY3030 - Projet de fin d'etudes
# Barabasi-Albert Network

import matplotlib.pyplot as plt
import networkx as nx
from itertools import groupby
import numpy as np

plt.style.use('patrice')

N = 250
m = 6

G=nx.barabasi_albert_graph(N,m)
degree_sequence=sorted(nx.degree(G).values(),reverse=True) # degree sequence
d = {x:degree_sequence.count(x) for x in set(degree_sequence)}

degree = np.array(list(d.keys()))
prob_degree = np.array(list(d.values()))

meandegree = sum([degree[i]*prob_degree[i] for i in range(len(degree))])/N

prob_degree = prob_degree / N
diameter = nx.diameter(G)

degree_power = np.zeros(len(degree))
degree_power = degree**(-2.5) * 26

plt.loglog(degree+1,prob_degree,'o')
plt.ylabel(r"distribution de degré $p_k$")
plt.xlabel(r"degré $k+1$")
plt.axvline(x=meandegree, ls='--', c='k')
plt.plot(degree,degree_power,'k:')
plt.text(max(degree)*0.7,max(prob_degree)*0.99, r'$N = $%d'%(N) )
plt.text(max(degree)*0.7,max(prob_degree)*0.75, r'$d_{max} = $%d '%diameter, fontsize=10)
plt.text(max(degree)*0.7,max(prob_degree)*0.55, r'$\langle k \rangle = $%.2f '%meandegree, fontsize = 10)
plt.axis([min(degree)-1,max(degree)+5,min(prob_degree)*0.75,max(prob_degree)*1.25])

plt.savefig("barabasi_albert_%d_degree_dist.png"%N)

nx.write_gexf(G,'barabasi_albert_%d.gexf'%N)

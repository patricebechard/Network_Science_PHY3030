# @Author: Patrice Bechard <Patrice>
# @Date:   2016-08-13T14:15:52-04:00
# @Email:  bechardpatrice@gmail.com
# @Last modified by:   Patrice
# @Last modified time: 2017-04-14T20:49:13-04:00

'''
Visualizes and analyzes social networks
Takes input data from Lost Circles for Facebook networks
'''
import json
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.optimize import curve_fit
from scipy.misc import factorial
plt.style.use('patrice')

"""
Types :
    0 : with main hub
    1 : without main hub
"""
type=1

def new_network():
    g = nx.Graph()
    return g

def build_facebook_net():
    fbg = new_network()

    file_path = '/Users/user/Desktop/Ecole/Hiver_2017_PHY/PHY3030-Projet_de_fin_d_etudes/04_Codes/facebook_friends/patrice_bechard.json'

    #Read data as dict and drop unnecessary elements
    with open(file_path) as json_data:
        data = json.load(json_data)
    data.pop('userName')
    data.pop('userId')

    fromList = []
    toList = []

    #Add connections between your friends
    for i in data['links']:
        fromList.append(i['source'])
        toList.append(i['target'])
    if type == 0:
        #Add connections between you and your friends
        for i in range(0, len(data['nodes'])):
            fromList.append(0)  #Let 0 be your node Id
            toList.append(i+1)

        #Populate nodes for everyone in your network plus yourself
        for i in range(0, len(data['nodes'])+1):
            fbg.add_node(i, label=i)
    else:
        #Populate nodes for everyone in your network plus yourself
        for i in range(0, len(data['nodes'])):
            fbg.add_node(i, label=i)

    #Populate edges for all connections in your network
    for i in range(0, len(fromList)):
        fbg.add_edge(fromList[i], toList[i], key=i)
    json_data.close()

    return fbg

def degree_distribution(facebook_net):

    degree_sequence=sorted(nx.degree(facebook_net).values(),reverse=True) # degree sequence
    d = {x:degree_sequence.count(x) for x in set(degree_sequence)}

    degree = np.array(list(d.keys()))							#all different degrees
    prob_degree = np.array(list(d.values()))					#and probability of choosing
    numNodes = nx.number_of_nodes(facebook_net)					#size of graph

    meandegree = sum([degree[i]*prob_degree[i] for i in range(len(degree))])/numNodes #mean degree
    clustering = nx.average_clustering(facebook_net)			#mean clustering coefficient
    avg_distance = nx.average_shortest_path_length(facebook_net)

    prob_degree = prob_degree / numNodes 						#normalized probability
    diameter = nx.diameter(facebook_net)						#diameter of graph

    #degree_power = degree**-2.5 * 26

    #domain = np.logspace(1,np.log10(750),30)		#split domain to have equally spaced bins in log log
    domain = np.logspace(np.log10(4),np.log10(750),40)		#split domain to have equally spaced bins in log log
    values = np.histogram(degree,bins=domain)		#number of points in each bin
    width = np.array([domain[i+1]-domain[i] for i in range(len(domain)-1)])	#width of each bin
    values = values[0]
    newdomain = []
    for i in range(len(domain)-1):
    	newdomain.append(0.5 * (domain[i] + domain[i+1]))
    domain = newdomain

    values = (values/width)
    values /= sum(values)



    #plt.loglog(domain,values,'o')
    domain = np.linspace(min(degree),max(degree),int((max(degree)-min(degree))/10))

    values = np.histogram(degree,bins=domain)		#number of points in each bin
    values = np.array(values[0])
    newdomain = []
    for i in range(len(domain)-1):
    	newdomain.append(0.5 * (domain[i] + domain[i+1]))
    domain = np.array(newdomain)

    values = values/np.sum(values)
    
    #prob_degree = [np.sum(prob_degree[i:i+5]) for i in range(len(degree))]

    degree_power = degree**-2.2 * 2500

    plt.loglog(domain,values,'.')

    plt.ylabel(r"distribution de degré $p_k$")
    plt.xlabel(r"degré $k+1$")
    plt.axvline(x=meandegree, ls=':', c='k', lw = 1)
    plt.plot(degree,degree_power,c='k',lw=1)
    plt.text(min(domain),1e-2, r'$N = $%d'%(numNodes),fontsize=10 )
    plt.text(min(domain),8e-3, r'$d_{max} = $%d '%diameter, fontsize=10)
    plt.text(min(domain),6.5e-3, r'$\langle d \rangle = $%.3f '%avg_distance, fontsize = 10)
    plt.text(min(domain),5.2e-3, r'$\langle k \rangle = $%.2f '%meandegree, fontsize = 10)
    plt.text(min(domain),4.2e-3, r'$\langle C \rangle = $%.3f '%clustering, fontsize = 10)
    plt.axis([4,1000,3e-3,1e-1])

    if type==0:
        plt.savefig('with/with_degree.png')
    elif type==1:
        plt.savefig('without/without_degree.png')
    plt.show()

#--------------------------Main------------------------------
for type in [type]:
    facebook_net = build_facebook_net()

    for i in range(facebook_net.number_of_nodes()):
    	if facebook_net.degree(i) == 0:
    		facebook_net.remove_node(i)

    degree_distribution(facebook_net)

    #distance(facebook_net)
    #Export data for use in Gephi
    if type==0:
        nx.write_gexf(facebook_net, "with/facebook_network_with.gexf")
    elif type==1:
        nx.write_gexf(facebook_net, "without/facebook_network_without.gexf")


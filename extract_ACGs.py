import networkx as nx
import operator
import os

def get_network(targetRelation):
    graph = nx.DiGraph()
    networkFile = open('./input_data/regulatory_signaling_network.txt','r')
    for line in networkFile:
        left, right, assoType, relation = line.strip().replace('"','').split('\t')
        if assoType not in targetRelation: continue
        graph.add_edge(left,right)
    networkFile.close()
    return graph

def get_DEGs():
    dicDEGs = {}
    degFile = open('./input_data/DEGs.txt','r')
    for line in degFile:
        gene, expression = line.strip().split('\t')
        dicDEGs[gene]=''
    degFile.close()
    # print ('the number of DEGs: ',len(list(dicDEGs.keys())))
    return dicDEGs

def extract_ACGs():
    print ('\n......extracting ACGs')
    sigGraph = get_network(['signaling'])
    dicDEGs = get_DEGs()

    dicNeighbors = {}
    for deg in list(dicDEGs.keys()):
        if deg not in sigGraph.nodes(): continue
        for neighbor in sigGraph.neighbors(deg):
            if neighbor in dicDEGs.keys(): continue
            if neighbor not in dicNeighbors.keys(): dicNeighbors[neighbor] = 0
            dicNeighbors[neighbor] += 1
    print ('the number of ACGs: ',len(list(dicNeighbors.keys())))

    if os.path.exists('./processed_data') == False: os.makedirs('./processed_data')
    acgFile = open('./processed_data/ACGs.txt','w+')
    acgFile.write('\n'.join(['\t'.join([acg,str(count)]) for acg, count in sorted(dicNeighbors.items(), key=operator.itemgetter(0))]))
    acgFile.close()

def get_ACGs():
    dicACGs = {}
    acgFile = open('./processed_data/ACGs.txt','r')
    for line in acgFile:
        acg, count = line.strip().replace('"','').split('\t')
        dicACGs[acg] = ''
    acgFile.close()
    # print ('the number of ACGs: ',len(list(dicACGs.keys())))
    return dicACGs
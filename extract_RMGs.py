from extract_ACGs import *
import gseapy as gp
import os

def extract_RMGs():
    print ('\n......extracting RMGs')
    moleGraph = get_network(['signaling','regulatory'])
    listDEG  = list(get_DEGs().keys())
    listACG = list(get_ACGs().keys())
    listPRG = listDEG+listACG
    moleGraph.remove_nodes_from(list(set(moleGraph.nodes())-set(listPRG)))

    ### get the largest connected component
    for nodes in sorted(list(nx.connected_components(moleGraph.to_undirected())))[1:]:
        moleGraph.remove_nodes_from(nodes)
    print ('the number of RMGs: ',len(list(moleGraph.nodes())))

    ### write RMGs
    rmgFile = open('./processed_data/RMGs.txt','w+')
    rmgFile.write('\n'.join(list(moleGraph.nodes())))
    rmgFile.close()

    perform_GO_enrichment_analysis(list(moleGraph.nodes()),'RMGs',0.05)

def get_RMGs():
    dicRMGs = {}
    fRMG = open('./processed_data/RMGs.txt')
    for line in fRMG: dicRMGs[line.strip()] = ''
    fRMG.close()
    # print ('the number of RMGs: ',len(list(dicRMGs.keys())))
    return dicRMGs

def perform_GO_enrichment_analysis(inputGenes,geneDescription,threshold):
    targetGeneSet = "GO_Biological_Process_2015" #"KEGG_2016","Reactome_2013","GO_Molecular_Function_2015","WikiPathways_2013"
    if os.path.exists('./gene_set_enrichment_analysis') == False: os.makedirs('./gene_set_enrichment_analysis')
    enr = gp.enrichr(gene_list=inputGenes, description=geneDescription, gene_sets=targetGeneSet, outdir = './gene_set_enrichment_analysis', cutoff = threshold)
    enr.res2d.head()
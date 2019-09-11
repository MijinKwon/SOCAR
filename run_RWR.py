from extract_ACGs import *
from extract_RMGs import *
import os
import glob
def get_approved_drugs():
    dicApprovedDrugs = {}
    fApprovedDrug = open('./input_data/approved_drugs.txt')
    for line in fApprovedDrug:
        dicApprovedDrugs[line.strip().replace('"','').replace("/",'').lower()]=''
    fApprovedDrug.close()
    return dicApprovedDrugs

def get_drug_target():
    dicApproved = get_approved_drugs()
    moleGraph = get_network(['signaling', 'regulatory'])
    dicDrugTarget = {}
    fDrugTaret = open('./input_data/drug_target_associations.txt','r',encoding='UTF8')
    for line in fDrugTaret:
        drug, targets = line.strip().replace('"','').split('\t')
        if drug.lower() not in dicApproved.keys(): continue
        dicDrugTarget[drug.replace("/",'').lower()] = [target for target in targets.split('|') if target in moleGraph.nodes()]
    fDrugTaret.close()

    fGoldStandard = open('./input_data/gold_standard_drugs.txt','r',encoding='UTF8')
    fGoldStandard.readline()
    for line in fGoldStandard:
        drug, targets = line.strip().replace('"','').split('\t')
    return dicDrugTarget

def process_rwr_input_file():
    print ('\n......generating RWR input files')
    geneFile = './processed_data/rwr_gene_dictionary.txt'
    pathwayFile = './processed_data/rwr_relation_dictionary.txt'
    targetFile = './processed_data/rwr_drug_target_dictionary.txt'

    moleGraph = get_network(['signaling','regulatory'])
    dicGeneID = {}
    int = 1
    for gene in moleGraph.nodes():
        dicGeneID[gene] = int
        int += 1

    foGeneFile = open(geneFile,'w+',encoding='UTF8')
    foGeneFile.write('\n'.join('\t'.join([gene, str(id)]) for gene, id in dicGeneID.items()))
    foGeneFile.close()

    foPathwayFile = open(pathwayFile,'w+',encoding='UTF8')
    for gene1, gene2 in moleGraph.edges():
        foPathwayFile.write('\t'.join([str(dicGeneID[gene1]),str(dicGeneID[gene2]),'->', str(1.0/len(moleGraph[gene1]))]) +'\n')
    foPathwayFile.close()

    dicDrugTarget = get_drug_target()
    foTargetFile = open(targetFile,'w+',encoding='UTF8')
    for drug, targets in dicDrugTarget.items():
        for target in targets:
            if target not in dicGeneID.keys(): continue
            foTargetFile.write('\t'.join([drug.title(), str(dicGeneID[target])])+'\n')
    foTargetFile.close()

    if os.path.exists('./random_walk_result') == False: os.mkdir('./random_walk_result')

def calculate_PSS(geneOption):
    print ('\n......calculating PSS')
    dicRMGs = {}
    if geneOption == 'DEGACG': dicRMGs = get_RMGs()
    elif geneOption == 'DEG': dicRMGs = get_DEGs()
    dicDrugTarget = get_drug_target()
    listDrugRWRFiles = glob.glob("./random_walk_result/*.txt")
    dicDrugScore = {}
    nTotal = len(listDrugRWRFiles)
    i = 1
    for drugFileName in listDrugRWRFiles:
        if i % 10 == 0: print('...... %s / %s' % (i, nTotal)); i += 1
        fDrugFile = open(drugFileName)
        drugName = os.path.basename(drugFileName).replace('.txt','')
        sum = 0.0
        for line in fDrugFile:
            gene, value = line.strip().split('\t')
            if not gene in dicRMGs.keys(): continue
            sum += float(value)
        numTarget = len(dicDrugTarget[drugName.lower()])
        scoreNormed = float(sum) / float(numTarget)
        fDrugFile.close()
        dicDrugScore[drugName] = scoreNormed

    listGoldStandard = get_gold_standard_drug()
    if os.path.exists('./model_evaluation') == False: os.makedirs('./model_evaluation')
    fScore = open('./model_evaluation/drug_score_%s.txt'%geneOption, 'w+', encoding='UTF8')
    fScore.write('\t'.join(['drugName', 'score','lable']) + '\n')
    for drug, score in sorted(dicDrugScore.items(), key=operator.itemgetter(1), reverse=True):
        if drug.lower() in listGoldStandard: lable = '1'
        else: lable = '0'
        fScore.write('\t'.join([drug, str(score),lable]) + '\n')
    fScore.close()

def get_gold_standard_drug():
    listGoldStandard = []
    fGoldStandard = open('./input_data/gold_standard_drugs.txt','r',encoding='UTF8')
    fGoldStandard.readline()
    listGoldStandard = [line.strip().split('\t')[0].lower() for line in fGoldStandard]
    fGoldStandard.close()
    return listGoldStandard

def get_drug_score(geneOption):
    dicDrugScore,dicDrugLable = {}, {}
    fScore = open('./model_evaluation/drug_score_%s.txt'%geneOption, encoding='UTF8')
    fScore.readline()
    for line in fScore:
        drug, score, lable = line.strip().split('\t')
        dicDrugScore[drug] = float(score)
        dicDrugLable[drug] = int(lable)
    fScore.close()
    return dicDrugScore, dicDrugLable
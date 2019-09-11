from run_RWR import *
from sklearn.metrics import roc_auc_score
from extract_RMGs import *

def calculate_AUROC(geneOption):
    print ('\n......calculating AUROC of model using %s'%geneOption)
    dicDrugScore, dicDrugLable = get_drug_score(geneOption)
    y_true = [int(lable) for drug, lable in sorted(dicDrugLable.items(), key=operator.itemgetter(0))]
    y_scores = [float(score) for drug, score in sorted(dicDrugScore.items(), key=operator.itemgetter(0))]

    auroc_macro = roc_auc_score(y_true,y_scores,'macro')
    print ('the number of test drugs: ',len(y_true))
    print ('the number of gold standard drugs: ',sum(y_true))
    print ('AUROC: ',round(auroc_macro,2))
    
    fAUROC = open('./model_evaluation/%s_AUROC.txt'%geneOption,'w+',encoding='UTF8')
    fAUROC.write(str(auroc_macro))
    fAUROC.close()

def predict_sensitization_mechanisms():
    print ('\n......predicting_sensitization_mechanisms')
    geneOption = 'DEGACG'
    numTopDrug, numTopRMG = 5, 50
    dicDrugScore = get_drug_score(geneOption)[0]
    dicRMGs = get_RMGs()
    listTopDrugs = sorted(dicDrugScore.items(), key=operator.itemgetter(1), reverse=True)[:numTopDrug]
    # print ('Top %s drugs: '%numTopDrug, listTopDrugs)

    for drug,score in listTopDrugs:
        fDrugFile = open('./random_walk_result/%s.txt'%drug,'r',encoding='UTF8')
        dicRMGScore = {}
        for line in fDrugFile:
            gene, score = line.strip().split('\t')
            if gene not in dicRMGs.keys(): continue
            dicRMGScore[gene] = float(score)
        fDrugFile.close()
        listTopRMGs = [gene for gene, score in sorted(dicRMGScore.items(), key=operator.itemgetter(1),reverse=True)][:numTopRMG]
        perform_GO_enrichment_analysis(listTopRMGs, '%s_top_%s'%(drug,numTopRMG), 0.05)


# -*- coding: utf-8 -*-
import fileinput
import os

snpimportancefile_ntree_1000 = "P:\\modeling\\new\\RFdata\\ntree1000_mtry315_all_chi_filtering_p005_data_class_rf_importance.dat"
snpimportancefile_ntree_800 = "P:\\modeling\\new\\RFdata\\ntree800_mtry319_all_chi_filtering_p005_data_class_rf_importance.dat"
snpimportancefile_ntree_600 = "P:\\modeling\\new\\RFdata\\ntree600_mtry316_all_chi_filtering_p005_data_class_rf_importance.dat"
gene_rs_file = "P:\\modeling\\new\\sum_gene_rs_id.dat"
#写入gene importance的文件
geneimportancefile_ntree_1000 = "P:\\modeling\\new\\RFdata\\ntree1000_mtry315_all_chi_filtering_p005_data_class_rf_importance_of_gene.dat"
geneimportancefile_ntree_800 = "P:\\modeling\\new\\RFdata\\ntree800_mtry319_all_chi_filtering_p005_data_class_rf_importance_of_gene.dat"
geneimportancefile_ntree_600 = "P:\\modeling\\new\\RFdata\\ntree600_mtry316_all_chi_filtering_p005_data_class_rf_importance_of_gene.dat"

def getimportanceDicfromfile(f):
    d = {}
    for line in fileinput.input(f):
        linepartDic={}
        if line.startswith("\"rs"):
            arr = line.strip().split(" ")
            id = arr[0].strip("\"")
            linepartDic["disease"] = arr[1]
            linepartDic["healthy"] = arr[2]
            linepartDic["MeanDecreaseAccuracy"] = arr[3]
            linepartDic["MeanDecreaseGini"] =arr[4]
            d[id] = linepartDic
        else:
            continue
    return d
def getrsIDdicforGene(f):
    d= {}
    for line in fileinput.input(f):
        rslst = []
        arr = line.strip().split("\t")
        geneid = arr[0]
        rslst = arr[1].split(",")
        d[geneid] = rslst

    return d

def getgeneimportancefile(snpimdic, gene2rsdic, geneimportancef):
    geneimportanceDic = {}
    nullrsgene = []
    for gene in gene2rsdic:
        rsIDlst = gene2rsdic[gene]
        candidatersforgenelst = []
        candidatersimportanceforgenelst = []
        for rs in rsIDlst:
            if snpimdic.has_key(rs):
                rsscoredic = snpimdic[rs]
                giniimportance = rsscoredic["MeanDecreaseGini"]
                candidatersforgenelst.append(rs)
                candidatersimportanceforgenelst.append(giniimportance)
            else:
                continue
        # sumimportance= sum(candidatersimportanceforgenelst)
        if len(candidatersforgenelst) >0:
            geneimportanceDic[gene] =[candidatersforgenelst,candidatersimportanceforgenelst,0]

        else:

            nullrsgene.append(gene)
    print  geneimportancef +": " + str(len(nullrsgene)) + " gene has no significant snp"


    fw= open(geneimportancef,"w")
    head = "geneID\trsIDinchiSet\tindividuleImportanceforEachsnp\tsumimportanceforgene\n"
    fw.write(head)
    for k in geneimportanceDic.keys():
        geneinfolst = geneimportanceDic[k]
        rsidset = geneinfolst[0]
        rsimportancelst = geneinfolst[1]
        sumimportance4gene = geneinfolst[2]
        newline =k + "\t"+",".join(rsidset)+"\t"+ ",".join(rsimportancelst)+"\t"+str(sumimportance4gene)+"\n"
        fw.write(newline)

    fw.close


snpimportanceDic_ntree_1000 = getimportanceDicfromfile(snpimportancefile_ntree_1000)
snpimportanceDic_ntree_800 = getimportanceDicfromfile(snpimportancefile_ntree_800)
snpimportanceDic_ntree_600 = getimportanceDicfromfile(snpimportancefile_ntree_600)
rsforgeneDic = getrsIDdicforGene(gene_rs_file)
getgeneimportancefile(snpimportanceDic_ntree_600,rsforgeneDic,geneimportancefile_ntree_600)
getgeneimportancefile(snpimportanceDic_ntree_800,rsforgeneDic,geneimportancefile_ntree_800)
getgeneimportancefile(snpimportanceDic_ntree_1000,rsforgeneDic,geneimportancefile_ntree_1000)

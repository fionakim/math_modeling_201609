# -*- coding: utf-8 -*-
import pandas
import fileinput
import math
import numpy as np
from scipy.stats import chisquare
import random

alldatafile = "P:\\modeling\\new\\multi_phenos\\all_phenotype_snp_file.dat"
phenotype_matrix_file_prefix = "P:\\modeling\\new\\multi_phenos\\all_phenotype_snp_file_"
file2writelst = []
mainpf = "P:\\modeling\\new\mainphenoype_logistic\\new_finalnumericgenotype.dat"
rsgenetypelogfile = "P:\\modeling\\rsnumericlog.dat"


# ======================================函数定义区
# 抽样函数区
def getrateDicfromfile(f):
    d = {}
    id = ''
    i = 0
    for line in fileinput.input(f):
        sd = {}
        i = i + 1
        if i > 2:
            arr = line.strip().split("\t")
            id = arr[0]
            onenumber = 0.0
            zeronumber = 0.0
            halfnumber = 0.0
            s = len(arr) - 1
            for e in arr:
                if e == '0':
                    zeronumber = zeronumber + 1
                if e == '1':
                    onenumber = onenumber + 1
                if e == '0.5':
                    halfnumber = halfnumber + 1
            if not s == onenumber + zeronumber + halfnumber:
                print 'not match'
            sd['0'] = zeronumber
            sd['1'] = onenumber
            sd['0.5'] = halfnumber
            sd['sum'] = s
        d[id] = sd

    return d


def getDicfromfile(f):
    d = {}
    for line in fileinput.input(f):
        subd = {}
        arr = line.strip().split("\t")
        subd['0'] = arr[1]
        subd['0.5'] = arr[2]
        subd['1'] = arr[3]
        d[arr[0]] = subd
    return d


def judgeifmatchhw(aa, ac, cc):  # aa是D6，ac是D7，cc是D8
    totalnumber = aa + ac + cc  # d9
    vdic = {}
    anumber = aa * 2 + ac  # D11
    cnumber = cc * 2 + ac  # D12
    totalloci = anumber + cnumber  # D13
    if totalloci == 0:
        return -1
    else:
        af = float(anumber) / float(totalloci)  # F11
        cf = float(cnumber) / float(totalloci)  # F12

        preaanumber = totalnumber * math.pow(af, 2)  # E6
        preacnumber = totalnumber * af * 2 * cf  # E7
        preccnumber = totalnumber * math.pow(cf, 2)  # E8

        # chisquare = math.pow((aa - preaanumber), 2) / preaanumber + math.pow((ac - preacnumber),
        #                                                                      2) / preacnumber + math.pow(
        #     (cc - preccnumber), 2) / preccnumber
        expected = [preaanumber, preacnumber, preccnumber]
        observed = [aa, ac, cc]
        ch, p = chisquare(observed, expected)
        print
        vdic['chi'] = ch
        vdic['pvalue'] = p
        return vdic


# ===============
def getDicfromfile(f):
    d = {}
    for line in fileinput.input(f):
        subd = {}
        arr = line.strip().split("\t")
        subd['0'] = arr[1]
        subd['0.5'] = arr[2]
        subd['1'] = arr[3]
        d[arr[0]] = subd
    return d


def getrateDicfromfile(f):
    d = {}
    i = 0
    errSNPset = set()
    for line in fileinput.input(f):
        sd = {}
        i = i + 1
        if i > 2:
            arr = line.strip().split("\t")
            id = arr[0]
            onenumber = 0.0
            zeronumber = 0.0
            halfnumber = 0.0
            s = len(arr) - 1
            for e in arr:
                if e == '0' or e == 0:
                    zeronumber = zeronumber + 1
                if e == '1' or e == 1:
                    onenumber = onenumber + 1
                if e == '0.5' or e == 0.5:
                    halfnumber = halfnumber + 1
            if not s == onenumber + zeronumber + halfnumber:
                print 'not match'
            if zeronumber == 0 or onenumber == 0 or halfnumber == 0:
                errSNPset.add(id)
            else:
                sd['0'] = zeronumber
                sd['1'] = onenumber
                sd['0.5'] = halfnumber
                sd['sum'] = s
                d[id] = sd

    return d, errSNPset


def writeDic2file(d, f):
    fw = open(f, "w")
    fw.write("rsID\tpValue\t-LOG10pValue\n")
    for k in d.keys():
        info = d[k]
        newline = k + "\t" + str(info[0]) + "\t" + str(info[1]) + "\n"
        fw.write(newline)
    fw.close()


def transposefile(srcfile, targetfile):
    f2 = open(targetfile, "w")
    newlst = []
    col = 0
    row = 0
    for line in fileinput.input(srcfile):
        row = row + 1
        patcharr = line.strip().split("\t")
        col = len(patcharr)
        for e in patcharr:
            newlst.append(e)
    for i in range(col):
        newlineitemlst = []
        for j in range(row):
            newlineitemlst.append(newlst[i + j * col])
        newline = "\t".join(newlineitemlst) + "\n"
        f2.write(newline)


def writelist2f(lst, f, h):
    with open(f, "w") as fw:
        fw.write(h + "\n")
        for e in lst:
            fw.write(e + "\n")
    fw.close()


def getsamplefile(file1, file2, dataclass):
    caselinelst = []
    controllinelst = []
    head = ""
    i = 0
    for line in fileinput.input(file1):
        i = i + 1
        if i == 1:
            head = line.strip()
        else:
            arr = line.strip().split("\t")
            if arr[1] == '0':
                caselinelst.append(line.strip())
            else:
                controllinelst.append(line.strip())
    if dataclass == "case":
        writelist2f(caselinelst, file2, head)
    if dataclass == "control":
        writelist2f(controllinelst, file2, head)

        # 下面从得到的file里随机抽样（待补充）


def reduceDimension(old, reduced, seed):
    fw = open(reduced, "w")
    i = 0
    for line in fileinput.input(old):
        i = i + 1
        if i == 1:
            fw.write(line.strip() + "\n")
        elif i == 2:
            arr = line.strip().split("\t")
            elst = []
            j=0
            for e in arr:
                j=j+1
                if j==1:
                    elst.append(e)
                elif str(e).find("1"):
                    elst.append("diseased")
                else:
                    elst.append("healthy")

            l = "\t".join(elst)
            fw.write(l + "\n")
        else:
            arr = line.strip().split("\t")
            if arr[0].strip() in seed:
                fw.write(line.strip() + "\n")
            else:
                continue



def getlstfromfile(o,l):
    lst = []
    for line in fileinput.input(o):
        if line.startswith("\"rs"):
            arr = line.strip().split("\t")
            id = arr[0].strip("\"")
            lst.append(id)
            if len(lst)>l:
                break
    fileinput.close()
    return lst


def gettimesinmultiplesbset(slst,e):
    i = 0
    lst = []
    for sub in slst:
        if e in sub:
            lst.append('1')
            i = i+1
        else:
            lst.append('0')
            continue
    return  i,lst

def writematrix2file(matrix, f):
    newlst = []
    for line in matrix:
        strline=[]
        for e in line:
            newe = str(e)
            strline.append(newe)
        newline = "\t".join(strline)
        newlst.append(newline)
    lines = "\n".join(newlst)
    fw = open(f,"w")
    fw.write(lines)
    fw.close()


# =====================================写十个性状的源数据文件

for i in range(1, 11):
    file2write = phenotype_matrix_file_prefix + str(i) + "th_phenotype.dat"
    file2writelst.append(file2write)
    # fw = open(file2write,"w")
    # for line in fileinput.input(alldatafile):
    #     arr = line.strip().split("\t")
    #     snplst = arr[11:len(arr)-1]
    #     newline = arr[0] +"\t"+arr[i]+"\t"+"\t".join(snplst)+"\n"
    #     fw.write(newline)
    # fw.close()

# =====================================================================================================#
# --------------每个文件进行分组

# file2writelst.append(mainpf)
pvalue_threthhold = 5*math.pow(10,-3)
phenotypeallfilelstDic = {}
for pf in file2writelst:
    fileDicForthisPhenotype = {}
    fallcase = pf + ".allcase.dat"
    fallcontrol = pf + ".allcontrol.dat"
    fileDicForthisPhenotype["allcase"] = fallcase
    fileDicForthisPhenotype["allcontrol"] = fallcontrol
    # getsamplefile(pf, fallcase, "case")
    # getsamplefile(pf, fallcontrol, "control")
    phenotypeallfilelstDic[pf] = fileDicForthisPhenotype

# 下面生成每个文件的转置文件

for pf in file2writelst:
    case_control_dic = phenotypeallfilelstDic[pf]
    trans_all_case_file = case_control_dic["allcase"] + ".transpose.dat"
    trans_all_control_file = case_control_dic["allcontrol"] + ".transpose.dat"
    case_control_dic["trans_allcase"] = trans_all_case_file
    case_control_dic["trans_allcontrol"] = trans_all_control_file
    # transposefile(case_control_dic["allcase"], trans_all_case_file)
    # transposefile(case_control_dic["allcontrol"], trans_all_control_file)

#下面生成每个性状的卡方检验

for pf in file2writelst:
    snpSignificancePvaluefile = pf + ".snpSignificancePvaluefile.dat"
    case_control_dic = phenotypeallfilelstDic[pf]
    case_control_dic["snpSignificancePvaluefile"] = snpSignificancePvaluefile
    # casersDic, caseErrsnp = getrateDicfromfile(case_control_dic["trans_allcase"])
    # controlDic, controlErrsnp = getrateDicfromfile(case_control_dic["trans_allcontrol"])
    # sumErrsnp = caseErrsnp | controlErrsnp
    # remainsnp = set(casersDic.keys()) - sumErrsnp
    # rs_pvalue_dic = {}
    # frecasedic = {}
    # frecontroldic = {}
    # for c in remainsnp:
    #     if not c.startswith('rs'):
    #         continue
    #     else:
    #         vcase = casersDic[c]
    #         vcontrol = controlDic[c]
    #         casearr = []
    #         controlarr = []
    #         if vcase.has_key('1') and vcase.has_key('0.5') and vcase.has_key('0'):
    #             casearr = [vcase['1'], vcase['0'], vcase['0.5']]
    #             frecasedic[c] = casearr
    #         else:
    #             print 'case: ' + c + " key error"
    #
    #         if vcontrol.has_key('1') and vcontrol.has_key('0.5') and vcontrol.has_key('0') and vcontrol:
    #             controlarr = [vcontrol['1'], vcontrol['0'], vcontrol['0.5']]
    #             frecontroldic[c] = controlarr
    #         else:
    #             print 'control: ' + c + " key error"
    #         if len(casearr) == 3 and len(controlarr) == 3:
    #             ch, p = chisquare(casearr, controlarr)
    #             minuslogPvalue = 0.0
    #             if p <= 0:
    #                 print "pVALUE IS ZERO :" + c + " " + str(p) + c + ":" + str(vcase) + ";" + str(vcontrol) + str(ch)
    #                 minuslogPvalue = 400
    #                 rs_pvalue_dic[c] = [p, minuslogPvalue]
    #             else:
    #                 minuslogPvalue = -math.log10(p)
    #                 rs_pvalue_dic[c] = [p, minuslogPvalue]
    #                 if p < math.pow(10, -8):
    #                     # print "strong significance: " +c + ":"+ str(vcase) + ";"+ str(vcontrol) +str(p) +" " +str(ch)
    #                     pass
    #
    #         else:
    #             print 'there is sth wrong in array '
    # writeDic2file(rs_pvalue_dic, snpSignificancePvaluefile)
    # candidatesnp_rf = []
    # for snp in rs_pvalue_dic:
    #     lst = rs_pvalue_dic[snp]
    #     if lst[0] <= pvalue_threthhold:
    #         candidatesnp_rf.append(snp)
    # case_control_dic["snp_for_rf"] = candidatesnp_rf
    # print len(candidatesnp_rf)
#==========批量跑随机森林step1:利用已选择的rsID对总文件进行降维，首先转置，再姜维，再转置回来


for pf in file2writelst:
    transposesumfile = pf + ".transpose.dat"
    reduced_snp_sumfile_by_chi_transpose = pf + ".reduced_snp.transpose.dat"
    reduced_snp_sumfile_by_chi = pf + ".reduced_snp.dat"
    case_control_dic["reduced_snp_sumfile_by_chi"] = reduced_snp_sumfile_by_chi
    case_control_dic = phenotypeallfilelstDic[pf]
    case_control_dic["sumfile_transpose"] = transposesumfile
    # transposefile(pf, transposesumfile)
    # reduceDimension(transposesumfile, reduced_snp_sumfile_by_chi_transpose, case_control_dic["snp_for_rf"])
    # transposefile(reduced_snp_sumfile_by_chi_transpose, reduced_snp_sumfile_by_chi)

# # ==========批量跑随机森林step2:调用R脚本选择了ntree=600,800,1000三个重复
#
#
# # ==== 每个重复选m个排名靠前的snp FREQUENCY代表

m = 50
FREQUENCY=3
sumset = set()
setlst = []
for pf in file2writelst:

    case_control_dic = phenotypeallfilelstDic[pf]
    importancefile_600 = case_control_dic["reduced_snp_sumfile_by_chi"] + ".ntree600.importance.dat"
    importancefile_800 = case_control_dic["reduced_snp_sumfile_by_chi"] + ".ntree800.importance.dat"
    importancefile_1000 = case_control_dic["reduced_snp_sumfile_by_chi"] + ".ntree1000.importance.dat"
    # USE R to order file
    # for (i in 1:30){filedir = paste(root,sep = "",filearr[i]); rdata <-read.table(filedir,header = TRUE);ordereddata<- rdata[order(rdata[,4],decreasing=T),];duishuzhi<-log10(ordereddata$MeanDecreaseGini);newdata<-cbind(ordereddata,duishuzhi);write.table(newdata,file = paste(root2,sep = "",filearr[i],"ordered_by_gene"),sep = "\t",row.names = TRUE,col.names = TRUE)}

    rs600_ordered = case_control_dic["reduced_snp_sumfile_by_chi"]  + " .ntree600.importance.datordered_by_gini.dat"
    rs800_ordered = case_control_dic["reduced_snp_sumfile_by_chi"]  + " .ntree800.importance.datordered_by_gini.dat"
    rs1000_ordered = case_control_dic["reduced_snp_sumfile_by_chi"]  + " .ntree1000.importance.datordered_by_gini.dat"
    rs600orderedlst=getlstfromfile(rs600_ordered,m)
    rs800orderedlst = getlstfromfile(rs800_ordered,m)
    rs1000orderedlst = getlstfromfile(rs1000_ordered,m)
    sharedset = set(rs600orderedlst) & set(rs800orderedlst) & set(rs1000orderedlst)
    case_control_dic["significant_snp_from_3_rf"] = sharedset
    sumset = sumset | sharedset
    setlst.append(sharedset)
    print "for " +pf +": " +str(len(sharedset)) +" snp significant"

n=0
data2writeDic = {}

for snp in sumset:
    frequency,yesornolst = gettimesinmultiplesbset(setlst,snp)
    #统计总集中RS在字迹出现的频率 结果有待商榷
    if frequency >=FREQUENCY:
        n=n+1
        data2writeDic[snp] = yesornolst
        print snp +": "+str(frequency)
print str(n) +" snp has more than " +str(FREQUENCY) +" times"


iset = sumset
for s in setlst:
    iset =iset & s
    print len(iset)
print "the last number is : " + str(len(iset)) +" the snp set is : " + str(iset)

#
# #写频次数据到文件里
# frequencyfile = "P:\\modeling\\new\\multi_phenos\\frequencyfile.dat"
# frequencyfile_trans = "P:\\modeling\\new\\multi_phenos\\frequencyfile_transpose.dat"
# fw = open(frequencyfile,"w")
# fw.write("snpID	pan	fang	jie	biye	huojiang	healthy	happy	growth	coding	love\n")
# for snp in data2writeDic.keys():
#     number = "\t".join(data2writeDic[snp])
#     newline = snp +"\t" +number+"\n"
#     fw.write(newline)
# fw.close()
# transposefile(frequencyfile,frequencyfile_trans)
#
#
# intersectionmatrix =[]
# numbermatrix = []
# for subset1 in setlst:
#     line=[]
#     noline = []
#     for subset2 in setlst:
#         inter = subset1 &subset2
#         line.append(inter)
#         noline.append(len(inter))
#     intersectionmatrix.append(line)
#     numbermatrix.append(noline)
# print numbermatrix
# subsetinterfile = "P:\\modeling\\new\\multi_phenos\\subsetinterfile.dat"
# writematrix2file(numbermatrix,subsetinterfile)

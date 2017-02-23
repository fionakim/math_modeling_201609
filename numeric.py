# -*- coding: utf-8 -*-
# __author__ = 'jinlinfang'
import fileinput

genoTypeDat = "P:\\modeling\\newgenotype.dat"
fwfile = "P:\\modeling\\numericgenotype.dat"
rsnumericlogfile = "P:\\modeling\\rsnumericlog.dat"
lineNo = 0
headline = ""
lociIDDic = {}

def getkeysforv(d, v):
    k=""
    for item in d.keys():
        if d[item] == v:
            k=item
    return k

for line in fileinput.input(genoTypeDat):
    lineNo = lineNo + 1
    typeSet = set()
    lociID = ""
    if lineNo == 1:
        headline = line.strip()
    else:
        vector = line.strip().split("\t")
        lineNo = lineNo + 1
        columnIndex = 0
        for field in vector:
            columnIndex = columnIndex + 1
            if columnIndex == 1:
                lociID = field.strip()
            if columnIndex > 1:
                typeSet.add(field.strip())
            # if len(typeSet) >= 3:
            #      break
        if len(typeSet) <3:
            print len(typeSet)

    lociIDDic[lociID] = typeSet

# numericLociDic = {}
# for rs in lociIDDic.keys():
#     tset = lociIDDic[rs]
#     numericDic = {}
#     i = 0
#     for e in tset:
#         numericDic[e] = i
#         i = i + 0.5
#     numericLociDic[rs] = numericDic  # 这里考虑将这个字典写到文件里
#
# #写下rsID与数值对应的文件
# logw= open(rsnumericlogfile,"w")
# logw.write("lociID\t0\t0.5\t1\n")
# for rsid in numericLociDic.keys():
#     newline=rsid
#     nlDic=numericLociDic[rsid]
#     lst=['a','b','c']
#     for num in nlDic.values():
#         if num == 0:
#             lst[0] = getkeysforv(nlDic,num)
#         if num == 0.5:
#             lst[1] = getkeysforv(nlDic, num)
#         if num == 1:
#             lst[2] = getkeysforv(nlDic, num)
#
#     newline= newline +"\t" + "\t".join(lst)+"\n"
#     logw.write(newline)
# logw.close()
#
#
# 开始写新文件
fw = open(fwfile, "w")
fw.write(headline + "\n")
linenumber= 0
for line in fileinput.input(genoTypeDat):
    linenumber = linenumber+1
    newline= ""
    if linenumber == 1:
        continue
    else:
        linearr = line.strip().split("\t")
        index = 0
        newline =linearr[0]
        numericarr=[]
        for e in linearr:
            index = index + 1
            if index > 1 and numericLociDic.has_key(linearr[0]):
                rsDic = numericLociDic[linearr[0]]
                numericarr.append(str(rsDic[e]))
        newline = newline + "\t"+"\t".join(numericarr) +"\n"
        fw.write(newline)
fw.close()

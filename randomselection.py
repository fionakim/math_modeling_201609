import os
import fileinput
import random

f = "P:\\modeling\\new_finalnumericgenotype.dat"
allcaselineDicList = []
allcontrollineDicList = []
head = ""
for line in fileinput.input(f):
    if line.startswith("sample	phenotype	rs3094315"):
        head = line.strip()
    else:
        linearr = line.strip().split("\t")
        lineDic = {}
        lineDic[linearr[0]] = line.strip()
        print linearr[1]
        if linearr[1] == '0':
            allcontrollineDicList.append(lineDic)
        if linearr[1] == '1':
            allcaselineDicList.append(lineDic)

traningcaserandomlst = random.sample(allcaselineDicList, 250)
traningcontrolrandomlst = random.sample(allcontrollineDicList, 250)


def getremaninglst(partlst, alllst):
    remainlst = []
    for e in alllst:
        if e not in partlst:
            remainlst.append(e)
    return remainlst


testcaserandomlst = getremaninglst(traningcaserandomlst, allcaselineDicList)
testcontrolrandomlst = getremaninglst(traningcontrolrandomlst, allcontrollineDicList)

trainingcontrolfile = "P:\\modeling\\new_trainingcontrolfile.dat"
testcontrolfile = "P:\\modeling\\new_testcontrolfile.dat"
traningcasefile = "P:\\modeling\\new_traningcasefile.dat"
testcasefile = "P:\\modeling\\new_testcasefile.dat"
allcontrolfile="P:\\modeling\\new_allcontrolfile.dat"
allcasefile= "P:\\modeling\\new_allcasefile.dat"

def writelist2file(lst, filew, h):
    outfile = open(filew, "w")
    outfile.write(h + "\n")
    sample_list = []
    for d in lst:
        v = "\n".join(d.values())
        sample_list.append(v)
    sample_list = [l + '\n' for l in sample_list]
    for ele in sample_list:
        outfile.write(ele)
    outfile.close()


writelist2file(traningcaserandomlst, traningcasefile, head)
writelist2file(testcaserandomlst, testcasefile, head)
writelist2file(traningcontrolrandomlst, trainingcontrolfile, head)
writelist2file(testcontrolrandomlst, testcontrolfile, head)
writelist2file(allcaselineDicList, allcasefile, head)
writelist2file(allcontrollineDicList, allcontrolfile, head)




import fileinput
from math import log
import math

from scipy.stats import chisquare

fcase = "P:\\modeling\\new\\del_new_allcasefile_transpose.dat"
fcontrol = "P:\\modeling\\new\\del_new_allcontrolfile_transpose.dat"
rsgenetypelogfile = "P:\\modeling\\rsnumericlog.dat"
snpSignificancePvaluefile = "P:\\modeling\\new\\snpSignificancePvaluefile.dat"
frecase = "P:\\modeling\\new\\new_allcasefile_transpose_frequency.dat"
frecontrol = "P:\\modeling\\new\\new_allcontrolfile_transpose_frequency.dat"


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
            if zeronumber ==0 or onenumber==0 or halfnumber ==0:
                print f + "'s " +id + " has unormal number:[0,1,0.5] " + str(zeronumber) + "," + str(onenumber)+"," +str(halfnumber)
            sd['0'] = zeronumber
            sd['1'] = onenumber
            sd['0.5'] = halfnumber
            sd['sum'] = s
        d[id] = sd

    return d


def writeDic2file(d, f):
    fw = open(f, "w")
    fw.write("rsID\tpValue\t-LOG10pValue\n")
    for k in d.keys():
        info = d[k]
        newline = k + "\t" + str(info[0]) + "\t" + str(info[1]) + "\n"
        fw.write(newline)
    fw.close()

rstypeDic = getDicfromfile(rsgenetypelogfile)
rs_case_fre_dic = ""

#
casersDic = getrateDicfromfile(fcase)
controlDic = getrateDicfromfile(fcontrol)
rs_pvalue_dic = {}
frecasedic = {}
frecontroldic ={}

for c in casersDic.keys():
    if not c.startswith('rs'):
        continue
    else:
        vcase = casersDic[c]
        vcontrol = controlDic[c]

        casearr = []
        controlarr = []
        if vcase.has_key('1') and vcase.has_key('0.5') and vcase.has_key('0'):
            casearr = [vcase['1'], vcase['0'], vcase['0.5']]
            frecasedic[c] = casearr
        else:
            print 'case: ' + c + " key error"

        if vcontrol.has_key('1') and vcontrol.has_key('0.5') and vcontrol.has_key('0') and vcontrol:
            controlarr = [vcontrol['1'], vcontrol['0'], vcontrol['0.5']]
            frecontroldic[c] = controlarr
        else:
            print 'control: ' + c + " key error"
        if len(casearr) == 3 and len(controlarr) == 3:
            ch, p = chisquare(casearr, controlarr)
            minuslogPvalue = 0.0
            if p<=0:
                print "pVALUE IS ZERO :" + c + " "+ str(p) +c + ":"+ str(vcase) + ";"+ str(vcontrol) +str(ch)
                minuslogPvalue = 400
                rs_pvalue_dic[c] = [p, minuslogPvalue]
            else:
                minuslogPvalue = -log(p, 10)
                rs_pvalue_dic[c] = [p, minuslogPvalue]
                if p < math.pow(10,-8):
                    print "strong significance: " +c + ":"+ str(vcase) + ";"+ str(vcontrol) +str(p) +" " +str(ch)

        else:
            print 'there is sth wrong in array '
writeDic2file(rs_pvalue_dic, snpSignificancePvaluefile)
# def writefreDic2file(d, f):
#     fw = open(f, "w")
#     fw.write("rsID\t1\n")
#     for k in d.keys():
#         info = d[k]
#         newline = k + "\t" + str(info[0]) + "\t" + str(info[1]) + "\n"
#         fw.write(newline)
#     fw.close()
#
#
#
# writefreDic2file(frecasedic, frecase)
# writefreDic2file(frecontroldic, frecontrol)

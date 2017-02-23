# -*- coding: utf-8 -*-
import pandas
import fileinput
import math
import numpy as np
from scipy.stats import chisquare

fcase = "P:\\modeling\\testcasefile_transpose.dat"
fcontrol = "P:\\modeling\\testcontrolfile_transpose.dat"
rsgenetypelogfile = "P:\\modeling\\rsnumericlog.dat"
ratethreshold = 0.9


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
    vdic={}
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
        expected = [preaanumber,preacnumber,preccnumber]
        observed = [aa,ac,cc]
        ch,p = chisquare(observed, expected)
        print
        vdic['chi'] = ch
        vdic['pvalue'] = p
        return vdic


rstypeDic = getDicfromfile(rsgenetypelogfile)

firstdiscardIDlst = []
chidiscardIDlist = []
qualifiedsnpIDlist = []
#
casersDic = getrateDicfromfile(fcase)
controlDic = getrateDicfromfile(fcontrol)
#

for c in casersDic.keys():
    if not c.startswith('rs'):
        continue
    else:

        vcase = casersDic[c]
        vcontrol = controlDic[c]

        vcasesum = vcase['sum']
        vcontrolsum = vcontrol['sum']

        casekeylst = vcase.keys()
        controlkeylst = vcontrol.keys()
        onecasenumber = vcase['1']
        halfcasenumber = vcase['0.5']
        zerocasenumber = vcase['0']
        onecontrolnumber = vcontrol['1']
        halfcontrolnumber = vcontrol['0.5']
        zerocontrolnumber = vcontrol['0']

        onecaserate = float(onecasenumber) / float(vcasesum)
        halfcaserate = float(halfcasenumber) / float(vcasesum)
        zerocaserate = float(zerocasenumber) / float(vcasesum)

        onecontrolrate = float(onecontrolnumber) / float(vcontrolsum)
        halfcontrolrate = float(halfcontrolnumber) / float(vcontrolsum)
        zerocontrolrate = float(zerocontrolnumber) / float(vcontrolsum)
        b1 = (onecaserate > ratethreshold) and (onecontrolrate > ratethreshold)
        b2 = (zerocaserate > ratethreshold) and (zerocontrolrate > ratethreshold)
        b3 = (halfcaserate > ratethreshold) and (halfcontrolrate > ratethreshold)
        if b1 or b2 or b3:
            firstdiscardIDlst.append(c)

        typedic = rstypeDic[c]
        onetype = typedic['1']
        halftype = typedic['0.5']
        zerotype = typedic['0']
        caseaanumber = 0
        caseacnumber = 0
        caseccnumber = 0
        controlaanumber = 0
        controlacnumber = 0
        controlccnumber = 0
        if vcase.has_key('1') and vcase.has_key('0.5') and vcase.has_key('0'):

            if len(set(onetype)) == 2:
                caseacnumber = vcase['1']
                caseaanumber = vcase['0']
                caseccnumber = vcase['0.5']
            elif len(set(halftype)) == 2:
                caseacnumber = vcase['0.5']
                caseaanumber = vcase['0']
                caseccnumber = vcase['1']
            else:
                caseacnumber = vcase['0']
                caseaanumber = vcase['0.5']
                caseccnumber = vcase['1']
        else:
            print 'case: ' + c
        if vcontrol.has_key('1') and vcontrol.has_key('0.5') and vcontrol.has_key('0'):

            if len(set(onetype)) == 2:
                controlacnumber = vcontrol['1']
                controlaanumber = vcontrol['0']
                controlccnumber = vcontrol['0.5']
            elif len(set(halftype)) == 2:
                controlacnumber = vcontrol['0.5']
                controlaanumber = vcontrol['0']
                controlccnumber = vcontrol['1']
            else:
                controlacnumber = vcontrol['0']
                controlaanumber = vcontrol['0.5']
                controlccnumber = vcontrol['1']
        else:
            print 'control: ' + c

        chicase = judgeifmatchhw(caseaanumber, caseacnumber, caseccnumber)
        chicontrol = judgeifmatchhw(controlaanumber, controlacnumber, controlccnumber)
        if chicase['pvalue'] < math.pow(10,-4) or chicontrol['pvalue'] < math.pow(10,-7):
            chidiscardIDlist.append(c)






print len(firstdiscardIDlst)
print len(chidiscardIDlist)


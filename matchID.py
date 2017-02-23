import fileinput

idfile = "P:\\modeling\\new\\chisquarefilteringsnpID_P0.005.dat"

allcasefile_transpose = "P:\\modeling\\new\\new_allcasefile_transpose.dat"
allcontrolfile_transpose = "P:\\modeling\\new\\new_allcontrolfile_transpose.dat"
testingcase_file = "P:\\modeling\\new\\new_testcasefile_transpose.dat"
trainingcase_file  = "P:\\modeling\\new\\new_traningcasefile_transpose.dat"
testingcontrol_file = "P:\\modeling\\new\\new_testingcontrol_transpose.dat"
trainingcontrol_file = "P:\\modeling\\new\\new_trainingcontrol_transpose.dat"
all_file = "P:\\modeling\\new\\new_finalnumericgenotype_transpose.dat"

select_rs_allcase_file = "P:\\modeling\\new\\allcase_chi_filtering_p005_data_transpose.dat"
select_rs_allcontrol_file = "P:\\modeling\\new\\allcontrol_chi_filtering_p005_data_transpose.dat"
select_rs_testingcase_file = "P:\\modeling\\new\\testingcase_chi_filtering_p005_data_transpose.dat"
select_rs_trainingcase_file = "P:\\modeling\\new\\trainingcase_chi_filtering_p005_data_transpose.dat"
select_rs_testingcontrol_file = "P:\\modeling\\new\\testingcontrol_chi_filtering_p005_data_transpose.dat"
select_rs_trainingcontrol_file = "P:\\modeling\\new\\trainingcontrol_chi_filtering_p005_data_transpose.dat"
select_rs_all_file = "P:\\modeling\\new\\all_chi_filtering_p005_data_transpose.dat"

def getdicfromfile(f):
    d = {}
    for line in fileinput.input(f):
        arr = line.strip().split("\t")
        d[arr[0]] = line.strip()
    return d

def writelstDic2file(lst, d, f):
    fw = open(f,"w")
    line1 = d["sample"] + "\n"
    fw.write(line1)
    line2 = d["phenotype"] + "\n"
    fw.write(line2)
    for e in lst:
        if d.has_key(e):
            fw.write(d[e] + "\n")
        else:
            print "key error"

idlst = []
for line in fileinput.input(idfile):
    arr = line.strip().split("\t")
    if not line.startswith("rsID"):
        i = arr[0]
        idlst.append(i)
allcaselineIDdic = getdicfromfile(allcasefile_transpose)
allcontrollineIDdic = getdicfromfile(allcontrolfile_transpose)
testingcaselineIDdic = getdicfromfile(testingcase_file)
trainingcaselineIDdic = getdicfromfile(trainingcase_file)
testingcontrollineIDdic = getdicfromfile(testingcontrol_file)
trainingcontrollineIDdic = getdicfromfile(trainingcontrol_file)
all_file_dic = getdicfromfile(all_file)

# writelstDic2file(idlst, allcaselineIDdic, select_rs_allcase_file)
# writelstDic2file(idlst, testingcaselineIDdic, select_rs_testingcase_file)
# writelstDic2file(idlst, trainingcaselineIDdic, select_rs_trainingcase_file)
# writelstDic2file(idlst, allcontrollineIDdic, select_rs_allcontrol_file)
# writelstDic2file(idlst, testingcontrollineIDdic, select_rs_testingcontrol_file)
# writelstDic2file(idlst, trainingcontrollineIDdic, select_rs_trainingcontrol_file)
writelstDic2file(idlst, all_file_dic, select_rs_all_file)
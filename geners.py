import fileinput
import os
rootdir = "P:\\modeling\\gene_info"
filearr = os.listdir(rootdir)
sumgeneidrs = "P:\\modeling\\new\\sum_gene_rs_id.dat"
sumrsidgene = "P:\\modeling\\new\\sum_rs_gene_id.dat"
def writeDic2file(d, f):
    fw = open(f, "w")
    fw.write("geneID\trsIDset\n")
    for k in d.keys():
        info = d[k]
        newline = k + "\t" + ",".join(info)+"\n"
        fw.write(newline)
    fw.close()
def getsetfromDic(de):
    s=set()
    for k in de.keys():
        lst = de[k]
        for i in lst:
            s.add(i.strip())
    return  s

def getDicfromfile(f):
    d = {}
    for line in fileinput.input(f):
        arr = line.strip().split("\t")
        d[arr[0]] = line.strip()
    return d

def getdicfromsetandfile(rss,sumf):
    rs_gene_dic = {}
    sumd = getDicfromfile(sumf)
    for e in rss:
        gene_for_rs_lst = []
        for t in sumd.keys():
            ta=sumd[t].split("\t")
            genecolumn = ta[1]
            linearr =genecolumn.split(",")
            if e in linearr:
                gene_for_rs_lst.append(t)
            else:
                continue
        rs_gene_dic[e] = gene_for_rs_lst
        if len(gene_for_rs_lst) ==0 :
            print e
    return rs_gene_dic

gene_rs_dic = {}
for fi in filearr:
    end = fi.find(".")
    geneID = fi[0:end]
    rslst = []
    abs_path = rootdir+"\\"+fi
    for line in fileinput.input(abs_path):
        rs = line.strip()
        rslst.append(rs)
    gene_rs_dic[geneID] = rslst

writeDic2file(gene_rs_dic,sumgeneidrs)

rsset = getsetfromDic(gene_rs_dic)
print len(rsset)

rs_gene_reference_dic = getdicfromsetandfile(rsset,sumgeneidrs)
writeDic2file(rs_gene_reference_dic,sumrsidgene)



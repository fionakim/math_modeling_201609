import Bio
import os
import shutil
import re
import fileinput
import subprocess
from Bio import SeqIO

candidateProteinSeqFile = "F:\\tmp\\at_seed_geneFromLiter_sequence.fa"
targetSpeciesFile = "F:\\tmp\\targetSpeciesList.txt"
rootProteinSeqFileDir = "F:\\tmp"


def getPathDicFromFile(targetFilePath, rootFileDir):

    pass #得到形如{“speciesA”: "protein_seq_fa_file_path"}


speciesProteinSeqFileDic = getPathDicFromFile(targetSpeciesFile, rootProteinSeqFileDir) #首先得到物种与其蛋白序列文件路径对应的字典

for currentSeqRecord in SeqIO.parse(candidateProteinSeqFile, "fasta"): # 迭代目标蛋白序列的记录

    currentSeqID= currentSeqRecord.id
    currentSeqDesc= currentSeqRecord.description
    currentSeq = currentSeqRecord.seq










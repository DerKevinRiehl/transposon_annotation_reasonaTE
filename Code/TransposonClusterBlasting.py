############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
from shutil import copyfile
import os
from os import path
from Bio import SeqIO
from shutil import rmtree

def doClusterBlasting(folderIn, folderOut, fastaFile, splitTranspLength=1000, splitParts=3, maxSeqPerFile=400):
    # Parameters
#    folderIn  = "transposonCandC"
#    folderOut = "transposonCandD"
#    fastaFile = "sequence.fasta"
#    splitTranspLength = 1000
#    splitParts = 3
#    maxSeqPerFile = 400
    
    # Copy blast transposons
    fileIn = os.path.join(folderIn,"ClusterSequences.fasta")
    copyfile(fileIn,os.path.join(folderOut,"ClusterSequences.fasta"))
    
    # Create additional blast sequences for split
    records = SeqIO.parse(os.path.join(folderOut,"ClusterSequences.fasta"), "fasta")
    f = open(os.path.join(folderOut,"ClusterSequencesFrag.fasta"),"w+")
    for r in records:
        label = r.id
        seq = str(r.seq)
        if(len(seq)>splitTranspLength):
            subSeq1 = seq[0:int(len(seq)/2)]
            subSeq2 = seq[int(len(seq)/4):int(len(seq)/4*3)]
            subSeq3 = seq[int(len(seq)/2):int(len(seq))]
            f.write(">"+label+"_frag1\n"+subSeq1+"\n")
            f.write(">"+label+"_frag2\n"+subSeq2+"\n")
            f.write(">"+label+"_frag3\n"+subSeq3+"\n")
    f.close()
    
    # Split files into parts of 400 sequences
    if(path.isdir(os.path.join(folderOut,"sequences"))):
        rmtree(os.path.join(folderOut,"sequences"))
    os.mkdir(os.path.join(folderOut,"sequences"))
    if(path.isdir(os.path.join(folderOut,"fragments"))):
        rmtree(os.path.join(folderOut,"fragments"))
    os.mkdir(os.path.join(folderOut,"fragments"))
    records = SeqIO.parse(os.path.join(folderOut,"ClusterSequences.fasta"), "fasta")
    counter = 1
    fileCounter = 1
    fw = open(os.path.join(folderOut,"sequences","Sequence"+str(fileCounter)+".fasta"),"w+")
    for r in records:
        if(counter<=maxSeqPerFile):
            counter += 1
            fw.write(">"+r.id+"\n"+str(r.seq)+"\n")
        else:
            counter = 0
            fileCounter += 1
            fw.close()
            fw = open(os.path.join(folderOut,"sequences","Sequence"+str(fileCounter)+".fasta"),"w+")
            fw.write(">"+r.id+"\n"+str(r.seq)+"\n")
    fw.close()
    records = SeqIO.parse(os.path.join(folderOut,"ClusterSequencesFrag.fasta"), "fasta")
    counter = 1
    fileCounter = 1
    fw = open(os.path.join(folderOut,"fragments","Fragment"+str(fileCounter)+".fasta"),"w+")
    for r in records:
        if(counter<=maxSeqPerFile):
            counter += 1
            fw.write(">"+r.id+"\n"+str(r.seq)+"\n")
        else:
            counter = 0
            fileCounter += 1
            fw.close()
            fw = open(os.path.join(folderOut,"fragments","Fragment"+str(fileCounter)+".fasta"),"w+")
            fw.write(">"+r.id+"\n"+str(r.seq)+"\n")
    fw.close()
    
    # Create BlastDBs
    filesS = os.listdir(os.path.join(folderOut,"sequences"))
    filesF = os.listdir(os.path.join(folderOut,"fragments"))
    for file in filesS:
        os.system("cd "+os.path.join(folderOut,"sequences")+" && makeblastdb -in "+file+" -parse_seqids -title "+file+" -dbtype nucl")
    for file in filesF:
        os.system("cd "+os.path.join(folderOut,"fragments")+" && makeblastdb -in "+file+" -parse_seqids -title "+file+" -dbtype nucl")
    
    # Run BLAST
    filesS = os.listdir(os.path.join(folderOut,"sequences"))
    filesF = os.listdir(os.path.join(folderOut,"fragments"))
    filesS2 = list()
    filesF2 = list()
    for f in filesS:
        if(f.endswith(".fasta")):
            filesS2.append(f)
    for f in filesF:
        if(f.endswith(".fasta")):
            filesF2.append(f)
    filesS2.sort()
    filesF2.sort()
    for file in filesS2:
        print("Blast ",file," / ",len(filesS2))
        os.system("blastn -db "+os.path.join(folderOut,"sequences",file)+" -out "+os.path.join(folderOut,"SequenceBlastResults_"+file+".txt")+" -outfmt \"7 sacc stitle qframe evalue bitscore qstart qend qlen sstart send\" -query "+fastaFile+" -evalue 0.1 -num_threads 10")
#    for file in filesF2:
#        print("Blast ",file," / ",len(filesF2))
#        os.system("blastn -db "+os.path.join(folderOut,"fragments",file)+" -out "+os.path.join(folderOut,"FragmentBlastResults_"+file+".txt")+" -outfmt \"7 sacc stitle qframe evalue bitscore qstart qend qlen sstart send\" -query "+fastaFile+" -evalue 0.1 -num_threads 10")

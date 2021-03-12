############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
import os

# Methods
def countAnnotationsBP1(seq,annotations,clss):
    if(seq not in annotations):
        return 0
    numTP = 0  
    for ann in annotations[seq]:
        if(clss=="all"):
            numTP += 1
        else:
            if(ann[2].startswith(clss)):
                numTP += 1
    return numTP
def countAnnotationsBP2(seq,annotations,clss):
    if(seq not in annotations):
        return 0
    sumBP = 0
    for ann in annotations[seq]:
        if(clss=="all"):
            sumBP += abs(ann[1]-ann[0])+1
        else:
            if(ann[2].startswith(clss)):
                sumBP += abs(ann[1]-ann[0])+1
    return sumBP
def countAnnotationsBP3(seq,annotations):
    if(seq not in annotations):
        return 0
    sumBP = 0
    for ann in annotations[seq]:
        sumBP += abs(ann[1]-ann[0])+1
    return sumBP
def createStatistics1(sequenceHeadsFile,fileTransposonAnn,fileW):
#    sequenceHeadsFile = "parsedAnnotations/sequence_heads.txt"
#    fileTransposonAnn = "finalResults/FinalAnnotations_Transposons.gff3"    
    # RFSB classes
    classes = "1,1/1,1/1/1,1/1/2,1/1/3,1/2,1/2/1,1/2/2,2,2/1,2/1/1,2/1/2,2/1/3,2/1/4,2/1/5,2/1/6,2/2,2/3".split(",")
    # load sequence heads
    heads = {}
    f = open(sequenceHeadsFile,"r")
    line = f.readline()
    while line!="":
        parts = line.replace("\n","").replace(">","").split("\t")
        heads[parts[0]] = parts[1]
        line = f.readline()
    f.close()
    # load annotations
    annotations = {}
    annotations["all"] = list()
    f = open(fileTransposonAnn,"r")
    line = f.readline()
    while line!="":
        parts = line.split("\t")
        start = int(parts[3])
        end   = int(parts[4])
        clss  = parts[8].split(";")[1].split("=")[1].replace("\n","").split("(")[0]
        chrom = parts[0]
        if(not chrom in annotations):
            annotations[chrom] = list()
        annotations["all"].append([start,end,clss])
        annotations[chrom].append([start,end,clss])
        line = f.readline()
    f.close()
    ## analyse and print results
    fileW.write("SeqID"+"\t"+"SeqName"+"\t"+"#Num_transposons by classes")
    fileW.write("\n")
    print("SeqID","\t","SeqName","\t","#Num_transposons by classes")
    fileW.write("SeqID"+"\t"+"SeqName"+"\t"+"all"+"\t"+"\t".join(classes))
    fileW.write("\n")
    print("SeqID","\t","SeqName","\t","all","\t","\t".join(classes))
    fileW.write("all"+"\t"+"all"+"\t"+str(countAnnotationsBP1("all",annotations,"all")))
    fileW.write("\t")
    print("all","\t","all","\t",countAnnotationsBP1("all",annotations,"all"),end="\t")
    for c in classes:
        fileW.write(str(countAnnotationsBP1("all",annotations,c)))
        fileW.write("\t")
        print(countAnnotationsBP1("all",annotations,c),end="\t")
    fileW.write("\n")
    print("")
    for seq in list(heads.keys()):
        fileW.write(seq+"\t"+heads[seq]+"\t"+str(countAnnotationsBP1(seq,annotations,"all")))
        fileW.write("\t")
        print(seq,"\t",heads[seq],"\t",countAnnotationsBP1(seq,annotations,"all"),end="\t")
        for c in classes:
            fileW.write(str(countAnnotationsBP1(seq,annotations,c)))
            fileW.write("\t")
            print(countAnnotationsBP1(seq,annotations,c),end="\t")
        fileW.write("\n")
        print("")
def createStatistics2(sequenceHeadsFile,fileTransposonAnn,fileW):
#    sequenceHeadsFile = "parsedAnnotations/sequence_heads.txt"
#    fileTransposonAnn = "finalResults/FinalAnnotations_Transposons.gff3"    
    # RFSB classes
    classes = "1,1/1,1/1/1,1/1/2,1/1/3,1/2,1/2/1,1/2/2,2,2/1,2/1/1,2/1/2,2/1/3,2/1/4,2/1/5,2/1/6,2/2,2/3".split(",")
    # load sequence heads
    heads = {}
    f = open(sequenceHeadsFile,"r")
    line = f.readline()
    while line!="":
        parts = line.replace("\n","").replace(">","").split("\t")
        heads[parts[0]] = parts[1]
        line = f.readline()
    f.close()
    # load annotations
    annotations = {}
    annotations["all"] = list()
    f = open(fileTransposonAnn,"r")
    line = f.readline()
    while line!="":
        parts = line.split("\t")
        start = int(parts[3])
        end   = int(parts[4])
        clss  = parts[8].split(";")[1].split("=")[1].replace("\n","").split("(")[0]
        chrom = parts[0]
        if(not chrom in annotations):
            annotations[chrom] = list()
        annotations["all"].append([start,end,clss])
        annotations[chrom].append([start,end,clss])
        line = f.readline()
    f.close()
    ## analyse and print results
    fileW.write("SeqID"+"\t"+"SeqName"+"\t"+"#BP_transposons by classes")
    fileW.write("\n")
    print("SeqID","\t","SeqName","\t","#BP_transposons by classes")
    fileW.write("SeqID"+"\t"+"SeqName"+"\t"+"all"+"\t"+"\t".join(classes))
    fileW.write("\n")
    print("SeqID","\t","SeqName","\t","all","\t","\t".join(classes))
    fileW.write("all"+"\t"+"all"+"\t"+str(countAnnotationsBP2("all",annotations,"all")))
    fileW.write("\t")
    print("all","\t","all","\t",countAnnotationsBP2("all",annotations,"all"),end="\t")
    for c in classes:
        fileW.write(str(countAnnotationsBP2("all",annotations,c)))
        fileW.write("\t")
        print(countAnnotationsBP2("all",annotations,c),end="\t")
    fileW.write("\n")
    print("")
    for seq in list(heads.keys()):
        fileW.write(seq+"\t"+heads[seq]+"\t"+str(countAnnotationsBP2(seq,annotations,"all")))
        fileW.write("\t")
        print(seq,"\t",heads[seq],"\t",countAnnotationsBP2(seq,annotations,"all"),end="\t")
        for c in classes:
            fileW.write(str(countAnnotationsBP2(seq,annotations,c)))
            fileW.write("\t")
            print(countAnnotationsBP2(seq,annotations,c),end="\t")
        fileW.write("\n")
        print("")
def createStatistics3(sequenceHeadsFile,fileTransposonMask,fileW):
#    sequenceHeadsFile  = "parsedAnnotations/sequence_heads.txt"
#    fileTransposonMask = "finalResults/FinalAnnotations_TransposonMask.gff3"    
    # load sequence heads
    heads = {}
    f = open(sequenceHeadsFile,"r")
    line = f.readline()
    while line!="":
        parts = line.replace("\n","").replace(">","").split("\t")
        heads[parts[0]] = parts[1]
        line = f.readline()
    f.close()
    # load annotations
    annotations = {}
    annotations["all"] = list()
    f = open(fileTransposonMask,"r")
    line = f.readline()
    while line!="":
        start = int(line.split("\t")[3])
        end   = int(line.split("\t")[4])
        chrom = line.split("\t")[0]
        if(not chrom in annotations):
            annotations[chrom] = list()
        annotations["all"].append([start,end])
        annotations[chrom].append([start,end])
        line = f.readline()
    f.close()
    # analyse and print results
    fileW.write("SeqID"+"\t"+"SeqName"+"\t"+"#BP_transposons")
    fileW.write("\n")
    print("SeqID","\t","SeqName","\t","#BP_transposons")
    fileW.write("all"+"\t"+"all"+str(countAnnotationsBP3("all",annotations)))
    fileW.write("\n")
    print("all","\t","all",countAnnotationsBP3("all",annotations))
    for seq in list(heads.keys()):
        fileW.write(seq+"\t"+heads[seq]+"\t"+str(countAnnotationsBP3(seq,annotations)))
        fileW.write("\n")
        print(seq,"\t",heads[seq],"\t",countAnnotationsBP3(seq,annotations))
def createFinalStatistics(sequenceHeadsFile,fileTransposonAnn,fileTransposonMask,fileOutput):
    fileW = open(fileOutput,"w+")
    createStatistics1(sequenceHeadsFile,fileTransposonAnn,fileW)
    fileW.write("\n\n")
    print("\n")
    createStatistics2(sequenceHeadsFile,fileTransposonAnn,fileW)
    fileW.write("\n\n")
    print("\n")
    createStatistics3(sequenceHeadsFile,fileTransposonMask,fileW)
    fileW.close()
    
def countAnns(file):
    ctr = 0
    f = open(file,"r")
    line = f.readline()
    while line!="":
        if(not line.startswith("#")):
            ctr += 1
        line = f.readline()
    f.close()
    return ctr
    
def createParsedAnnotationStatistics(parsedAnnFolder):
    files = os.listdir(parsedAnnFolder)
    files.sort()
    for file in files:
        if(file.endswith(".gff3")):
            print(countAnns(os.path.join(parsedAnnFolder,file)),"\t",file)
  
#createParsedAnnotationStatistics("parsedAnnotations")

#import sys
#args = sys.argv
#createParsedAnnotationStatistics(args[1])
#createStatistics("parsedAnnotations/sequence_heads.txt","finalResults/ToolAnnotations_Transposons.gff3","finalResults/ToolAnnotations_TransposonMask.gff3","ResultTool.txt")
#createStatistics("parsedAnnotations/sequence_heads.txt","finalResults/FinalAnnotations_Transposons.gff3","finalResults/FinalAnnotations_TransposonMask.gff3","ResultFinal.txt")

#python3 Statistics.py projects_RHIZIPHAGUS_IRR/Japan/parsedAnnotations/
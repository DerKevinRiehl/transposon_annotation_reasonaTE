############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
import os
from Bio import SeqIO
from operator import itemgetter
from GFFTools import generateMask

# Methods
def correctPositions(start, end):
    start = int(start)
    end = int(end)
    if(start>end):
        c = end
        end = start
        start = c
#    start = str(start)
#    end = str(end)
    return start, end

def getChromosomeLengths(genomeFile):
    lengths = {}
    f = open(genomeFile,"r")
    line = f.readline().replace("\n","")
    seq = 0
    chrom = ""
    while line!="":
        if(line.startswith(">")):
            if(chrom!=""):
                lengths[chrom] = seq
            chrom = line.replace(">","").replace("\n","")
            seq = 0
        else:        
            seq += len(line.replace("\n",""))
        line = f.readline().replace("\n","")
    f.close()
    if(chrom!=""):
        lengths[chrom] = seq
    return lengths

def parseBlastOutput(lengths, folderIn, suffix, folderOut, resultFileOut):
    transposonAnnotations = {}
    files = os.listdir(folderIn)
    for file in files:
        if(file.startswith(suffix)):
            print(file)
            f = open(os.path.join(folderIn,file),"r")
            line = f.readline()
            chrom = f.readline().replace("\n","").split(":")[1].rstrip().lstrip()
            if(not chrom in transposonAnnotations):
                transposonAnnotations[chrom] = list()
            line = f.readline()
            while True:
                while line!="" and line.startswith("#") and not line.startswith("# BLAST processed") and not line.startswith("# RPSTBLASTN"):
                    line = f.readline()
                while line!="" and not line.startswith("#"):
                    transposons = line.replace("\n","").split("\t")
                    annoSoftware = "blastn" 
                    features     = transposons[0]
                    attributes = ""
                    transpNr = "."
                    if(int(transposons[5])<int(transposons[6])):
                        strand = "+"
                    else:
                        strand = "-"
                    start  = int(transposons[5])
                    end    = int(transposons[6])
                    start,end = correctPositions(start, end)
                    score  = "."
                    phase  = "."
                    if(start<lengths[chrom] and end<lengths[chrom]):
                        transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
                    line = f.readline()
                if(line=="" or line.startswith("# BLAST processed")):
                    break
                else:
                    line = f.readline()
                    # print(line)
                    chrom = line.replace("\n","").split(":")[1].rstrip().lstrip().split("=")[0].replace("length","").lstrip().rstrip().replace(">","")
                    if(not chrom in transposonAnnotations):
                        transposonAnnotations[chrom] = list()
                    line = f.readline()
            f.close()
    # Print results to Annotation file
    f = open(os.path.join(folderOut, resultFileOut),"w+")
    keys = list(transposonAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()
    
def getAnnotation(parts):
    start = int(parts[3])
    end   = int(parts[4])
    if(start>=end):
        c = start
        start = end
        end = c
    return [parts[0],start,end,parts[2],end-start,parts[6]]
    # chrom, start, end, cluster, length, strand
    #   0  ,   1  ,  2 ,    3   ,   4   ,   5

def loadAnnotations(fileGFF):
    print(fileGFF,"\t",end="")
    annotationsTransposons = {}
    file = open(fileGFF,"r")
    line = file.readline().replace("\n","")
    while line!="":
        if(line.startswith("#")):
            line = file.readline().replace("\n","")
            continue
        parts = line.split("\t")
        annotation = getAnnotation(parts)
        chrom = annotation[0]
        if(not chrom in annotationsTransposons):
            annotationsTransposons[chrom] = list()
        annotationsTransposons[chrom].append(annotation)
        line = file.readline().replace("\n","")
    file.close()
    return annotationsTransposons

def loadAnnotationsB(fileGFF):
    print(fileGFF,"\t",end="")
    annotationsTransposons = {}
    file = open(fileGFF,"r")
    line = file.readline().replace("\n","")
    while line!="":
        if(line.startswith("#")):
            line = file.readline().replace("\n","")
            continue
        parts = line.split("\t")
        annotation = getAnnotation(parts)
        chrom = annotation[0]
        if(not chrom in annotationsTransposons):
            annotationsTransposons[chrom] = {}
        if(not annotation[3] in annotationsTransposons[chrom]):
            annotationsTransposons[chrom][annotation[3]] = list()
        annotationsTransposons[chrom][annotation[3]].append(annotation)
        line = file.readline().replace("\n","")
    file.close()
    return annotationsTransposons




def doAnalysis(genomeFile, folderInCB,folderInCC,folderInCD,folderOut, maskMergeTol=5, lenTol=0.1, tolBP=0.1):
#    folderInCB = "transposonCandB"
#    folderInCC = "transposonCandC"
#    folderInCD = "transposonCandD"
#    folderOut  = "transposonCandE"
#    maskMergeTol = 5
#    lenTol       = 0.1
#    tolBP        = 0.1
    
    ## ###############################################################################################################################
    ## Filter Full Annotations (they can only occur if length is similar to full transposons) and save results
    ## ###############################################################################################################################
    # Get Chromosome Lengths
    lengths = getChromosomeLengths(genomeFile)
    # Parse BlastOutputs
    parseBlastOutput(lengths,os.path.join(folderInCD), "SequenceBlastResults_Sequence", folderOut, "SequenceBlastResults.gff3")
    #    parseBlastOutput(os.path.join(folderInCD), "FragmentBlastResults_Fragment", folderOut, "FragmentBlastResults.gff3")
    # Load Full Annotations
    annotationsTransposons = loadAnnotations(os.path.join(folderOut,"SequenceBlastResults.gff3"))
    # Load Cluster Sequence Length
    records = SeqIO.parse(os.path.join(folderInCD,"ClusterSequences.fasta"), "fasta")
    clusterLenght = {}
    for r in records:
        clusterLenght[r.id] = len(r.seq)  
    # Filter Full Annotations (they can only occur if length is similar to full transposons) and save results
    f = open(os.path.join(folderOut,"SequenceBlastResults_filtered.gff3"),"w")
    for chrom in annotationsTransposons:
        for ann in annotationsTransposons[chrom]:
            if(abs(ann[4]-clusterLenght[ann[3]])<lenTol*clusterLenght[ann[3]]):
                f.write(chrom+"\t"+"blastn"+"\t"+ann[3]+"\t"+str(ann[1])+"\t"+str(ann[2])+"\t"+"."+"\t"+ann[5]+"\t"+"."+"\t"+"."+"\t")
                f.write("\n")
    f.close()
    
    ## ###############################################################################################################################
    ## Filter Full Annotations that are intersecting with mask from tool annotations
    ## ###############################################################################################################################
    # Generate tool annotations
    files = os.listdir(folderInCB)
    filesN = list()
    for f in files:
        if(f.startswith("CandidatesB_it")):
            filesN.append(int(f.replace("CandidatesB_it","").replace(".gff3","")))
    filesN.sort()
    generateMask(os.path.join(folderInCB,"CandidatesB_it"+str(filesN[-1])+".gff3"), os.path.join(folderOut,"MaskTools.gff3"), "transposon")
    # Load Full Annotations again
    annotationsTransposons = loadAnnotations(os.path.join(folderOut,"SequenceBlastResults_filtered.gff3"))
    # Load tool annotations
    annotationsMask = loadAnnotations(os.path.join(folderOut,"MaskTools.gff3"))
    # Filter Full Annotations that are intersecting with mask from tool annotations
    for chrom in list(annotationsTransposons.keys()):
        if(chrom in annotationsMask):
            print(chrom,"\tIteration","\t",len(annotationsTransposons[chrom]))
            annotationsMask[chrom] = sorted(annotationsMask[chrom], key=itemgetter(0))
            delIndex = list()
            # Search for simple duplicates
            numAnns = len(annotationsTransposons[chrom])
            for ix1 in range(0,numAnns):
                if(ix1%1000==0):
                    print(chrom,"\t",ix1,"\t",numAnns)
                annT = annotationsTransposons[chrom][ix1]
                startIx2 = 0
                for ix2 in range(0,len(annotationsMask[chrom])):
                    annM = annotationsMask[chrom][ix2]
                    if(annM[2]>=annT[1]):
                        startIx2 = ix2
                        break
                    if(annT[2]<annM[1]):
                        break
                for ix2 in range(startIx2,len(annotationsMask[chrom])):
                    annM = annotationsMask[chrom][ix2]
                    if(annT[2]<annM[1]):
                        break
                    if((annM[1]+maskMergeTol>=annT[1] and annM[1]<=annT[2]-maskMergeTol) or (annT[2]<=annM[2]+maskMergeTol and annT[2]>=annM[1]-maskMergeTol)):
                        delIndex.append(ix1)
            # Delete them
            delIndex = list(set(delIndex))
            for index in sorted(delIndex, reverse=True):
                del annotationsTransposons[chrom][index]
            print(chrom,"\tIteration","\t",len(annotationsTransposons[chrom]))
    f = open(os.path.join(folderOut,"SequenceBlastResults_filtered2.gff3"),"w")
    for chrom in list(annotationsTransposons.keys()):
        for ann in annotationsTransposons[chrom]:
            f.write(chrom+"\t"+"blastn"+"\t"+ann[3]+"\t"+str(ann[1])+"\t"+str(ann[2])+"\t"+"."+"\t"+ann[5]+"\t"+"."+"\t"+"."+"\t")
            f.write("\n")
    f.close()
    
    ## ###############################################################################################################################
    ## Filter Full Annotations they can only pass if their length is similar to the cluster
    ## ###############################################################################################################################
    # Load Full Annotations again
    annotationsTransposons = loadAnnotations(os.path.join(folderOut,"SequenceBlastResults_filtered2.gff3"))
    # Load Cluster Sequence Lengths
    clusterLengths = {}
    f = open(os.path.join(folderInCC,"ClusterSequences.fasta"),"r")
    line = f.readline()
    while line!="":
        if(line.startswith(">")):
            name = line.replace(">","").replace("\n","")
            line = f.readline().replace("\n","")
            clusterLengths[name] = len(line)
        line = f.readline()
    f.close()
    # Filter Full Annotations they can only pass if their length is similar to the cluster
    f = open(os.path.join(folderOut,"SequenceBlastResults_filtered3.gff3"),"w+")
    for chrom in list(annotationsTransposons.keys()):
        for ann in annotationsTransposons[chrom]:
            if(int(float(ann[4])*(1.0+lenTol))>=clusterLengths[ann[3]]):
                f.write(chrom+"\t"+"blastn"+"\t"+ann[3]+"\t"+str(ann[1])+"\t"+str(ann[2])+"\t"+"."+"\t"+ann[5]+"\t"+"."+"\t"+"."+"\t")
                f.write("\n")
    f.close()
    
    ## ###############################################################################################################################
    ## Filter Full Annotations (for each cluster they can only occur once at the same position) and save results
    ## ###############################################################################################################################
    # Load Full Annotations again
    annotationsTransposons = loadAnnotationsB(os.path.join(folderOut,"SequenceBlastResults_filtered3.gff3"))
    # Filter Full Annotations (for each cluster they can only occur once at the same position) and save results
    for chrom in list(annotationsTransposons.keys()):
        for cluster in list(annotationsTransposons[chrom].keys()):
            iteration = 1
            while True:
                print(chrom,"\t",cluster,"\tIteration\t",iteration,"\t",len(annotationsTransposons[chrom][cluster]))
                delIndex = list()
                # Search for simple duplicates
                numAnns = len(annotationsTransposons[chrom][cluster])
                for ix1 in range(0,numAnns):
                    if(ix1%1000==0):
                        print(ix1,"/",numAnns,"...")
                    for ix2 in range(0,numAnns):
                        if(ix1!=ix2 and not(ix2 in delIndex) and not(ix1 in delIndex)):
                            ann1 = annotationsTransposons[chrom][cluster][ix1]
                            ann2 = annotationsTransposons[chrom][cluster][ix2]
                            if(abs(ann2[1]-ann1[1])<tolBP*ann1[4]  and  abs(ann2[2]-ann1[2])<tolBP*ann1[4]):
                                delIndex.append(ix2)
                # Delete them
                delIndex = list(set(delIndex))
                if(len(delIndex)==0):
                    break
                for index in sorted(delIndex, reverse=True):
                    del annotationsTransposons[chrom][cluster][index]
                iteration += 1
    f = open(os.path.join(folderOut,"SequenceBlastResults_filtered4.gff3"),"w")
    for chrom in list(annotationsTransposons.keys()):
        for cluster in list(annotationsTransposons[chrom].keys()):
            for ann in annotationsTransposons[chrom][cluster]:
                f.write(chrom+"\t"+"blastn"+"\t"+ann[3]+"\t"+str(ann[1])+"\t"+str(ann[2])+"\t"+"."+"\t"+ann[5]+"\t"+"."+"\t"+"."+"\t")
                f.write("\n")
    f.close()
    
    ## ###############################################################################################################################
    ## Merge annotations (of different clusters this time) if at same position (overlap 95%)
    ## ###############################################################################################################################
    # Load Filtered Annotations2
    annotationsTransposons = loadAnnotations(os.path.join(folderOut,"SequenceBlastResults_filtered4.gff3"))
    # Merge annotaitons (of different clusters this time) if at same position (overlap 95%)
    for chrom in list(annotationsTransposons.keys()):
        iteration = 1
        while True:
            annotationsTransposons[chrom] = sorted(annotationsTransposons[chrom], key=itemgetter(1))
            print(chrom,"\tIteration\t",iteration,"\t",len(annotationsTransposons[chrom]))
            delIndex = list()
            # Search for simple duplicates
            numAnns = len(annotationsTransposons[chrom])
            for ix1 in range(0,numAnns):
                if(ix1%1000==0):
                    print(ix1,"/",numAnns,"...")
                if(ix1 not in delIndex):
                    for ix2 in range(max(0,ix1-100),min(numAnns,ix1+100)):
                        if(ix1!=ix2):
                            if(not(ix2 in delIndex)):
                                ann1 = annotationsTransposons[chrom][ix1]
                                ann2 = annotationsTransposons[chrom][ix2]
                                if(abs(ann2[1]-ann1[1])<tolBP*ann1[4]  and  abs(ann2[2]-ann1[2])<tolBP*ann1[4]):
                                    delIndex.append(ix2)
            # Delete them
            delIndex = list(set(delIndex))
            if(len(delIndex)==0):
                break
            for index in sorted(delIndex, reverse=True):
                del annotationsTransposons[chrom][index]
            iteration += 1
    f = open(os.path.join(folderOut,"SequenceBlastResults_filtered5.gff3"),"w")
    for chrom in list(annotationsTransposons.keys()):
        for ann in annotationsTransposons[chrom]:
            f.write(chrom+"\t"+"blastn"+"\t"+ann[3]+"\t"+str(ann[1])+"\t"+str(ann[2])+"\t"+"."+"\t"+ann[5]+"\t"+"."+"\t"+"."+"\t")
            f.write("\n")
    f.close()

    ## ###############################################################################################################################
    ## Merge annotations (of different clusters this time) if same length and almost same position (overlap 95%)
    ## ###############################################################################################################################
    # Load Filtered Annotations2
    annotationsTransposons = loadAnnotations(os.path.join(folderOut,"SequenceBlastResults_filtered5.gff3"))
    # Merge annotaitons (of different clusters this time) if at same position (overlap 95%)
    for chrom in list(annotationsTransposons.keys()):
        iteration = 1
        while True:
            annotationsTransposons[chrom] = sorted(annotationsTransposons[chrom], key=itemgetter(1))
            print(chrom,"\tIteration\t",iteration,"\t",len(annotationsTransposons[chrom]))
            delIndex = list()
            # Search for simple duplicates
            numAnns = len(annotationsTransposons[chrom])
            for ix1 in range(0,numAnns):
                if(ix1%1000==0):
                    print(ix1,"/",numAnns,"...")
                if(ix1 not in delIndex):
                    for ix2 in range(max(0,ix1-100),min(numAnns,ix1+100)):
                        if(ix1!=ix2):
                            if(not(ix2 in delIndex)):
                                ann1 = annotationsTransposons[chrom][ix1]
                                ann2 = annotationsTransposons[chrom][ix2]
                                if(abs(ann1[4]-ann2[4])<ann1[4]*4*lenTol):
                                    if(abs(ann2[1]-ann1[1])<(1-tolBP)*ann1[4]  and  abs(ann2[2]-ann1[2])<(1-tolBP)*ann1[4]):
                                        delIndex.append(ix2)
            # Delete them
            delIndex = list(set(delIndex))
            if(len(delIndex)==0):
                break
            for index in sorted(delIndex, reverse=True):
                del annotationsTransposons[chrom][index]
            iteration += 1
    f = open(os.path.join(folderOut,"SequenceBlastResults_filtered6.gff3"),"w")
    for chrom in list(annotationsTransposons.keys()):
        for ann in annotationsTransposons[chrom]:
            f.write(chrom+"\t"+"blastn"+"\t"+ann[3]+"\t"+str(ann[1])+"\t"+str(ann[2])+"\t"+"."+"\t"+ann[5]+"\t"+"."+"\t"+"."+"\t")
            f.write("\n")
    f.close()


############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
from operator import itemgetter
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

# Methods
def generateMask(fileIn, fileOut, filterSuffix, maskMergeTol=10):
    # Parameters
#    fileIn = "transposonCandB/CandidatesB_it16.gff3"
#    fileOut = "MaskResult.gff3"
#    filterSuffix = "transposon"
#    maskMergeTol = 10
    # Load Annotation File
    annotations = {}
    f = open(fileIn,"r")
    line = f.readline()
    while line!="":
        if(not line.startswith("#")):
            parts = line.split("\t")
            if(parts[2]==filterSuffix):
                chrom = parts[0]
                if(not chrom in annotations):
                    annotations[chrom] = list()
                start = int(parts[3])
                end   = int(parts[4])
                if(end<start):
                    c = start
                    start = end
                    end = c
                annotations[chrom].append([start,end])
        line = f.readline()
    f.close() 
    # Sort Annotations by start point
    for chrom in list(annotations.keys()):
        annotations[chrom] = sorted(annotations[chrom], key=itemgetter(0))
    # Merge Annotations
    for chrom in list(annotations.keys()):
        iteration = 0
        while True:
            annotations[chrom] = sorted(annotations[chrom], key=itemgetter(0))
            coupleL = list()
            foundXL = list()
            annsNo = len(annotations[chrom])
            iteration += 1
            for x in range(0,annsNo):
                for y in range(x,annsNo):
                    if(x!=y):
                        annX = annotations[chrom][x]
                        annY = annotations[chrom][y]
                        
                        if(annX[1]<annY[0]-maskMergeTol):
                            break
                        else:
                            if((not x in foundXL) and (not y in foundXL)):
                                if((annY[0]+maskMergeTol>=annX[0] and annY[0]<=annX[1]-maskMergeTol) or (annX[1]<=annY[1]+maskMergeTol and annX[1]>=annY[0]-maskMergeTol)):
                                    foundXL.append(x)
                                    foundXL.append(y)
                                    coupleL.append([min(annX[0],annY[0]),max(annX[1],annY[1])])
            if(len(foundXL)==0):
                break
            else:
                foundXL = list(set(foundXL))
                for index in sorted(foundXL, reverse=True):
                    del annotations[chrom][index]
                for ann in coupleL:
                    annotations[chrom].append(ann)
    # Sort Annotations by start point
    for chrom in list(annotations.keys()):
        annotations[chrom] = sorted(annotations[chrom], key=itemgetter(0))
    # Save Annotations
    f = open(fileOut,"w+")
    for chrom in list(annotations.keys()):
        for ann in annotations[chrom]:
            f.write(chrom+"\t"+"mask"+"\t"+"mask"+"\t"+str(ann[0])+"\t"+str(ann[1])+"\t"+"."+"\t"+"+"+"\t"+"."+"\t"+"."+"\n")
    f.close()

def mergeAnnotations(files, fileOut):
    fW = open(fileOut,"w+")
    for fileIn in files:
        f = open(fileIn,"r")
        line = f.readline()
        while line!="":
            if(not line.startswith("#")):
                fW.write(line)
            line = f.readline()
        f.close() 
    fW.close()
    
def applyMask(maskIn, fileIn, fileOut):
    # Parameters
#    maskIn = "finalResults/ToolAnnotations_TransposonMask.gff3"
#    fileIn = "parsedAnnotations/transposonPSI.gff3"
#    fileOut = "finalResults/Proteins.gff3"
    # Load Annotation File
    maskAnnotations = {}
    f = open(maskIn,"r")
    line = f.readline()
    while line!="":
        if(not line.startswith("#")):
            parts = line.split("\t")
            chrom = parts[0]
            if(not chrom in maskAnnotations):
                maskAnnotations[chrom] = list()
            start = int(parts[3])
            end   = int(parts[4])
            if(end<start):
                c = start
                start = end
                end = c
            maskAnnotations[chrom].append([start,end])
        line = f.readline()
    f.close() 
    # Sort Annotations by start point
    for chrom in list(maskAnnotations.keys()):
        maskAnnotations[chrom] = sorted(maskAnnotations[chrom], key=itemgetter(0))
    # Apply Mask
    fW = open(fileOut,"w+")
    fR = open(fileIn,"r")
    line = fR.readline()
    while line!="":
        if(line.startswith("#")):
            fW.write(line)
        else:
            parts = line.split("\t")
            chrom = parts[0]
            start = int(parts[3])
            end   = int(parts[4])
            if(start>end):
                c = start
                start = end
                end = c
            if chrom in maskAnnotations:
                for x in range(0,len(maskAnnotations[chrom])):
                    if(maskAnnotations[chrom][x][0]>end):
                        break
                    else:
                        maskStart = maskAnnotations[chrom][x][0]
                        maskEnd   = maskAnnotations[chrom][x][1]
                        if((start>=maskStart and start<=maskEnd) or (end>=maskStart and end<=maskEnd)):
                            fW.write(line)
        line = fR.readline()
    fR.close()
    fW.close()
    
def cleanLine(line):
    parts = line.replace("\n","").split("\t")
    parts[8] = "transposon="+parts[8]
    if(len(parts)>9):
        parts[8] = parts[8]+" ;description="+";".join(parts[9:])
    parts = parts[:9]
    return "\t".join(parts)

def cleanLine2(line):
    parts = line.replace("\n","").split("\t")
    parts[8] = ""
    if(len(parts)>9):
        parts[8] = "description="+";".join(parts[9:])
    parts = parts[:9]
    return "\t".join(parts)

def getAnnotation(parts):
    if(int(parts[3])>=int(parts[4])):
        print("ERROR")
    return [parts[0],int(parts[3]),int(parts[4]),parts[2],parts[6],parts[8].split(";")[0].split("transposon=")[1]]
    # chrom, start, end, type  , strand, transpNr
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

def extractFastaFromGFF(gffFile, fastaSource, fastaTarget):
#    gffFile = "finalResults/ToolAnnotations_Transposons.gff3"
#    fastaSource = "parsedAnnotations/sequence.fasta"
#    fastaTarget = "finalResults/ToolAnnotations_TransposonSequences.fasta"
    annotationsTransposons = loadAnnotations(gffFile)
    f = open(fastaTarget,"w+")
    records = SeqIO.parse(fastaSource, "fasta")
    for r in records:
        for chromosome in annotationsTransposons:
            if(r.id==chromosome):
                for ann in annotationsTransposons[chromosome]:
                    start  = int(ann[1])
                    end    = int(ann[2])
                    if(end==start):
                        end += 1
                    elif(start>end):
                        c = end
                        end = start
                        start = c
                    start  = str(start)
                    end    = str(end)
                    strand = ann[4]
                    strandNum = 1
                    if(strand=="-"):
                        strandNum = -1
                    transposonFeature = SeqFeature(FeatureLocation(int(start), int(end), strand=strandNum), type="CDS")
                    extractedSequence  = transposonFeature.extract(r.seq)
                    f.write(">Transposon"+ann[5]+"\t"+chromosome+":"+"["+str(strand)+"]"+str(start)+"-"+str(end)+" ("+str(len(extractedSequence))+")")
                    f.write("\n")
                    f.write(str(extractedSequence))
                    f.write("\n")
    f.close()
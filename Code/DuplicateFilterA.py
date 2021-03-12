############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
import os

# Methods
def getAnnotation(parts):
    if(int(parts[3])>=int(parts[4])):
        print("ERROR")
    return [parts[0],int(parts[3]),int(parts[4]),parts[6],int(parts[8])]
    # chrom, start, end, strand, transpNR
    #   0  ,   1  ,  2 ,    3  ,    4

def loadAnnotations(fileGFF):
    print(fileGFF,"\t",end="")
    counter = 0
    annotationsTransposons = {}
    annotationsFeatures    = {}
    file = open(fileGFF,"r")
    line = file.readline().replace("\n","")
    while line!="":
        if(line.startswith("#")):
            line = file.readline().replace("\n","")
            continue
        parts = line.split("\t")
        if(parts[2]=="transposon"):
            annotation = getAnnotation(parts)
            chrom = annotation[0]
            if(not chrom in annotationsTransposons):
                annotationsTransposons[chrom] = list()
            annotationsTransposons[chrom].append(annotation)
            counter += 1
        elif(parts[1]=="repeatMasker"):
            pass
        elif(parts[1]=="repeatmodel"):
            pass
        else:      
            annotation = getAnnotation(parts)
            chrom = annotation[0]
            if(not chrom in annotationsFeatures):
                annotationsFeatures[chrom] = list()
            annotationsFeatures[chrom].append(annotation)
        line = file.readline().replace("\n","")
    file.close()
    print(counter)
    return annotationsTransposons, annotationsFeatures
    
def saveAnnotations(fileGFFfrom, fileGFFto, nrList):
    file1 = open(fileGFFfrom,"r")
    file2 = open(fileGFFto, "w+")
    line = file1.readline()
    while line!="":
        if(line.startswith("#")):
            file2.write(line)
        else:
            parts = line.split("\t")
            if(not parts[1]=="ltrHarvest"):
                if(int(parts[8]) in nrList):
                    file2.write(line) 
            else:
                if(int(parts[8]) in nrList and parts[6]=="+"):
                    file2.write(line) 
        line = file1.readline()
    file1.close()
    file2.close()
    
def getTransposonNrList(annotationsTransposons):
    nrList = list()
    for chrom in list(annotationsTransposons.keys()):
        for elem in annotationsTransposons[chrom]:
            if(not elem[4] in nrList):
                nrList.append(elem[4])
    nrList.sort()
    return nrList

def removeIrrelvantFeatures(annotationsFeatures, nrList):
    annotationsFeatures2 = {}
    for chrom in list(annotationsFeatures.keys()):
        annotationsFeatures2[chrom] = list()
    for chrom in list(annotationsFeatures.keys()):
        for elem in annotationsFeatures[chrom]:
            if(elem[4] in nrList):
                annotationsFeatures2[chrom].append(elem)
    return annotationsFeatures2
    
# Start Filtering
def doFiltering(folderIn, folderOut, tolBP=0.05):
#    folderIn = "parsedAnnotations"
#    folderOut = "transposonCandA"
    for tool in ["helitronScanner", "ltrHarvest", "mitefind", "mitetracker", "must", "repeatmodel", "repMasker", "sinefind", "sinescan", "tirvish", "ltrPred"]:
        # Load Annotations
        fileIn = tool+".gff3"
        fileOut = tool+".gff3"
        if(os.path.isfile(os.path.join(folderIn,fileIn))):
            annotationsTransposons, annotationsFeatures = loadAnnotations(os.path.join(folderIn,fileIn))
            # Start Filtering
            for chrom in list(annotationsTransposons.keys()):
                iteration = 1
                while True:
                    print(tool,"\t",chrom,"\tIteration\t",iteration,"\t",len(annotationsTransposons[chrom]))
                    delIndex = list()
                    # Search for simple duplicates
                    for ix1 in range(0,len(annotationsTransposons[chrom])):
                        for ix2 in range(0,len(annotationsTransposons[chrom])):
                            if(ix1!=ix2 and not(ix2 in delIndex) and not(ix1 in delIndex)):
                                ann1 = annotationsTransposons[chrom][ix1]
                                ann2 = annotationsTransposons[chrom][ix2]
                                ann1Length = abs(ann1[2]-ann1[1])
                                if(abs(ann2[1]-ann1[1])<tolBP*ann1Length  and  abs(ann2[2]-ann1[2])<tolBP*ann1Length):
                                    delIndex.append(ix2)
                    # Delete them
                    delIndex = list(set(delIndex))
                    if(len(delIndex)==0):
                        break
                    for index in sorted(delIndex, reverse=True):
                        del annotationsTransposons[chrom][index]
                    iteration += 1
                    # chrom, start, end, strand, transpNR
                    #   0  ,   1  ,  2 ,    3  ,    4
            # Save results
            nrList = getTransposonNrList(annotationsTransposons)
            saveAnnotations(os.path.join(folderIn,fileIn), os.path.join(folderOut,fileOut), nrList)
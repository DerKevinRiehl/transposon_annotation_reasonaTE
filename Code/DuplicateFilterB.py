############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
import os
from shutil import copyfile

# Methods
def writeGFFhead(f):
    f.write("# Merged Transposon and Feature Annotations")
    f.write("\n")
    f.write("# Structure of this file: ")
    f.write("\n")
    f.write("# 1)\tChromosome")
    f.write("\n")
    f.write("# 2)\tAnnotation Software")
    f.write("\n")
    f.write("# 3)\tAnnotation Feature")
    f.write("\n")
    f.write("# 4)\tStart Position")
    f.write("\n")
    f.write("# 5)\tEnd Position")
    f.write("\n")
    f.write("# 6)\tScore / Probability")
    f.write("\n")
    f.write("# 7)\tStrand")
    f.write("\n")
    f.write("# 8)\tPhase (not used)")
    f.write("\n")
    f.write("# 9)\tTransposonNr")
    f.write("\n")
    f.write("# 10)\tAttributes separated by ;")
    f.write("\n")

def saveAnnotationsFilter(fileGFFfrom, fileGFFto, nrList):
    file1 = open(fileGFFfrom,"r")
    file2 = open(fileGFFto, "w+")
    line = file1.readline()
    while line!="":
        if(line.startswith("#")):
            file2.write(line)
        else:
            parts = line.split("\t")
            if(int(parts[8]) in nrList):
                file2.write(line) 
        line = file1.readline()
    file1.close()
    file2.close()
    
def mergeFiles(folderIn, folderOut):
    fileOut   = "CandidatesB_0.gff3"    
    counter = 0
    fileW  = open(os.path.join(folderOut,fileOut),"w+")
    writeGFFhead(fileW)
    for tool in ["helitronScanner", "ltrHarvest", "mitefind", "mitetracker", "must", "repeatmodel", "repMasker", "sinefind", "sinescan", "tirvish", "ltrPred"]:
        fileIn = tool+".gff3"
        if(os.path.isfile(os.path.join(folderIn,fileIn))):
            fileR = open(os.path.join(folderIn,fileIn),"r")
            line = fileR.readline()
            while line!="":
                if(not line.startswith("#")):
                    parts = line.replace("\n","").replace("\r","").split("\t")
                    if(parts[2]=="transposon"):
                        counter += 1
                    originID = int(parts[8])
                    parts[8] = str(counter)
                    if(len(parts)==10):
                        parts[9] = parts[9].replace("transposon "+str(originID),"transposon "+str(counter))
                    resultLine = "\t".join(parts)
                    fileW.write(resultLine)
                    fileW.write("\n")
                line = fileR.readline()
            fileR.close()
    fileW.close()



def doFilteringB(folderIn, folderOut, minTransposonLenght=25, tolBP=0.05):
#    folderIn  = "transposonCandA"
#    folderOut = "transposonCandB"

    # Merge Annotations
    mergeFiles(folderIn, folderOut)
    
    # DO NOT FILTER, seems okay no trash, MITEs are short anyway!
    # Filter small ones except for repeatmasker and repeatmodeler
    fileGFFfrom = os.path.join(folderOut,"CandidatesB_0.gff3")
    fileGFFto   = os.path.join(folderOut,"CandidatesB_1.gff3")
    nrList = list()
    file = open(fileGFFfrom,"r")
    line = file.readline()
    while line!="":
        if(not line.startswith("#")):
            parts = line.split("\t")
            if(parts[2]=="transposon" and (not parts[1]=="repeatMasker") and (not parts[1]=="repeatmodel")):
                if(abs(int(parts[4])-int(parts[3])+1) > minTransposonLenght):
                    print(parts)
                    nrList.append(int(parts[8]))
        line = file.readline()
    file.close()
    saveAnnotationsFilter(fileGFFfrom, fileGFFto, nrList)
    
    # Merge extremely similar annotations
    fileGFFfrom = os.path.join(folderOut,"CandidatesB_1.gff3")
    file = open(fileGFFfrom,"r")
    line = file.readline()
    transposons = {}
    while line!="":
        if(not line.startswith("#")):
            parts = line.split("\t")
            if(parts[2]=="transposon"):
                if(parts[0] not in transposons):
                    transposons[parts[0]] = list()
                transposons[parts[0]].append([int(parts[3]),int(parts[4]),int(parts[8])])
        line = file.readline()
    file.close()
    
    copyfile(fileGFFfrom,os.path.join(folderOut,"CandidatesB_it0.gff3"))
        
    iteration = 1
    for chrom in list(transposons.keys()):
        while True:
            fileGFFfrom = fileGFFfrom+".it0"
            couplesA = list()
            couplesB = list()
            for ix1 in range(0,len(transposons[chrom])):
                for ix2 in range(0,len(transposons[chrom])):
                    if(ix1!=ix2):
                        ann1 = transposons[chrom][ix1]
                        ann2 = transposons[chrom][ix2]
                        if(not(ann1[2] in couplesB) and not(ann1[2] in couplesA) and not(ann2[2] in couplesA) and not(ann2[2] in couplesB)):
                            ann1Length = abs(ann1[1]-ann1[0])
                            if(abs(ann1[0]-ann2[0])<tolBP*ann1Length  and  abs(ann1[1]-ann2[1])<tolBP*ann1Length):
                                couplesA.append(ann1[2])
                                couplesB.append(ann2[2])
            fileGFFfrom = os.path.join(folderOut,"CandidatesB_it"+str(iteration-1)+".gff3")
            fileGFFto   = os.path.join(folderOut,"CandidatesB_it"+str(iteration)+".gff3")
            fileR = open(fileGFFfrom,"r")
            fileW = open(fileGFFto, "w+")
            writeGFFhead(fileW)
            line = fileR.readline()
            while line!="":
                if(not line.startswith("#")):
                    parts = line.split("\t")
                    transpID = int(parts[8])
                    if(transpID in couplesB):
                        newTranspID = couplesA[couplesB.index(transpID)]
                        if(parts[2]=="transposon"):
                            pass
                        else:
                            parts[8] = str(newTranspID)
                            if(len(parts)==10):
                                parts[9] = parts[9].replace("transposon "+str(transpID),"transposon "+str(newTranspID))
                            fileW.write("\t".join(parts))
                    else:
                        fileW.write(line)
                line = fileR.readline()
            fileR.close()
            fileW.close()
            print(chrom,iteration,len(couplesA))
            print("before",len(transposons[chrom]))
            iteration += 1
            if(len(couplesA)==0):
                break
            delList = list()
            for x in range(0,len(transposons[chrom])):
                if(transposons[chrom][x][2] in couplesB):
                    delList.append(x)
            delList = list(set(delList))
            for index in sorted(delList, reverse=True):
                del transposons[chrom][index]
            print("after",len(transposons[chrom]))
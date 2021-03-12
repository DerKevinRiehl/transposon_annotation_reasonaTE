############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
from GFFTools import generateMask, mergeAnnotations, applyMask, cleanLine, cleanLine2, extractFastaFromGFF
from shutil import copyfile
import os
import os.path

# Methods
def mergeFiles(fileA,fileB,fileC):
    filenames = [fileA,fileB]
    with open(fileC, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
                
def createToolAnnotation_Files(folderProject,folderFinalResults, folderParsedAnnotations, folderTranspCandB, folderTranspCandF):
#    folderFinalResults = "finalResults"
#    folderParsedAnnotations = "parsedAnnotations"
#    folderTranspCandB = "transposonCandB"
#    folderTranspCandE = "transposonCandE"
#    folderTranspCandF = "transposonCandF"
    # Create Mask of Annotation Tools
    files = os.listdir(folderTranspCandB)
    filesN = list()
    for f in files:
        if(f.startswith("CandidatesB_it")):
            filesN.append(int(f.replace("CandidatesB_it","").replace(".gff3","")))
    filesN.sort()
    generateMask(os.path.join(folderTranspCandB,"CandidatesB_it"+str(filesN[-1])+".gff3"), os.path.join(folderFinalResults,"ToolAnnotations_TransposonMask.gff3"), "transposon")   
    # Create GFF of Transposon Annotations
    fR = open(os.path.join(folderTranspCandB,"CandidatesB_it"+str(filesN[-1])+".gff3"), "r")
    fW = open(os.path.join(folderFinalResults,"ToolAnnotations_Transposons.gff3"), "w+")
    line = fR.readline()
    while line!="":
        if(not line.startswith("#")):
            if(line.split("\t")[2]=="transposon"):
                newline = cleanLine(line)
                fW.write(newline+"\n")
        line = fR.readline()
    fR.close()
    fW.close()   
    # Create GFF of Related Structural Features
    fR = open(os.path.join(folderTranspCandB,"CandidatesB_it"+str(filesN[-1])+".gff3"), "r")
    fW = open(os.path.join(folderFinalResults,"ToolAnnotations_StructuralFeatures.gff3"), "w+")
    line = fR.readline()
    while line!="":
        if(not line.startswith("#")):
            if(not line.split("\t")[2]=="transposon"):
                fW.write(cleanLine(line)+"\n")
        line = fR.readline()
    fR.close()
    fW.close()
    # Create GFF of Related Protein Features
    mergeAnnotations([os.path.join(folderParsedAnnotations,"transposonPSI.gff3"),os.path.join(folderParsedAnnotations,"NCBICDD1000.gff3")], os.path.join(folderParsedAnnotations,"proteinfeatures.gff3"))
    applyMask(os.path.join(folderFinalResults,"ToolAnnotations_TransposonMask.gff3"), os.path.join(folderParsedAnnotations,"proteinfeatures.gff3"), os.path.join(folderParsedAnnotations,"proteinfeatures_masked.gff3"))
    fR = open(os.path.join(folderParsedAnnotations,"proteinfeatures_masked.gff3"), "r")
    fW = open(os.path.join(folderFinalResults,"ToolAnnotations_ProteinFeatures.gff3"), "w+")
    line = fR.readline()
    while line!="":
        if(not line.startswith("#")):
            if(not line.split("\t")[2]=="transposon"):
                fW.write(cleanLine2(line)+"\n")
        line = fR.readline()
    fR.close()
    fW.close()
    # Create Transposon Sequences Fasta
    extractFastaFromGFF(os.path.join(folderFinalResults,"ToolAnnotations_Transposons.gff3"), os.path.join(folderProject,"sequence.fasta"), os.path.join(folderFinalResults,"ToolAnnotations_TransposonSequences.fasta"))
    # Classify Sequences
    os.system("transposon_classifier_RFSB -mode classify -fastaFile "+os.path.join(folderFinalResults,"ToolAnnotations_TransposonSequences.fasta")+" -outputPredictionFile "+os.path.join(folderTranspCandF,"ToolAnnotations_TransposonSequencesClasses.txt"))
    # Load classification
    classesData = {}
    f = open(os.path.join(folderTranspCandF,"ToolAnnotations_TransposonSequencesClasses.txt"),"r")
    line = f.readline()
    while line!="":
        if(line.startswith(">")):
            transpNr = line.split("\t")[0].split(">Transposon")[1].rstrip()
            newline = f.readline()
            transpClassDescr = newline.split(" ")[1]+"("+newline.split(" ")[0]+")"
            classesData[transpNr] = transpClassDescr
        line = f.readline()
    f.close()
    # Create GFF of Transposon Annotations
    fR = open(os.path.join(folderTranspCandB,"CandidatesB_it"+str(filesN[-1])+".gff3"), "r")
    fW = open(os.path.join(folderFinalResults,"ToolAnnotations_Transposons.gff3"), "w+")
    line = fR.readline()
    while line!="":
        if(not line.startswith("#")):
            if(line.split("\t")[2]=="transposon"):
                newline = cleanLine(line)
                parts = newline.split("\t")
                partsD = parts[8].split(";")
                transpNr = partsD[0].split("=")[1].replace(" ","")
                partsD.insert(1,"class="+classesData[transpNr])
                parts[1] = "reasonaTE"
                parts[8] = ";".join(partsD)
                fW.write("\t".join(parts)+"\n")
        line = fR.readline()
    fR.close()
    fW.close()   

def createPipelineAnnotation_Files(folderProject, folderFinalResults, folderParsedAnnotations, folderTranspCandE, folderTranspCandF):
#    folderFinalResults = "finalResults"
#    folderParsedAnnotations = "parsedAnnotations"
#    folderTranspCandE = "transposonCandE"
#    folderTranspCandF = "transposonCandF"
    numTransposons = 0
    # get number of transposons
    f = open(os.path.join(folderFinalResults,"ToolAnnotations_Transposons.gff3"),"r")
    line = f.readline()
    while line!="":
        info = line.split("\t")[8].split(";")[0].split("=")[1].replace(" ","")
        if(int(info)>numTransposons):
            numTransposons = int(info)
        line = f.readline()
    f.close()
    # start to create final transposons
    transpNr = numTransposons+1
    fW = open(os.path.join(folderFinalResults,"PipelineAnnotations_Transposons.gff3"),"w+")
    fR = open(os.path.join(folderTranspCandE,"SequenceBlastResults_filtered6.gff3"),"r")
    line = fR.readline()
    while line!="":
        parts = line.replace("\n","").split("\t")
        fW.write(parts[0]+"\t"+"reasonaTE"+"\t"+"transposon"+"\t"+parts[3]+"\t"+parts[4]+"\t"+parts[5]+"\t"+parts[6]+"\t"+parts[7]+"\t"+"transposon="+str(transpNr)+"\n")
        transpNr += 1
        line = fR.readline()
    fR.close()    
    fW.close()
    # Extract Fasta of additional transposons
    extractFastaFromGFF(os.path.join(folderFinalResults,"PipelineAnnotations_Transposons.gff3"), os.path.join(folderProject,"sequence.fasta"), os.path.join(folderFinalResults,"PipelineAnnotations_TransposonSequences.fasta"))
    # Classify Sequences
    os.system("transposon_classifier_RFSB -mode classify -fastaFile "+os.path.join(folderFinalResults,"PipelineAnnotations_TransposonSequences.fasta")+" -outputPredictionFile "+os.path.join(folderTranspCandF,"PipelineAnnotations_TransposonSequencesClasses.txt"))
    # Load classification
    classesData = {}
    f = open(os.path.join(folderTranspCandF,"PipelineAnnotations_TransposonSequencesClasses.txt"),"r")
    line = f.readline()
    while line!="":
        if(line.startswith(">")):
            transpNr = line.split("\t")[0].split(">Transposon")[1].replace("\n","")
            newline = f.readline()
            transpClassDescr = newline.split(" ")[1]+"("+newline.split(" ")[0]+")"
            classesData[transpNr] = transpClassDescr
        line = f.readline()
    f.close()
    # start to create final transposons
    fW = open(os.path.join(folderFinalResults,"PipelineAnnotations_Transposons2.gff3"),"w+")
    fR = open(os.path.join(folderFinalResults,"PipelineAnnotations_Transposons.gff3"),"r")
    line = fR.readline()
    while line!="":
        parts = line.split("\t")
        transpNr = parts[8].split("=")[1].replace(" ","").replace("\n","")
        parts[8] = parts[8].replace("\n","")+";"+"class="+classesData[transpNr]
        parts[1] = "reasonaTE"
        fW.write("\t".join(parts)+"\n")
        line = fR.readline()
    fR.close()    
    fW.close()
    # Remove 1, rename 2
    os.remove(os.path.join(folderFinalResults,"PipelineAnnotations_Transposons.gff3"))
    os.rename(os.path.join(folderFinalResults,"PipelineAnnotations_Transposons2.gff3"),os.path.join(folderFinalResults,"PipelineAnnotations_Transposons.gff3"))
    generateMask(os.path.join(folderFinalResults,"PipelineAnnotations_Transposons.gff3"), os.path.join(folderFinalResults,"PipelineAnnotations_TransposonMask.gff3"), "transposon")   
    # Create mask and protein features
    applyMask(os.path.join(folderFinalResults,"PipelineAnnotations_TransposonMask.gff3"), os.path.join(folderParsedAnnotations,"proteinfeatures.gff3"), os.path.join(folderParsedAnnotations,"proteinfeatures_masked2.gff3"))
    fR = open(os.path.join(folderParsedAnnotations,"proteinfeatures_masked2.gff3"), "r")
    fW = open(os.path.join(folderFinalResults,"PipelineAnnotations_ProteinFeatures.gff3"), "w+")
    line = fR.readline()
    while line!="":
        if(not line.startswith("#")):
            if(not line.split("\t")[2]=="transposon"):
                fW.write(cleanLine2(line)+"\n")
        line = fR.readline()
    fR.close()
    fW.close()
    
def createFinalAnnotation_Files(folderFinalResults, folderParsedAnnotations):
#    folderFinalResults = "finalResults"
#    folderParsedAnnotations = "parsedAnnotations"
    copyfile(os.path.join(folderFinalResults,"ToolAnnotations_StructuralFeatures.gff3"), os.path.join(folderFinalResults,"FinalAnnotations_StructuralFeatures.gff3"))
    mergeFiles(os.path.join(folderFinalResults,"ToolAnnotations_TransposonSequences.fasta"),os.path.join(folderFinalResults,"PipelineAnnotations_TransposonSequences.fasta"),os.path.join(folderFinalResults,"FinalAnnotations_TransposonSequences.fasta"))
    mergeFiles(os.path.join(folderFinalResults,"ToolAnnotations_Transposons.gff3"),os.path.join(folderFinalResults,"PipelineAnnotations_Transposons.gff3"),os.path.join(folderFinalResults,"FinalAnnotations_Transposons.gff3"))
    generateMask(os.path.join(folderFinalResults,"FinalAnnotations_Transposons.gff3"), os.path.join(folderFinalResults,"FinalAnnotations_TransposonMask.gff3"), "transposon")   
    applyMask(os.path.join(folderFinalResults,"FinalAnnotations_TransposonMask.gff3"), os.path.join(folderParsedAnnotations,"proteinfeatures.gff3"), os.path.join(folderParsedAnnotations,"proteinfeatures_masked3.gff3"))
    fR = open(os.path.join(folderParsedAnnotations,"proteinfeatures_masked3.gff3"), "r")
    fW = open(os.path.join(folderFinalResults,"FinalAnnotations_ProteinFeatures.gff3"), "w+")
    line = fR.readline()
    while line!="":
        if(not line.startswith("#")):
            if(not line.split("\t")[2]=="transposon"):
                fW.write(cleanLine2(line)+"\n")
        line = fR.readline()
    fR.close()
    fW.close()

















#def createToolAnnotation_Files():
#    # Create Mask of Annotation Tools
#    files = os.listdir("transposonCandB")
#    filesN = list()
#    for f in files:
#        if(f.startswith("CandidatesB_it")):
#            filesN.append(int(f.replace("CandidatesB_it","").replace(".gff3","")))
#    filesN.sort()
#    generateMask(os.path.join("transposonCandB","CandidatesB_it"+str(filesN[-1])+".gff3"), "finalResults/ToolAnnotations_TransposonMask.gff3", "transposon")   
#    # Create GFF of Transposon Annotations
#    fR = open(os.path.join("transposonCandB","CandidatesB_it"+str(filesN[-1])+".gff3"), "r")
#    fW = open("finalResults/ToolAnnotations_Transposons.gff3", "w+")
#    line = fR.readline()
#    while line!="":
#        if(not line.startswith("#")):
#            if(line.split("\t")[2]=="transposon"):
#                newline = cleanLine(line)
#                fW.write(newline+"\n")
#        line = fR.readline()
#    fR.close()
#    fW.close()   
#    # Create GFF of Related Structural Features
#    fR = open(os.path.join("transposonCandB","CandidatesB_it"+str(filesN[-1])+".gff3"), "r")
#    fW = open("finalResults/ToolAnnotations_StructuralFeatures.gff3", "w+")
#    line = fR.readline()
#    while line!="":
#        if(not line.startswith("#")):
#            if(not line.split("\t")[2]=="transposon"):
#                fW.write(cleanLine(line)+"\n")
#        line = fR.readline()
#    fR.close()
#    fW.close()
#    # Create GFF of Related Protein Features
#    mergeAnnotations(["parsedAnnotations/transposonPSI.gff3","parsedAnnotations/NCBICDD1000.gff3"], "parsedAnnotations/proteinfeatures.gff3")
#    applyMask("finalResults/ToolAnnotations_TransposonMask.gff3", "parsedAnnotations/proteinfeatures.gff3", "parsedAnnotations/proteinfeatures_masked.gff3")
#    fR = open("parsedAnnotations/proteinfeatures_masked.gff3", "r")
#    fW = open("finalResults/ToolAnnotations_ProteinFeatures.gff3", "w+")
#    line = fR.readline()
#    while line!="":
#        if(not line.startswith("#")):
#            if(not line.split("\t")[2]=="transposon"):
#                fW.write(cleanLine2(line)+"\n")
#        line = fR.readline()
#    fR.close()
#    fW.close()
#    # Create Transposon Sequences Fasta
#    extractFastaFromGFF("finalResults/ToolAnnotations_Transposons.gff3", "parsedAnnotations/sequence.fasta", "finalResults/ToolAnnotations_TransposonSequences.fasta")
#    # Classify Sequences
#    os.system("transposon_classifier_RFSB -mode classify -fastaFile finalResults/ToolAnnotations_TransposonSequences.fasta -outputPredictionFile transposonCandF/ToolAnnotations_TransposonSequencesClasses.txt")
#    # Load classification
#    classesData = {}
#    f = open("transposonCandF/ToolAnnotations_TransposonSequencesClasses.txt","r")
#    line = f.readline()
#    while line!="":
#        if(line.startswith(">")):
#            transpNr = line.split("\t")[0].split(">Transposon")[1]
#            newline = f.readline()
#            transpClassDescr = newline.split(" ")[1]+"("+newline.split(" ")[0]+")"
#            classesData[transpNr] = transpClassDescr
#        line = f.readline()
#    f.close()
#    # Create GFF of Transposon Annotations
#    fR = open(os.path.join("transposonCandB","CandidatesB_it"+str(filesN[-1])+".gff3"), "r")
#    fW = open("finalResults/ToolAnnotations_Transposons.gff3", "w+")
#    line = fR.readline()
#    while line!="":
#        if(not line.startswith("#")):
#            if(line.split("\t")[2]=="transposon"):
#                newline = cleanLine(line)
#                parts = newline.split("\t")
#                partsD = parts[8].split(";")
#                transpNr = partsD[0].split("=")[1].replace(" ","")
#                partsD.insert(1,"class="+classesData[transpNr])
#                parts[1] = "reasonaTE"
#                parts[8] = ";".join(partsD)
#                fW.write("\t".join(parts)+"\n")
#        line = fR.readline()
#    fR.close()
#    fW.close()   
#
#def createPipelineAnnotation_Files():
#    numTransposons = 0
#    # get number of transposons
#    f = open("finalResults/ToolAnnotations_Transposons.gff3","r")
#    line = f.readline()
#    while line!="":
#        info = line.split("\t")[8].split(";")[0].split("=")[1].replace(" ","")
#        if(int(info)>numTransposons):
#            numTransposons = int(info)
#        line = f.readline()
#    f.close()
#    # start to create final transposons
#    transpNr = numTransposons+1
#    fW = open("finalResults/PipelineAnnotations_Transposons.gff3","w+")
#    fR = open("transposonCandE/SequenceBlastResults_filtered6.gff3","r")
#    line = fR.readline()
#    while line!="":
#        parts = line.replace("\n","").split("\t")
#        fW.write(parts[0]+"\t"+"reasonaTE"+"\t"+"transposon"+"\t"+parts[3]+"\t"+parts[4]+"\t"+parts[5]+"\t"+parts[6]+"\t"+parts[7]+"\t"+"transposon="+str(transpNr)+"\n")
#        transpNr += 1
#        line = fR.readline()
#    fR.close()    
#    fW.close()
#    # Extract Fasta of additional transposons
#    extractFastaFromGFF("finalResults/PipelineAnnotations_Transposons.gff3", "parsedAnnotations/sequence.fasta", "finalResults/PipelineAnnotations_TransposonSequences.fasta")
#    # Classify Sequences
#    os.system("transposon_classifier_RFSB -mode classify -fastaFile finalResults/PipelineAnnotations_TransposonSequences.fasta -outputPredictionFile transposonCandF/PipelineAnnotations_TransposonSequencesClasses.txt")
#    # Load classification
#    classesData = {}
#    f = open("transposonCandF/PipelineAnnotations_TransposonSequencesClasses.txt","r")
#    line = f.readline()
#    while line!="":
#        if(line.startswith(">")):
#            transpNr = line.split("\t")[0].split(">Transposon")[1].replace("\n","")
#            newline = f.readline()
#            transpClassDescr = newline.split(" ")[1]+"("+newline.split(" ")[0]+")"
#            classesData[transpNr] = transpClassDescr
#        line = f.readline()
#    f.close()
#    # start to create final transposons
#    fW = open("finalResults/PipelineAnnotations_Transposons2.gff3","w+")
#    fR = open("finalResults/PipelineAnnotations_Transposons.gff3","r")
#    line = fR.readline()
#    while line!="":
#        parts = line.split("\t")
#        transpNr = parts[8].split("=")[1].replace(" ","").replace("\n","")
#        parts[8] = parts[8].replace("\n","")+";"+"class="+classesData[transpNr]
#        parts[1] = "reasonaTE"
#        fW.write("\t".join(parts)+"\n")
#        line = fR.readline()
#    fR.close()    
#    fW.close()
#    # Remove 1, rename 2
#    os.remove("finalResults/PipelineAnnotations_Transposons.gff3") 
#    os.rename("finalResults/PipelineAnnotations_Transposons2.gff3","finalResults/PipelineAnnotations_Transposons.gff3")
#    generateMask("finalResults/PipelineAnnotations_Transposons.gff3", "finalResults/PipelineAnnotations_TransposonMask.gff3", "transposon")   
#    # Create mask and protein features
#    applyMask("finalResults/PipelineAnnotations_TransposonMask.gff3", "parsedAnnotations/proteinfeatures.gff3", "parsedAnnotations/proteinfeatures_masked2.gff3")
#    fR = open("parsedAnnotations/proteinfeatures_masked2.gff3", "r")
#    fW = open("finalResults/PipelineAnnotations_ProteinFeatures.gff3", "w+")
#    line = fR.readline()
#    while line!="":
#        if(not line.startswith("#")):
#            if(not line.split("\t")[2]=="transposon"):
#                fW.write(cleanLine2(line)+"\n")
#        line = fR.readline()
#    fR.close()
#    fW.close()
#    
#def createFinalAnnotation_Files():
#    copyfile("finalResults/ToolAnnotations_StructuralFeatures.gff3", "finalResults/FinalAnnotations_StructuralFeatures.gff3")
#    mergeFiles("finalResults/ToolAnnotations_TransposonSequences.fasta","finalResults/PipelineAnnotations_TransposonSequences.fasta","finalResults/FinalAnnotations_TransposonSequences.fasta")
#    mergeFiles("finalResults/ToolAnnotations_Transposons.gff3","finalResults/PipelineAnnotations_Transposons.gff3","finalResults/FinalAnnotations_Transposons.gff3")
#    generateMask("finalResults/FinalAnnotations_Transposons.gff3", "finalResults/FinalAnnotations_TransposonMask.gff3", "transposon")   
#    applyMask("finalResults/FinalAnnotations_TransposonMask.gff3", "parsedAnnotations/proteinfeatures.gff3", "parsedAnnotations/proteinfeatures_masked3.gff3")
#    fR = open("parsedAnnotations/proteinfeatures_masked3.gff3", "r")
#    fW = open("finalResults/FinalAnnotations_ProteinFeatures.gff3", "w+")
#    line = fR.readline()
#    while line!="":
#        if(not line.startswith("#")):
#            if(not line.split("\t")[2]=="transposon"):
#                fW.write(cleanLine2(line)+"\n")
#        line = fR.readline()
#    fR.close()
#    fW.close()

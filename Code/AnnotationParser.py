############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import os
from AnnotationChecker import checkHelitronScanner, checkLtrHarvest, checkLtrPred, checkMiteFind, checkMiteTracker, checkMust, checkNCBICDD1000, checkRepeatModeler, checkRepeatMasker, checkSineFind, checkSineScan, checkTirVish, checkTransposonPSI 

# Methods
def getArgument(args, title):
    for i in range(0, len(args)):
        if(args[i].startswith("-"+title)):
            if(i<len(args)-1):
                return args[i+1]
            else:
                return ""
    return ""

def correctPositions(start, end):
    start = int(start)
    end = int(end)
    if(start>end):
        c = end
        end = start
        start = c
    start = str(start)
    end = str(end)
    return start, end

def writeGFFhead(f):
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
    
def parseTirvishLine(line):
    parts = line.split("\t")
    chromosome = parts[0].replace(">","")
    annotationSoftware = "TIRvish"
    feature = ""
    if(parts[2].startswith("repeat_regio")):
        feature = "transposon"
    elif(parts[2].startswith("target_site")):
        feature = "tsd"
    elif(parts[2].startswith("terminal_inverted_repeat_element")):
        feature = "ignore"
    elif(parts[2].startswith("terminal_inverted_repeat")):
        feature = "tir"
    start = parts[3]
    end = parts[4]
    start,end    =  correctPositions(start, end)
    score = "."
    strand = "+"
    phase = "."
    return chromosome, annotationSoftware, feature, start, end, score, strand, phase

def parseTirvishLineRC(line, lengthDict):
    parts = line.split("\t")
    chromosome = parts[0].replace(">","")
    annotationSoftware = "TIRvish"
    feature = ""
    if(parts[2].startswith("repeat_regio")):
        feature = "transposon"
    elif(parts[2].startswith("target_site")):
        feature = "tsd"
    elif(parts[2].startswith("terminal_inverted_repeat_element")):
        feature = "ignore"
    elif(parts[2].startswith("terminal_inverted_repeat")):
        feature = "tir"
    start = int(parts[3])
    end = int(parts[4])
    start = lengthDict[chromosome] - int(start)
    end = lengthDict[chromosome] - int(end)  
    start,end    =  correctPositions(start, end)      
    score = "."
    strand = "-"
    phase = "."
    return chromosome, annotationSoftware, feature, start, end, score, strand, phase

def exportFastaFile(targetFastaFile, sourceFastaFile, chromosomes, transposonAnnotations):
    f = open(targetFastaFile,"w+")
    records = SeqIO.parse(sourceFastaFile, "fasta")
    for r in records:   
        for key in chromosomes:
            chromosome = r.id
            if(chromosome==key):
                sequences  = r.seq
                for i in range(0,len(transposonAnnotations[key])):
                    parts  = transposonAnnotations[key][i].split("\t")
                    typ = parts[2]
                    if(typ=="transposon"):
                        start  = parts[3]
                        end    = parts[4]
                        strand = parts[6]
                        number = parts[8]
                        strandNum = 1
                        if(strand=="-"):
                            strandNum = -1
                        transposonFeature = SeqFeature(FeatureLocation(int(start), int(end), strand=strandNum), type="CDS")
                        extractedSequence  = transposonFeature.extract(sequences)
                        f.write(">Transposon "+str(number)+"\t"+chromosome+":"+"["+str(strand)+"]"+str(start)+"-"+str(end)+" ("+str(len(extractedSequence))+")")
                        f.write("\n")
                        f.write(str(extractedSequence))
                        f.write("\n")
    f.close()
    
def loadRepBaseTaxonomyScheme():
    script_dir = os.path.dirname(__file__)
    f = open(os.path.join(script_dir,"ncbi_cdd_candidates","TaxonomyScheme.csv"),"r")
    linesA = f.readlines()
    f.close()
    f = open(os.path.join(script_dir,"ncbi_cdd_candidates","TaxonomyScheme2.csv"),"r")
    linesB = f.readlines()
    f.close()
    f = open(os.path.join(script_dir,"ncbi_cdd_candidates","TaxonomyScheme3.csv"),"r")
    linesC = f.readlines()
    f.close()
    for i in range(0,len(linesA)):
        linesA[i] = linesA[i].replace("\n","")
    for i in range(0,len(linesB)):
        linesB[i] = linesB[i].replace("\n","")
    for i in range(0,len(linesC)):
        linesC[i] = linesC[i].replace("\n","")
    scheme = {}
    for l in linesA:
        parts = l.split("\t")
        if(parts[0]=="RepBase23.08"):
            key = parts[1]
            val = parts[2]
            scheme[key] = val
    for l in linesB:
        parts = l.split("\t")
        for key in list(scheme.keys()):
            if(scheme[key]==parts[0]):
                scheme[key] = parts[1]
    for l in linesC:
        parts = l.split("\t")
        if(not parts[0] in scheme):
            scheme[parts[0]] = parts[1]
    return scheme

def getTransposonClass(label, scheme):
    for key in list(scheme.keys()):
        if(label==key):
            return scheme[key]
    if(label.startswith("LINE")):
        return "1.2.1"
    if(label.startswith("SINE")):
        return "1.2.2"
    if(label.startswith("SINE")):
        return "1.2.2"
    if(label.startswith("DNA")):
        return "2"
    if(label.startswith("LTR/ERV")):
        return "1.3.1"
    if(label.startswith("LTR")):
        return "1.1"
    else:
        return ""
    
def getFeatureType(label, scheme):
    if(not getTransposonClass(label,scheme)==""):
        return "transposon"
    else:
        if(label=="Satellite"):
            return "satellite"
        elif(label=="Simple_repeat" or label=="Low_complexity"):
            return "simple_repeat"
        else:
            return ""
    
def loadCandidates():
    script_dir = os.path.dirname(__file__)
    files = os.listdir(os.path.join(script_dir,"ncbi_cdd_candidates"))
    candidates = {}
    for file in files:
        category = file.split(".")[0].replace("cdd_result_","").split("_")[0]
        f = open(os.path.join(script_dir,"ncbi_cdd_candidates", file),"r")
        lines = f.readlines()
        f.close()
        for l in lines:
            if(l.startswith("Accession")):
                code = l.split(":")[1].replace("ID","").replace(" ","")
                if(not code in candidates):
                    candidates[code] = ""
                if(not category in candidates[code]):
                    candidates[code] += category+","
    return candidates

def getRepeatModelTransposonType(description):
    parts = description.split(";")
    if(parts[2].startswith("DNA_Transposon")):
        if("Helitr" in description):
            return "2.2"
        elif(parts[3].startswith("Terminal_Inverted_Rep")):
            if("Tc1-Mariner-like" in description):
                return "2.1.1"
            else:
                return "2.1"
    elif(parts[2].startswith("Retrotransposed_Element")):
        if("LINE" in description):
            return "1.2.1"
        elif("SINE" in description):
            return "1.2.2"
        elif("Long_Terminal_Repeat" in description):
            if("Gypsy" in description):
                return "1.1.2"
            elif("Copy" in description):
                return "1.1.1"
            else:
                return "1.1"
    else:
        return "?"
    return "?"

def parseHelitronScanner(pathResDir, pathResDirRC, fastaFile, targetGFFFile, targetFastaFile):
    print("Parse HelitronScanner outputs...")
    # Get Records lengths
    lengthDict = {}
    for record in SeqIO.parse(fastaFile, "fasta"):
        print("%s %i" % (record.id, len(record)))
        lengthDict[record.id] = len(record)
    # Read Transposon Annotations from file
    transposonAnnotations = {}
    f = open(os.path.join(pathResDir,"result.txt"),"r")
    line = f.readline()
    counter = 0
    while line!="":
        if(line.startswith(">")):
            chrom = line.replace(">","").replace("\n","")
            line = f.readline().replace("\n","")
            if(line==""):
                break
            helitrons = line.split("]")
            if(not chrom in transposonAnnotations):
                transposonAnnotations[chrom] = list()
            annoSoftware = "helitronScanner"
            features     = "transposon"
            for i in range(0,len(helitrons)):
                counter += 1
                if(helitrons[i].startswith("[")):
                    continue
                if(helitrons[i]==" "):
                    break
                start        = helitrons[i].split(":")[0].rstrip().lstrip()
                end          = helitrons[i].split(":")[1].split(" ")[0].rstrip().lstrip()
                start,end    =  correctPositions(start, end)
                score        = helitrons[i].split("[")[1].rstrip().lstrip()
                strand       = "+"
                phase        = "."
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+start+"\t"+end+"\t"+score+"\t"+strand+"\t"+phase+"\t")
        line = f.readline()
    f.close()
    f = open(os.path.join(pathResDirRC,"result.txt"),"r")
    line = f.readline()
    while line!="":
        if(line.startswith(">")):
            chrom = line.replace(">","").replace("\n","")
            line  = f.readline().replace("\n","")
            if(line==""):
                break
            helitrons = line.split("]")
            if(not chrom in transposonAnnotations):
                transposonAnnotations[chrom] = list()
            annoSoftware = "helitronScanner"
            features     = "transposon"
            for i in range(0,len(helitrons)):
                counter += 1
                if(helitrons[i].startswith("[")):
                    continue
                if(helitrons[i]==" "):
                    break
                start        = helitrons[i].split(":")[0].rstrip().lstrip()
                start        = lengthDict[chrom] - int(start)
                end          = helitrons[i].split(":")[1].split(" ")[0].rstrip().lstrip()
                end          = lengthDict[chrom] - int(end)
                start,end    =  correctPositions(start, end)
                score        = helitrons[i].split("[")[1].rstrip().lstrip()
                strand       = "-"
                phase        = "."
                start = str(start)
                end   = str(end)
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+start+"\t"+end+"\t"+score+"\t"+strand+"\t"+phase+"\t")
        line = f.readline()
    f.close()
    # Sort Annotations by Chromosome and Start Coordinate
    for key in transposonAnnotations:
        startCoords = [int(x.split("\t")[3]) for x in transposonAnnotations[key]]
        transposons = transposonAnnotations[key]
        transposons = [x for _,x in sorted(zip(startCoords,transposons))]
        transposonAnnotations[key] = transposons
    # Print results to Annotation file
    f = open(targetGFFFile,"w+")
    f.write("# HelitronScanner Annotation")
    writeGFFhead(f)
    keys = list(transposonAnnotations.keys())
    keys.sort()
    counter = 0
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            counter += 1
            transposonAnnotations[key][i] = transposonAnnotations[key][i]+str(counter)
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()
    # Get Transposon Sequences
    exportFastaFile(targetFastaFile, fastaFile, keys, transposonAnnotations)
    
def getFile(folder, suffix):
    files = os.listdir(folder)
    for f in files:
        if(f.endswith(suffix)):
            return f
    return ""

def parseLtrHarvest(pathResDir, fastaFile, targetGFFFile, targetFastaFile):
    print("Parse ltrHarvest outputs...")
    # Read Order of Chromosomes in Fasta file
    chromosomeList = list()
    f = open(fastaFile,"r")
    line = f.readline()
    while line!="":
        if(line.startswith(">")):
            if(not line.replace(">","").replace("\n","") in  chromosomeList):
                chromosomeList.append(line.replace(">","").replace("\n",""))
        line = f.readline()
    f.close()
    # Read Transposon Annotations from files
    transposonAnnotations = {}
    f = open(os.path.join(pathResDir,getFile(pathResDir,".txt")),"r")
    line = f.readline()
    counter = 0
    while line!="":
        if(not line.startswith("#")):
            transposons = line.replace("\n","").replace("  "," ").split(" ")
            chrom = chromosomeList[int(transposons[10])].replace(">","")
            if(not chrom in transposonAnnotations):
                transposonAnnotations[chrom] = list()
            counter += 1
            annoSoftware = "ltrHarvest"
            features     = "transposon"
            start  = int(transposons[0])
            end    = int(transposons[1])
            start,end    =  correctPositions(start, end)
            score  = float(transposons[9])
            strand = "+"
            phase  = "."
            transpNr = str(counter)
            attributes = ""
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+"-"+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            features = "ltr"
            start  = int(transposons[3])
            end    = int(transposons[4])
            start,end    =  correctPositions(start, end)
            score  = "."
            strand = "+"
            phase  = "."
            transpNr = str(counter)
            attributes = "Left LTR of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+"-"+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            features = "ltr"
            start  = int(transposons[6])
            end    = int(transposons[7])
            start,end    =  correctPositions(start, end)
            score  = "."
            strand = "+"
            phase  = "."
            transpNr = str(counter)
            attributes = "Right LTR of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+"-"+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        line = f.readline()
    f.close()
    # Print results to Annotation file
    f = open(targetGFFFile,"w+")
    f.write("# ltrHarvest Annotation")
    writeGFFhead(f)
    keys = list(transposonAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()
    # Get Transposon Sequences
    exportFastaFile(targetFastaFile, fastaFile, keys, transposonAnnotations)

def parseTirvish(pathResDir, pathResDirRC, fastaFile, targetGFFFile, targetFastaFile):
    print("Parse tirvish...")
    # Get Records lengths
    lengthDict = {}
    for record in SeqIO.parse(fastaFile, "fasta"):
        print("%s %i" % (record.id, len(record)))
        lengthDict[record.id] = len(record)
    # Parse Annotations
    transposonAnnotations = {}
    counter = 0
    f = open(os.path.join(pathResDir,getFile(pathResDir,"result.txt")),"r")
    line = f.readline()
    while line!="":
        if(not line.startswith("#")):
            lLines = list()
            lLines.append(line.replace("\n",""))
            line = f.readline()
            while line!="" and not line.startswith("#"):
                lLines.append(line.replace("\n",""))
                line = f.readline()
            counter += 1
            chromosome, software, feature, start, end, score, strand, phase = parseTirvishLine(lLines[0])
            if(not chromosome in transposonAnnotations):
                transposonAnnotations[chromosome] = list()
            chrom, software, feature, start, end, score, strand, phase = parseTirvishLine(lLines[0])
            transposonAnnotations[chromosome].append(chrom+"\t"+software+"\t"+feature+"\t"+start+"\t"+end+"\t"+score+"\t"+strand+"\t"+phase+"\t"+str(counter))
            chrom, software, feature, start, end, score, strand, phase = parseTirvishLine(lLines[1])
            transposonAnnotations[chromosome].append(chrom+"\t"+software+"\t"+feature+"\t"+start+"\t"+end+"\t"+score+"\t"+strand+"\t"+phase+"\t"+str(counter)+"\tLeft TSD of transposon "+str(counter))
            chrom, software, feature, start, end, score, strand, phase = parseTirvishLine(lLines[3])
            transposonAnnotations[chromosome].append(chrom+"\t"+software+"\t"+feature+"\t"+start+"\t"+end+"\t"+score+"\t"+strand+"\t"+phase+"\t"+str(counter)+"\tLeft TIR of transposon "+str(counter))
            chrom, software, feature, start, end, score, strand, phase = parseTirvishLine(lLines[4])
            transposonAnnotations[chromosome].append(chrom+"\t"+software+"\t"+feature+"\t"+start+"\t"+end+"\t"+score+"\t"+strand+"\t"+phase+"\t"+str(counter)+"\tRight TIR of transposon "+str(counter))
            chrom, software, feature, start, end, score, strand, phase = parseTirvishLine(lLines[5])
            transposonAnnotations[chromosome].append(chrom+"\t"+software+"\t"+feature+"\t"+start+"\t"+end+"\t"+score+"\t"+strand+"\t"+phase+"\t"+str(counter)+"\tRight TSD of transposon "+str(counter))
        line = f.readline()
    f.close()
    f = open(os.path.join(pathResDirRC,getFile(pathResDirRC,"result.txt")),"r")
    line = f.readline()
    while line!="":
        if(not line.startswith("#")):
            lLines = list()
            lLines.append(line.replace("\n",""))
            line = f.readline()
            while line!="" and not line.startswith("#"):
                lLines.append(line.replace("\n",""))
                line = f.readline()
            counter += 1
            chromosome, software, feature, start, end, score, strand, phase = parseTirvishLineRC(lLines[0], lengthDict)
            if(not chromosome in transposonAnnotations):
                transposonAnnotations[chromosome] = list()
            chrom, software, feature, start, end, score, strand, phase = parseTirvishLineRC(lLines[0], lengthDict)
            transposonAnnotations[chromosome].append(chrom+"\t"+software+"\t"+feature+"\t"+start+"\t"+end+"\t"+score+"\t"+strand+"\t"+phase+"\t"+str(counter))
            chrom, software, feature, start, end, score, strand, phase = parseTirvishLineRC(lLines[1], lengthDict)
            transposonAnnotations[chromosome].append(chrom+"\t"+software+"\t"+feature+"\t"+start+"\t"+end+"\t"+score+"\t"+strand+"\t"+phase+"\t"+str(counter)+"\tLeft TSD of transposon "+str(counter))
            chrom, software, feature, start, end, score, strand, phase = parseTirvishLineRC(lLines[3], lengthDict)
            transposonAnnotations[chromosome].append(chrom+"\t"+software+"\t"+feature+"\t"+start+"\t"+end+"\t"+score+"\t"+strand+"\t"+phase+"\t"+str(counter)+"\tLeft TIR of transposon "+str(counter))
            chrom, software, feature, start, end, score, strand, phase = parseTirvishLineRC(lLines[4], lengthDict)
            transposonAnnotations[chromosome].append(chrom+"\t"+software+"\t"+feature+"\t"+start+"\t"+end+"\t"+score+"\t"+strand+"\t"+phase+"\t"+str(counter)+"\tRight TIR of transposon "+str(counter))
            chrom, software, feature, start, end, score, strand, phase = parseTirvishLineRC(lLines[5], lengthDict)
            transposonAnnotations[chromosome].append(chrom+"\t"+software+"\t"+feature+"\t"+start+"\t"+end+"\t"+score+"\t"+strand+"\t"+phase+"\t"+str(counter)+"\tRight TSD of transposon "+str(counter))
        line = f.readline()
    f.close()
    # Print results to Annotation file
    f = open(targetGFFFile,"w+")
    f.write("# tirvish Annotation")
    writeGFFhead(f)
    keys = list(transposonAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()
    # Get Transposon Sequences
    exportFastaFile(targetFastaFile, fastaFile, keys, transposonAnnotations)

def parseMust(pathResDir, fastaFile, targetGFFFile, targetFastaFile):
    print("Parse must...")
    transposonAnnotations = {}
    f = open(os.path.join(pathResDir,"result.txt"),"r")
    line = f.readline()
    counter = 0
    while line!="":
        if(not line.startswith("#")):
            transposons = line.replace("\n","").split("\t")
            chrom = transposons[0].replace(">","")
            if(not chrom in transposonAnnotations):
                transposonAnnotations[chrom] = list()
            counter += 1
            annoSoftware = "must"
            features     = "transposon"
            strand = transposons[5]
            start  = int(transposons[3])
            end    = int(transposons[4])
            start,end    =  correctPositions(start, end)
            score  = float(transposons[22])
            phase  = "."
            transpNr = str(counter)
            if(transposons[23]=="Full Copy"):
                attributes = "Complete transposon"
            else:
                attributes = "Incomplete transposon"
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            if(float(transposons[7])!=0):
                features     = "dr"
                start  = int(transposons[3])-int(transposons[7])
                end    = int(transposons[3])
                start,end    =  correctPositions(start, end)
                score  = float(transposons[8])
                strand = transposons[5]
                phase  = "."
                transpNr = str(counter)
                attributes = "Left direct repeat of transposon "+str(transpNr)
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
                features     = "dr"
                start  = int(transposons[4])
                end    = int(transposons[4])+int(transposons[7])
                start,end    =  correctPositions(start, end)
                score  = float(transposons[8])
                strand = transposons[5]
                phase  = "."
                transpNr = str(counter)
                attributes = "Right direct repeat of transposon "+str(transpNr)
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            if(float(transposons[12])!=0):
                features     = "tir"
                start  = int(transposons[3])
                end    = int(transposons[3])+int(transposons[7])
                start,end    =  correctPositions(start, end)
                score  = float(transposons[8])
                strand = transposons[5]
                phase  = "."
                transpNr = str(counter)
                attributes = "Left TIR of transposon "+str(transpNr)
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
                features     = "tir"
                start  = int(transposons[4])-int(transposons[7])
                end    = int(transposons[4])
                start,end    =  correctPositions(start, end)
                score  = float(transposons[8])
                strand = transposons[5]
                phase  = "."
                transpNr = str(counter)
                attributes = "Right TIR of transposon "+str(transpNr)
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        line = f.readline()
    f.close()
    # Print results to Annotation file
    f = open(targetGFFFile,"w+")
    f.write("# must Annotation")
    writeGFFhead(f)
    keys = list(transposonAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()
    # Get Transposon Sequences
    exportFastaFile(targetFastaFile, fastaFile, keys, transposonAnnotations)

def parseMiteFind(pathResDir, pathResDirRC, fastaFile, targetGFFFile, targetFastaFile):
    print("Parse miteFinderII outputs...")
    # Get Records lengths
    lengthDict = {}
    for record in SeqIO.parse(fastaFile, "fasta"):
        print("%s %i" % (record.id, len(record)))
        lengthDict[record.id] = len(record)
    # Read Order of Chromosomes in Fasta file
    chromosomeList = list()
    f = open(fastaFile,"r")
    line = f.readline()
    while line!="":
        if(line.startswith(">")):
            if(not line.replace("\n","").replace(">","") in chromosomeList):
                chromosomeList.append(line.replace("\n","").replace(">",""))
        line = f.readline()
    f.close()
    transposonAnnotations = {}
    f = open(os.path.join(pathResDir,"result.txt"),"r")
    line = f.readline()
    counter = 0
    while line!="":
        if(line.startswith(">")):
            transposons = line.replace("\n","").split("|")
            chrom = chromosomeList[int(transposons[1])-1]
            if(not chrom in transposonAnnotations):
                transposonAnnotations[chrom] = list()
            counter += 1
            annoSoftware = "miteFinderII"
            features     = "transposon"
            strand = "+"
            start  = int(transposons[2])
            end    = int(transposons[5])
            start,end    =  correctPositions(start, end)
            score  = float(transposons[9].split(":")[1])
            phase  = "."
            transpNr = str(counter)
            attributes = ""
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            features     = "tir"
            start  = int(transposons[2])
            end    = int(transposons[3])
            start,end    =  correctPositions(start, end)
            transpNr = str(counter)
            attributes = "Left TIR of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            features     = "tir"
            start  = int(transposons[4])
            end    = int(transposons[5])
            start,end    =  correctPositions(start, end)
            transpNr = str(counter)
            attributes = "Right TIR of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        line = f.readline()
    f.close()
    f = open(os.path.join(pathResDirRC,"result.txt"),"r")
    line = f.readline()
    while line!="":
        if(line.startswith(">")):
            transposons = line.replace("\n","").split("|")
            chrom = chromosomeList[int(transposons[1])-1]
            if(not chrom in transposonAnnotations):
                transposonAnnotations[chrom] = list()
            counter += 1
            annoSoftware = "miteFinderII"
            features     = "transposon"
            strand = "-"
            start  = int(transposons[2])
            start  = lengthDict[chrom] - int(start)
            end    = int(transposons[5])
            end    = lengthDict[chrom] - int(end)
            start,end    =  correctPositions(start, end)
            score  = float(transposons[9].split(":")[1])
            phase  = "."
            transpNr = str(counter)
            attributes = ""
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            features     = "tir"
            start  = int(transposons[2])
            start  = lengthDict[chrom] - int(start)
            end    = int(transposons[3])
            end    = lengthDict[chrom] - int(end)
            start,end    =  correctPositions(start, end)
            transpNr = str(counter)
            attributes = "Left TIR of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            features     = "tir"
            start  = int(transposons[4])
            start  = lengthDict[chrom] - int(start)
            end    = int(transposons[5])
            end    = lengthDict[chrom] - int(end)
            start,end    =  correctPositions(start, end)
            transpNr = str(counter)
            attributes = "Right TIR of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        line = f.readline()
    f.close()
    # Print results to Annotation file
    f = open(targetGFFFile,"w+")
    f.write("# miteFinderII Annotation ")
    writeGFFhead(f)
    keys = list(transposonAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()
    # Get Transposon Sequences
    exportFastaFile(targetFastaFile, fastaFile, keys, transposonAnnotations)

def parseSinefind(pathResDir, pathResDirRC, fastaFile, targetGFFFile, targetFastaFile):
    print("Parse sineFinder outputs...")
    # Get Records lengths
    lengthDict = {}
    for record in SeqIO.parse(fastaFile, "fasta"):
        print("%s %i" % (record.id, len(record)))
        lengthDict[record.id] = len(record)
    # Load Annotations
    transposonAnnotations = {}
    counter = 0
    f = open(os.path.join(pathResDir,getFile(pathResDir,"-matches.fasta")),"r")
    line = f.readline()
    while line!="":
        if(line.startswith(">")):
            transposons = line.replace("\n","").split(" ")
            chrom = transposons[0].replace(">","")
            if(not chrom in transposonAnnotations):
                transposonAnnotations[chrom] = list()
            counter += 1
            annoSoftware = "sineFinder"
            features     = "transposon"
            strand = "+"
            start  = int(transposons[2].split(":")[0])
            end    = int(transposons[2].split(":")[1])
            start,end    =  correctPositions(start, end)
            score  = "."
            phase  = "."
            transpNr = str(counter)
            attributes = ""
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            features     = "tsd"
            start  = int(transposons[2].split(":")[0])
            end    = int(transposons[2].split(":")[0])+int(transposons[3].split(";")[0].split("=")[1])
            start,end    =  correctPositions(start, end)
            score  = int(transposons[3].split(";")[1].split("=")[1])
            phase  = "."
            transpNr = str(counter)
            attributes = "Left TSD of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            features     = "tsd"
            start  = int(transposons[2].split(":")[1])-int(transposons[3].split(";")[0].split("=")[1])
            end    = int(transposons[2].split(":")[1])
            start,end = correctPositions(start, end)
            score  = int(transposons[3].split(";")[1].split("=")[1])
            phase  = "."
            transpNr = str(counter)
            attributes = "Right TSD of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        line = f.readline()
    f.close()
    f = open(os.path.join(pathResDirRC,getFile(pathResDirRC,"-matches.fasta")),"r")
    line = f.readline()
    while line!="":
        if(line.startswith(">")):
            transposons = line.replace("\n","").split(" ")
            chrom = transposons[0].replace(">","")
            if(not chrom in transposonAnnotations):
                transposonAnnotations[chrom] = list()
            counter += 1
            annoSoftware = "sineFinder"
            features     = "transposon"
            strand = "-"
            start  = int(transposons[2].split(":")[0])
            start    = lengthDict[chrom] - int(start)
            end    = int(transposons[2].split(":")[1])
            end    = lengthDict[chrom] - int(end)
            start,end = correctPositions(start, end)
            score  = "."
            phase  = "."
            transpNr = str(counter)
            attributes = ""
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            features     = "tsd"
            start  = int(transposons[2].split(":")[0])
            start    = lengthDict[chrom] - int(start)
            end    = int(transposons[2].split(":")[0])+int(transposons[3].split(";")[0].split("=")[1])
            end    = lengthDict[chrom] - int(end)
            start,end = correctPositions(start, end)
            score  = int(transposons[3].split(";")[1].split("=")[1])
            phase  = "."
            transpNr = str(counter)
            attributes = "Left TSD of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            features     = "tsd"
            start  = int(transposons[2].split(":")[1])-int(transposons[3].split(";")[0].split("=")[1])
            start    = lengthDict[chrom] - int(start)
            end    = int(transposons[2].split(":")[1])
            end    = lengthDict[chrom] - int(end)
            start,end = correctPositions(start, end)
            score  = int(transposons[3].split(";")[1].split("=")[1])
            phase  = "."
            transpNr = str(counter)
            attributes = "Right TSD of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        line = f.readline()
    f.close()
    # Print results to Annotation file
    f = open(targetGFFFile,"w+")
    f.write("# sineFinder Annotation")
    writeGFFhead(f)
    keys = list(transposonAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()
    # Get Transposon Sequences
    exportFastaFile(targetFastaFile, fastaFile, keys, transposonAnnotations)

def parseLtrpred(pathResDir, fastaFile, targetGFFFile, targetFastaFile):
    print("Parse ltrPred outputs...")
    transposonAnnotations = {}
    folder = os.listdir(pathResDir)[0]
    files  = os.listdir(os.path.join(pathResDir,folder))
    selFile = ""
    for f in files:
        if(f.endswith(".gff")):
            selFile = os.path.join(pathResDir,folder,f)
    f = open(selFile,"r")
    line = f.readline()
    counter = 0
    while line!="":
        transposons = line.replace("\n","").split("\t")
        chrom = transposons[0].replace(">","")
        if(not chrom in transposonAnnotations):
            transposonAnnotations[chrom] = list()
        counter += 1
        annoSoftware = "LTRpred"
        features     = "transposon"
        strand = transposons[6]
        start  = int(transposons[3])
        end    = int(transposons[4])
        start,end = correctPositions(start, end)
        score  = "."
        phase  = "."
        transpNr = str(counter)
        comments = transposons[8].split(";")
        attributes = comments[31]+";"+comments[32]
        transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        if(comments[3].split("=")[1]!="NA"):
            features     = "ltr"
            strand = transposons[6]
            start  = int(comments[3].split("=")[1])
            end    = int(comments[4].split("=")[1])
            start,end = correctPositions(start, end)
            attributes = "Left LTR of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        if(comments[6].split("=")[1]!="NA"):
            features     = "ltr"
            strand = transposons[6]
            start  = int(comments[6].split("=")[1])
            end    = int(comments[7].split("=")[1])
            start,end = correctPositions(start, end)
            attributes = "Right LTR of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes) 
        if(comments[9].split("=")[1]!="NA"):
            features     = "tsd"
            strand = transposons[6]
            start  = int(comments[9].split("=")[1])
            end    = int(comments[10].split("=")[1])
            start,end = correctPositions(start, end)
            attributes = "Left TSD of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        if(comments[12].split("=")[1]!="NA"):
            features     = "tsd"
            strand = transposons[6]
            start  = int(comments[12].split("=")[1])
            end    = int(comments[13].split("=")[1])
            start,end = correctPositions(start, end)
            attributes = "Right TSD of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        if(comments[15].split("=")[1]!="NA"):
            features     = "ppt"
            strand = transposons[6]
            start  = int(comments[15].split("=")[1])
            end    = int(comments[16].split("=")[1])
            start,end = correctPositions(start, end)
            attributes = "Polypurine tract of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        if(comments[21].split("=")[1]!="NA"):
            features     = "pbs"
            strand = transposons[6]
            start  = int(comments[21].split("=")[1])
            end    = int(comments[22].split("=")[1])
            start,end = correctPositions(start, end)
            attributes = "Primer binding site of transposon "+str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        line = f.readline()
    f.close()
    # Print results to Annotation file
    f = open(targetGFFFile,"w+")
    f.write("# ltrPred Annotation")
    writeGFFhead(f)
    keys = list(transposonAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()
    # Get Transposon Sequences
    exportFastaFile(targetFastaFile, fastaFile, keys, transposonAnnotations)

def parseRepeatMasker(pathResDir, fastaFile, targetGFFFile, targetRepeatFile, targetFastaFile):
    print("Parse repeatMasker...")
    scheme = loadRepBaseTaxonomyScheme()
    transposonAnnotations = {}
    repeatAnnotations = {}
    f = open(os.path.join(pathResDir,getFile(pathResDir,".out")),"r")
    line = f.readline()
    line = f.readline()
    line = f.readline()
    line = f.readline()
    counter = 0
    while line!="":
        while "  " in line:
            line = line.replace("  "," ")
        line = line.lstrip().rstrip()
        transposons = line.replace("\n","").split(" ")
        chrom = transposons[4].replace(">","")
        if(not chrom in transposonAnnotations):
            transposonAnnotations[chrom] = list()
        if(not chrom in repeatAnnotations):
            repeatAnnotations[chrom] = list()
        annoSoftware = "repeatMasker"
        label = transposons[10]
        Lclass = getTransposonClass(label, scheme)
        Ltype  = getFeatureType(label, scheme)
        
        if(not Lclass==""):
            features     = "transposon"
            attributes = "class "+Lclass+" ("+label+")"
            transpNr = str(counter)
            counter += 1
        else:
            features     = Ltype
            attributes   = "repmasker ("+label+")"
            transpNr = "."
        strand = transposons[8]
        if(strand=="C"):
            strand = "-"
        start  = int(transposons[5])
        end    = int(transposons[6])
        start,end = correctPositions(start, end)
        attributes = attributes.replace("\n","").replace("\r","")
        score  = "."
        phase  = "."
        if(features=="transposon"):
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        else:
            repeatAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        line = f.readline()
    f.close()
    # Print results to Annotation file
    f = open(targetGFFFile,"w+")
    f.write("# repeatMasker Annotation")
    writeGFFhead(f)
    keys = list(transposonAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()
    # Print results to Repeat file
    f = open(targetRepeatFile,"w+")
    f.write("# repeatModeler Annotation")
    writeGFFhead(f)
    keys = list(repeatAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(repeatAnnotations[key])):
            f.write(repeatAnnotations[key][i])
            f.write("\n")
    f.close()
    # Get Transposon Sequences
    exportFastaFile(targetFastaFile, fastaFile, keys, transposonAnnotations)

def parseTransposonPSI(pathResDir, fastaFile, targetGFFFile, targetFastaFile):
    print("Parse transposonPSI outputs...")
    transposonAnnotations = {}
    f = open(os.path.join(pathResDir,getFile(pathResDir,".TPSI.allHits")),"r")
    line = f.readline()
    counter = 0
    while line!="":
        transposons = line.replace("\n","").split("\t")
        chrom = transposons[4].split("/")[-1].split(".")[0].replace("_pilon","").replace(">","")
        if(not chrom in transposonAnnotations):
            transposonAnnotations[chrom] = list()
        annoSoftware = "transposonPSI" 
        features     = "protein/domain"
        attributes = "protein("+transposons[0]+")"
        transpNr = "."
        counter += 1
        sStrand = transposons[17]
        if(sStrand=="Plus"):
            strand = "+"
        else:
            strand = "-"
        start  = int(transposons[8])
        end    = int(transposons[9])
        start,end = correctPositions(start, end)
        score  = "."
        phase  = "."
        transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        line = f.readline()
    f.close()
    # Print results to Annotation file
    f = open(targetGFFFile,"w+")
    f.write("# transposonPSI Annotation")
    writeGFFhead(f)
    keys = list(transposonAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()

def parseNCBICDD1000(pathResDir, fastaFile, targetGFFFile, targetFastaFile):
    print("Parse NCBICDD1000 outputs...")
    candidates = loadCandidates()
    transposonAnnotations = {}
    files = os.listdir(os.path.join(pathResDir,"temp"))
    for file in files:
        f = open(os.path.join(pathResDir,"temp",file),"r")
        line = f.readline()
        chrom = f.readline().replace("\n","").split(":")[1].rstrip().lstrip().split("=")[0].replace("length","").lstrip().rstrip().replace(">","")
        if(not chrom in transposonAnnotations):
            transposonAnnotations[chrom] = list()
        line = f.readline()
        while True:
            while line!="" and line.startswith("#") and not line.startswith("# BLAST processed") and not line.startswith("# RPSTBLASTN"):
                line = f.readline()
            while line!="" and not line.startswith("#"):
                transposons = line.replace("\n","").split("\t")
                code = transposons[0].split(",")[0].rstrip().lstrip()
                annoSoftware = "NCBI_CDD_Selection" 
                features     = "protein"
                attributes = "protein("+code+"), "+candidates[code]+" "+transposons[0].split(",")[1].rstrip().lstrip()
                transpNr = "."
                if(int(transposons[2])<int(transposons[3])):
                    strand = "+"
                else:
                    strand = "-"
                start  = int(transposons[2])
                end    = int(transposons[3])
                start,end = correctPositions(start, end)
                score  = "."
                phase  = "."
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
    f = open(targetGFFFile,"w+")
    f.write("# NCBIsel Annotation")
    writeGFFhead(f)
    keys = list(transposonAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()
    
def parseMiteTracker(pathResDir, pathResDirRC, fastaFile, targetGFFFile, targetFastaFile):
    print("Parse miteTracker...")
    # Get Records lengths
    lengthDict = {}
    for record in SeqIO.parse(fastaFile, "fasta"):
        print("%s %i" % (record.id, len(record)))
        lengthDict[record.id] = len(record)
    # Load Annotations
    transposonAnnotations = {}
    counter = 0
    folder = os.listdir(pathResDir)[0]
    folder2 = os.listdir(os.path.join(pathResDir,folder))[0]
    f = open(os.path.join(pathResDir,folder,folder2,"all.gff3"),"r")
    line = f.readline()
    while line!="":
        if(not line.startswith("#")):
            transposons = line.replace("\n","").split("\t")
            chrom = transposons[0].replace(">","")
            if(not chrom in transposonAnnotations):
                transposonAnnotations[chrom] = list()
            annoSoftware = "miteTracker"
            features     = "transposon"
            counter += 1
            strand = transposons[6]
            start  = int(transposons[3])
            end    = int(transposons[4])
            start,end = correctPositions(start, end)
            score  = "."
            phase  = "."
            attributes = ""
            transpNr = str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        line = f.readline()
    f.close()
    folder = os.listdir(pathResDirRC)[0]
    folder2 = os.listdir(os.path.join(pathResDirRC,folder))[0]
    f = open(os.path.join(pathResDirRC,folder,folder2,"all.gff3"),"r")
    line = f.readline()
    while line!="":
        if(not line.startswith("#")):
            transposons = line.replace("\n","").split("\t")
            chrom = transposons[0].replace(">","")
            if(not chrom in transposonAnnotations):
                transposonAnnotations[chrom] = list()
            annoSoftware = "miteTracker"
            features     = "transposon"
            counter += 1
            strand = "-"
            start  = int(transposons[3])
            start  = lengthDict[chrom] - int(start)
            end    = int(transposons[4])
            end    = lengthDict[chrom] - int(end)
            start,end = correctPositions(start, end)
            score  = "."
            phase  = "."
            attributes = ""
            transpNr = str(counter)
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        line = f.readline()
    f.close()
    # Print results to Annotation file
    f = open(targetGFFFile,"w+")
    f.write("# miteTracker Annotation")
    writeGFFhead(f)
    keys = list(transposonAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()
    # Get Transposon Sequences
    exportFastaFile(targetFastaFile, fastaFile, keys, transposonAnnotations)

def parseSineScan(pathResDir, fastaFile, targetGFFFile, targetFastaFile):
    print("Parse sineScan...")
    transposonAnnotations = {}
#    files = ["chromosomes-matches.fasta","chromosomes-5smatches.fasta"]
    files = [getFile(os.path.join(pathResDir,"result"),"-matches.fasta"), getFile(os.path.join(pathResDir,"result"),"-5smatches.fasta")]
    counter = 0
    for file in files:
        f = open(os.path.join(pathResDir,"result",file),"r")
        line = f.readline()
        while line!="":
            if(line.startswith(">")):
                transposons = line.replace("\n","").split(" ")
                chrom = transposons[0].replace(">","")
                if(not chrom in transposonAnnotations):
                    transposonAnnotations[chrom] = list()
                annoSoftware = "sineScan"
                features     = "transposon"
                counter += 1
                sStrand = transposons[1]
                if(sStrand=="R"):
                    strand = "+"
                elif(sStrand=="F"):
                    strand = "-"
                start  = int(transposons[2].split(":")[0])
                end    = int(transposons[2].split(":")[1])
                start,end = correctPositions(start, end)
                score  = "."
                phase  = "."
                attributes = ""
                transpNr = str(counter)
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
                features     = "tsd"
                sStrand = transposons[1]
                if(sStrand=="R"):
                    strand = "+"
                elif(sStrand=="F"):
                    strand = "-"
                start  = int(transposons[2].split(":")[0]) - int(transposons[3].split(";")[0].split("=")[1])
                end    = int(transposons[2].split(":")[0])
                start,end = correctPositions(start, end)
                score  = "."
                phase  = "."
                attributes = "Left TSD of transposon "+str(counter)
                transpNr = str(counter)
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
                features     = "tsd"
                sStrand = transposons[1]
                if(sStrand=="R"):
                    strand = "+"
                elif(sStrand=="F"):
                    strand = "-"
                start  = int(transposons[2].split(":")[1])
                end    = int(transposons[2].split(":")[1]) + int(transposons[3].split(";")[0].split("=")[1])
                start,end = correctPositions(start, end)
                score  = "."
                phase  = "."
                attributes = "Right TSD of transposon "+str(counter)
                transpNr = str(counter)
                transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
            line = f.readline()
        f.close()
    # Print results to Annotation file
    f = open(targetGFFFile,"w+")
    f.write("# sineScan Annotation")
    writeGFFhead(f)
    keys = list(transposonAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()
    # Get Transposon Sequences
    exportFastaFile(targetFastaFile, fastaFile, keys, transposonAnnotations)

def parseRepeatModeler(pathResDir, fastaFile, targetGFFFile, targetRepeatFile, targetFastaFile):
    print("Parse repeatModeler...")
    transposonAnnotations = {}
    repeatAnnotations = {}
    counter = 0
    f = open(os.path.join(pathResDir,getFile(pathResDir,"-families.stk")),"r")
    line = f.readline()
    seqName = ""
    seqType = ""
    while line!="":
        if(line.startswith("//")):
            line = f.readline()
        if(line.startswith("# STOCKHOLM")):
            line = f.readline()
            head = list()
            while line.startswith("#"):
                head.append(line)
                line = f.readline()
            for h in head:
                if(h.startswith("#=GF ID")):
                    seqName = h.replace("\n","").split(" ")[-1]
                if(h.startswith("#=GF TP")):
                    seqType = h.replace("\n","").split(" ")[-1]    
        transposons = line.replace("\n","").split(" ")
        chrom = transposons[0].split(":")[0].replace(">","")
        if(not chrom in transposonAnnotations):
            transposonAnnotations[chrom] = list()
        if(not chrom in repeatAnnotations):
            repeatAnnotations[chrom] = list()
        if(line==""):
            break    
        annoSoftware = "repeatmodel"
        strand = "+"
        counter += 1
        start  = int(transposons[0].split(":")[1].split("-")[0])
        end    = int(transposons[0].split(":")[1].split("-")[1])
#        start,end = correctPositions(start, end)
        if(end<start):
            strand = "-"
            c = start
            start = end
            end = c
        score  = "."
        phase  = "."
        transpNr = str(counter)
        attributes = seqName+"("+seqType+")"    
        seqTypeLabelA = seqType.split(";")[1]
        if(seqTypeLabelA.startswith("Transposable")):
            features = "transposon"
            classT  = getRepeatModelTransposonType(seqType)
            attributes = "class "+classT+" "+seqName+"("+seqType+")"   
            attributes = attributes.replace("\n","").replace("\r","")
            transposonAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        else:
            attributes = ""
            features = "interspersed_repeat"
            repeatAnnotations[chrom].append(chrom+"\t"+annoSoftware+"\t"+features+"\t"+str(start)+"\t"+str(end)+"\t"+str(score)+"\t"+strand+"\t"+phase+"\t"+transpNr+"\t"+attributes)
        line = f.readline()
    f.close()
    # Print results to Annotation file
    f = open(targetGFFFile,"w+")
    f.write("# repeatModeler Annotation")
    writeGFFhead(f)
    keys = list(transposonAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(transposonAnnotations[key])):
            f.write(transposonAnnotations[key][i])
            f.write("\n")
    f.close()
    # Print results to Repeat file
    f = open(targetRepeatFile,"w+")
    f.write("# repeatModeler Annotation")
    writeGFFhead(f)
    keys = list(repeatAnnotations.keys())
    keys.sort()
    for key in keys:
        for i in range(0,len(repeatAnnotations[key])):
            f.write(repeatAnnotations[key][i])
            f.write("\n")
    f.close()
    # Get Transposon Sequences
    exportFastaFile(targetFastaFile, fastaFile, keys, transposonAnnotations)

# Main Code
def parseAvailableResults(projectFolderPath):    
    if(checkHelitronScanner(projectFolderPath)):
        pathResDir      = os.path.join(projectFolderPath,"helitronScanner")
        pathResDirRC    = os.path.join(projectFolderPath,"helitronScanner_rc")
        fastaFile       = os.path.join(projectFolderPath,"sequence.fasta")
        targetGFFFile   = os.path.join(projectFolderPath,"parsedAnnotations","helitronScanner.gff3")
        targetFastaFile = os.path.join(projectFolderPath,"parsedAnnotations","helitronScanner.fasta")
        parseHelitronScanner(pathResDir, pathResDirRC, fastaFile, targetGFFFile, targetFastaFile)
    if(checkLtrHarvest(projectFolderPath)):
        pathResDir      = os.path.join(projectFolderPath,"ltrHarvest")
        fastaFile       = os.path.join(projectFolderPath,"sequence.fasta")
        targetGFFFile   = os.path.join(projectFolderPath,"parsedAnnotations","ltrHarvest.gff3")
        targetFastaFile = os.path.join(projectFolderPath,"parsedAnnotations","ltrHarvest.fasta")
        parseLtrHarvest(pathResDir, fastaFile, targetGFFFile, targetFastaFile)
    if(checkLtrPred(projectFolderPath)):
        pathResDir      = os.path.join(projectFolderPath,"ltrPred")
        fastaFile       = os.path.join(projectFolderPath,"sequence.fasta")
        targetGFFFile   = os.path.join(projectFolderPath,"parsedAnnotations","ltrPred.gff3")
        targetFastaFile = os.path.join(projectFolderPath,"parsedAnnotations","ltrPred.fasta")
        parseLtrpred(pathResDir, fastaFile, targetGFFFile, targetFastaFile)
    if(checkMiteFind(projectFolderPath)):
        pathResDir      = os.path.join(projectFolderPath,"mitefind")
        pathResDirRC    = os.path.join(projectFolderPath,"mitefind_rc")
        fastaFile       = os.path.join(projectFolderPath,"sequence.fasta")
        targetGFFFile   = os.path.join(projectFolderPath,"parsedAnnotations","mitefind.gff3")
        targetFastaFile = os.path.join(projectFolderPath,"parsedAnnotations","mitefind.fasta")
        parseMiteFind(pathResDir, pathResDirRC, fastaFile, targetGFFFile, targetFastaFile)
    if(checkMiteTracker(projectFolderPath)):
        pathResDir      = os.path.join(projectFolderPath,"mitetracker")
        pathResDirRC    = os.path.join(projectFolderPath,"mitetracker_rc")
        fastaFile       = os.path.join(projectFolderPath,"sequence.fasta")
        targetGFFFile   = os.path.join(projectFolderPath,"parsedAnnotations","mitetracker.gff3")
        targetFastaFile = os.path.join(projectFolderPath,"parsedAnnotations","mitetracker.fasta")
        parseMiteTracker(pathResDir, pathResDirRC, fastaFile, targetGFFFile, targetFastaFile)
    if(checkMust(projectFolderPath)):
        pathResDir      = os.path.join(projectFolderPath,"must")
        fastaFile       = os.path.join(projectFolderPath,"sequence.fasta")
        targetGFFFile   = os.path.join(projectFolderPath,"parsedAnnotations","must.gff3")
        targetFastaFile = os.path.join(projectFolderPath,"parsedAnnotations","must.fasta")
        parseMust(pathResDir, fastaFile, targetGFFFile, targetFastaFile)
    if(checkNCBICDD1000(projectFolderPath)):
        pathResDir      = os.path.join(projectFolderPath,"NCBICDD1000")
        fastaFile       = os.path.join(projectFolderPath,"sequence.fasta")
        targetGFFFile   = os.path.join(projectFolderPath,"parsedAnnotations","NCBICDD1000.gff3")
        targetFastaFile = os.path.join(projectFolderPath,"parsedAnnotations","NCBICDD1000.fasta")   
        parseNCBICDD1000(pathResDir, fastaFile, targetGFFFile, targetFastaFile)
    if(checkRepeatModeler(projectFolderPath)):
        pathResDir      = os.path.join(projectFolderPath,"repeatmodel")
        fastaFile       = os.path.join(projectFolderPath,"sequence.fasta")
        targetGFFFile   = os.path.join(projectFolderPath,"parsedAnnotations","repeatmodel.gff3")
        targetGFFrepe   = os.path.join(projectFolderPath,"parsedAnnotations","repeatmodel_repeats.gff3")
        targetFastaFile = os.path.join(projectFolderPath,"parsedAnnotations","repeatmodel.fasta")   
        parseRepeatModeler(pathResDir, fastaFile, targetGFFFile, targetGFFrepe, targetFastaFile)
    if(checkRepeatMasker(projectFolderPath)):
        pathResDir      = os.path.join(projectFolderPath,"repMasker")
        fastaFile       = os.path.join(projectFolderPath,"sequence.fasta")
        targetGFFFile   = os.path.join(projectFolderPath,"parsedAnnotations","repMasker.gff3")
        targetGFFrepe   = os.path.join(projectFolderPath,"parsedAnnotations","repMasker_repeats.gff3")
        targetFastaFile = os.path.join(projectFolderPath,"parsedAnnotations","repMasker.fasta")   
        parseRepeatMasker(pathResDir, fastaFile, targetGFFFile, targetGFFrepe, targetFastaFile)
    if(checkSineFind(projectFolderPath)):
        pathResDir      = os.path.join(projectFolderPath,"sinefind")
        pathResDirRC    = os.path.join(projectFolderPath,"sinefind_rc")
        fastaFile       = os.path.join(projectFolderPath,"sequence.fasta")
        targetGFFFile   = os.path.join(projectFolderPath,"parsedAnnotations","sinefind.gff3")
        targetFastaFile = os.path.join(projectFolderPath,"parsedAnnotations","sinefind.fasta")
        parseSinefind(pathResDir, pathResDirRC, fastaFile, targetGFFFile, targetFastaFile)
    if(checkSineScan(projectFolderPath)):
        pathResDir      = os.path.join(projectFolderPath,"sinescan")
        fastaFile       = os.path.join(projectFolderPath,"sequence.fasta")
        targetGFFFile   = os.path.join(projectFolderPath,"parsedAnnotations","sinescan.gff3")
        targetFastaFile = os.path.join(projectFolderPath,"parsedAnnotations","sinescan.fasta")   
        parseSineScan(pathResDir, fastaFile, targetGFFFile, targetFastaFile)
    if(checkTirVish(projectFolderPath)):
        pathResDir      = os.path.join(projectFolderPath,"tirvish")
        pathResDirRC    = os.path.join(projectFolderPath,"tirvish_rc")
        fastaFile       = os.path.join(projectFolderPath,"sequence.fasta")
        targetGFFFile   = os.path.join(projectFolderPath,"parsedAnnotations","tirvish.gff3")
        targetFastaFile = os.path.join(projectFolderPath,"parsedAnnotations","tirvish.fasta")
        parseTirvish(pathResDir, pathResDirRC, fastaFile, targetGFFFile, targetFastaFile)
    if(checkTransposonPSI(projectFolderPath)):
        pathResDir      = os.path.join(projectFolderPath,"transposonPSI")
        fastaFile       = os.path.join(projectFolderPath,"sequence.fasta")
        targetGFFFile   = os.path.join(projectFolderPath,"parsedAnnotations","transposonPSI.gff3")
        targetFastaFile = os.path.join(projectFolderPath,"parsedAnnotations","transposonPSI.fasta")   
        parseTransposonPSI(pathResDir, fastaFile, targetGFFFile, targetFastaFile)
        
#projectFolder = "projects_C_Elegans"
#projectName = "PRJEB28388"
#projectFolderPath = os.path.abspath(os.path.join(projectFolder,projectName))
#parseAvailableResults(projectFolderPath)
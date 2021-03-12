############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
import os
from os import path
import os.path

# Methods
def getArgument(args, title):
    for i in range(0, len(args)):
        if(args[i].startswith("-"+title)):
            if(i<len(args)-1):
                return args[i+1]
            else:
                return ""
    return ""

def getFile(folder, suffix):
    files = os.listdir(folder)
    for f in files:
        if(f.endswith(suffix)):
            return f
    return ""

def checkHelitronScanner(projectFolderPath):
    if path.isdir(os.path.join(projectFolderPath, "helitronScanner")) and path.isdir(os.path.join(projectFolderPath, "helitronScanner_rc")):
        if(path.isfile(os.path.join(projectFolderPath, "helitronScanner", "result.txt"))) and path.isfile(os.path.join(projectFolderPath, "helitronScanner_rc", "result.txt")):
            return True
    return False
    
def checkLtrHarvest(projectFolderPath):
    if path.isdir(os.path.join(projectFolderPath, "ltrHarvest")):
        if(getFile(os.path.join(projectFolderPath, "ltrHarvest"),".txt")!=""):
            return True
#        if(path.isfile(os.path.join(projectFolderPath, "ltrHarvest", "result.txt"))):
#            return True
    return False

def checkLtrPred(projectFolderPath):
    if path.isdir(os.path.join(projectFolderPath, "ltrPred")):
        folders = os.listdir(os.path.join(projectFolderPath, "ltrPred"))
        found = False
        folder = ""
        for f in folders:
            if(path.isdir(os.path.join(projectFolderPath, "ltrPred", f))):
                found = True
                folder = f
        if(found):
            files = os.listdir(os.path.join(projectFolderPath, "ltrPred",folder))
            for f in files:
                if(f.endswith(".gff")):
                    return True
    return False

def checkMiteFind(projectFolderPath):
    if path.isdir(os.path.join(projectFolderPath, "mitefind")) and path.isdir(os.path.join(projectFolderPath, "mitefind_rc")):
        if(path.isfile(os.path.join(projectFolderPath, "mitefind","result.txt")) and path.isfile(os.path.join(projectFolderPath, "mitefind_rc","result.txt"))):
            return True
    return False
    
def checkMiteTracker(projectFolderPath):
    if path.isdir(os.path.join(projectFolderPath, "mitetracker")) and path.isdir(os.path.join(projectFolderPath, "mitetracker_rc")):
        if path.isdir(os.path.join(projectFolderPath, "mitetracker", "results")) and path.isdir(os.path.join(projectFolderPath, "mitetracker_rc", "results")):
            if path.isdir(os.path.join(projectFolderPath, "mitetracker", "results", "job")) and path.isdir(os.path.join(projectFolderPath, "mitetracker_rc", "results", "job")):
                if path.isfile(os.path.join(projectFolderPath, "mitetracker", "results", "job", "all.gff3")) and path.isfile(os.path.join(projectFolderPath, "mitetracker_rc", "results", "job", "all.gff3")):
                    return True
    return False

def checkMust(projectFolderPath):
    if path.isdir(os.path.join(projectFolderPath, "must")):
        if path.isfile(os.path.join(projectFolderPath, "must", "result.txt")):
            return True
    return False

def checkNCBICDD1000(projectFolderPath):
    if path.isdir(os.path.join(projectFolderPath, "NCBICDD1000")):
        if path.isdir(os.path.join(projectFolderPath, "NCBICDD1000", "temp")):
            files = os.listdir(os.path.join(projectFolderPath, "NCBICDD1000", "temp"))
            if(len(files)==50):
                return True
    return False

def checkRepeatModeler(projectFolderPath):
    if path.isdir(os.path.join(projectFolderPath, "repeatmodel")):
        files = os.listdir(os.path.join(projectFolderPath, "repeatmodel"))
        for f in files:
            if(f.endswith("-families.stk")):
                return True
    return False

def checkRepeatMasker(projectFolderPath):
    if path.isdir(os.path.join(projectFolderPath, "repMasker")):
        files = os.listdir(os.path.join(projectFolderPath, "repMasker"))
        for f in files:
            if(f.endswith(".out")):
                return True
    return False

def checkSineFind(projectFolderPath):
    if path.isdir(os.path.join(projectFolderPath, "sinefind")) and path.isdir(os.path.join(projectFolderPath, "sinefind_rc")):
        filesA = os.listdir(os.path.join(projectFolderPath, "sinefind"))
        filesB = os.listdir(os.path.join(projectFolderPath, "sinefind_rc"))
        foundA = False
        foundB = False
        for f in filesA:
            if(f.endswith("-matches.fasta")):
                foundA=True
                break
        for f in filesB:
            if(f.endswith("-matches.fasta")):
                foundB=True
                break
        if(foundA and foundB):
            return True
    return False

def checkSineScan(projectFolderPath):
    if path.isdir(os.path.join(projectFolderPath, "sinescan")):
        if path.isdir(os.path.join(projectFolderPath, "sinescan", "result")):
            files = os.listdir(os.path.join(projectFolderPath, "sinescan", "result"))
            foundA = False
            foundB = False
            for f in files:
                if(f.endswith("-matches.fasta")):
                    foundA = True
                if(f.endswith("-5smatches.fasta")):
                    foundB = True
            if(foundA and foundB):
                return True
    return False

def checkTirVish(projectFolderPath):
    if path.isdir(os.path.join(projectFolderPath, "tirvish")) and path.isdir(os.path.join(projectFolderPath, "tirvish_rc")):
        if(path.isfile(os.path.join(projectFolderPath, "tirvish","result.txt")) and path.isfile(os.path.join(projectFolderPath, "tirvish_rc","result.txt"))):
            return True
    return False

def checkTransposonPSI(projectFolderPath):
    if path.isdir(os.path.join(projectFolderPath, "transposonPSI")):
        files = os.listdir(os.path.join(projectFolderPath, "transposonPSI"))
        for f in files:
            if(f.endswith(".TPSI.allHits")):
                return True
    return False

def checkAnnotations(projectFolder, projectName):
    projectFolderPath = os.path.abspath(os.path.join(projectFolder,projectName))
    listOkay = list()
    for tool in ["helitronScanner", "ltrHarvest", "ltrPred", "mitefind", "mitetracker", "must", "repeatmodel", "repMasker", "sinefind", "sinescan", "tirvish", "transposonPSI", "NCBICDD1000"]:
        status = False
        if(tool=="helitronScanner"):
            status = checkHelitronScanner(projectFolderPath)
        elif(tool=="ltrHarvest"):
            status = checkLtrHarvest(projectFolderPath)
        elif(tool=="mitefind"):
            status = checkMiteFind(projectFolderPath)
        elif(tool=="mitetracker"):
            status = checkMiteTracker(projectFolderPath)
        elif(tool=="must"):
            status = checkMust(projectFolderPath)
        elif(tool=="repeatmodel"):
            status = checkRepeatModeler(projectFolderPath)
        elif(tool=="repMasker"):
            status = checkRepeatMasker(projectFolderPath)
        elif(tool=="sinefind"):
            status = checkSineFind(projectFolderPath)
        elif(tool=="sinescan"):
            status = checkSineScan(projectFolderPath)
        elif(tool=="tirvish"):
            status = checkTirVish(projectFolderPath)
        elif(tool=="transposonPSI"):
            status = checkTransposonPSI(projectFolderPath)
        elif(tool=="NCBICDD1000"):
            status = checkNCBICDD1000(projectFolderPath)
        elif(tool=="ltrPred"):
            status = checkLtrPred(projectFolderPath)
        if(status):
            print("Checking "+tool+"\t... completed")
            listOkay.append(tool)
        else:
            print("Checking "+tool+"\t... not completed")
    return listOkay

#args = sys.argv
#projectFolder = getArgument(args, "projectFolder")
#projectName   = getArgument(args, "projectName")
#checkAnnotations(projectFolder, projectName)
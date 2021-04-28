from __future__ import print_function
############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
import os
from shutil import copyfile, rmtree
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

def runHelitronScanner(projectFolderPath):
    os.system("helitronscanner scanHead -g "+os.path.join(projectFolderPath,"sequence.fasta")+" -bs 0 -o "+os.path.join(projectFolderPath,"helitronScanner","scanHead.txt"))
    os.system("helitronscanner scanTail -g "+os.path.join(projectFolderPath,"sequence.fasta")+" -bs 0 -o "+os.path.join(projectFolderPath,"helitronScanner","scanTail.txt"))
    os.system("helitronscanner pairends -hs "+os.path.join(projectFolderPath,"helitronScanner","scanHead.txt")+" -ts "+os.path.join(projectFolderPath,"helitronScanner","scanTail.txt")+" -o "+os.path.join(projectFolderPath,"helitronScanner","result.txt"))
    os.system("helitronscanner scanHead -g "+os.path.join(projectFolderPath,"sequence_rc.fasta")+" -bs 0 -o "+os.path.join(projectFolderPath,"helitronScanner_rc","scanHead.txt"))
    os.system("helitronscanner scanTail -g "+os.path.join(projectFolderPath,"sequence_rc.fasta")+" -bs 0 -o "+os.path.join(projectFolderPath,"helitronScanner_rc","scanTail.txt"))
    os.system("helitronscanner pairends -hs "+os.path.join(projectFolderPath,"helitronScanner_rc","scanHead.txt")+" -ts "+os.path.join(projectFolderPath,"helitronScanner_rc","scanTail.txt")+" -o "+os.path.join(projectFolderPath,"helitronScanner_rc","result.txt"))

def runMust(projectFolderPath, addCommand):
    os.mkdir(os.path.join(projectFolderPath,"must","temp"))
#    os.system("mustv2 "+os.path.join(projectFolderPath,"sequence.fasta")+" "+os.path.join(projectFolderPath,"must","result.txt")+" "+os.path.join(projectFolderPath,"must","temp"))
    os.system("cd "+os.path.join(projectFolderPath,"must")+" && mustv2 "+os.path.join(projectFolderPath,"sequence.fasta")+" result.txt temp "+addCommand)

def runSineFinder(projectFolderPath, addCommand):
    copyfile(os.path.join(projectFolderPath,"sequence.fasta"),os.path.join(projectFolderPath,"sinefind","sequence.fasta"))
    copyfile(os.path.join(projectFolderPath,"sequence_rc.fasta"),os.path.join(projectFolderPath,"sinefind_rc","sequence_rc.fasta"))
    if(addCommand==""):
        os.system("sine_finder -V "+os.path.join(projectFolderPath,"sinefind","sequence.fasta"))
        os.system("sine_finder -V "+os.path.join(projectFolderPath,"sinefind_rc","sequence_rc.fasta"))
    else:
        os.system("sine_finder "+addCommand+" "+os.path.join(projectFolderPath,"sinefind","sequence.fasta"))
        os.system("sine_finder "+addCommand+" "+os.path.join(projectFolderPath,"sinefind_rc","sequence_rc.fasta"))
    os.remove(os.path.join(projectFolderPath,"sinefind","sequence.fasta")) 
    os.remove(os.path.join(projectFolderPath,"sinefind_rc","sequence_rc.fasta")) 
    
def runMiteTracker(projectFolderPath, addCommand):
    if(path.isdir(os.path.join(projectFolderPath,"mitetracker","results"))):
        rmtree(os.path.join(projectFolderPath,"mitetracker","results"))
    if(path.isdir(os.path.join(projectFolderPath,"mitetracker_rc","results"))):
        rmtree(os.path.join(projectFolderPath,"mitetracker_rc","results"))
    os.mkdir(os.path.join(projectFolderPath,"mitetracker","results"))
    if(addCommand==""):
        os.system("cd "+os.path.join(projectFolderPath,"mitetracker")+" && mitetracker -g "+os.path.join(projectFolderPath,"sequence.fasta")+" -j job -w 4 ")
    else:
        os.system("cd "+os.path.join(projectFolderPath,"mitetracker")+" && mitetracker -g "+os.path.join(projectFolderPath,"sequence.fasta")+" -j job "+addCommand)
    os.mkdir(os.path.join(projectFolderPath,"mitetracker_rc","results"))
    if(addCommand==""):
        os.system("cd "+os.path.join(projectFolderPath,"mitetracker_rc")+" && mitetracker -g "+os.path.join(projectFolderPath,"sequence_rc.fasta")+" -j job -w 4 ")
    else:
        os.system("cd "+os.path.join(projectFolderPath,"mitetracker_rc")+" && mitetracker -g "+os.path.join(projectFolderPath,"sequence_rc.fasta")+" -j job -w 4 "+addCommand)

def runMiteFinder(projectFolderPath, addCommand):
    if(addCommand==""):
        os.system("miteFinderII -input "+os.path.join(projectFolderPath,"sequence.fasta")+" -output "+os.path.join(projectFolderPath,"mitefind","result.txt")+" "+addCommand)
        os.system("miteFinderII -input "+os.path.join(projectFolderPath,"sequence_rc.fasta")+" -output "+os.path.join(projectFolderPath,"mitefind_rc","result.txt")+" "+addCommand)
    else:
        os.system("miteFinderII -input "+os.path.join(projectFolderPath,"sequence.fasta")+" -output "+os.path.join(projectFolderPath,"mitefind","result.txt")+" -threshold 0.5 "+addCommand)
        os.system("miteFinderII -input "+os.path.join(projectFolderPath,"sequence_rc.fasta")+" -output "+os.path.join(projectFolderPath,"mitefind_rc","result.txt")+" -threshold 0.5 "+addCommand)

def runSineScan(projectFolderPath, addCommand):
    if(path.isdir(os.path.join(projectFolderPath,"sinescan","result"))):
        rmtree(os.path.join(projectFolderPath,"sinescan","result"))
    if(path.isdir(os.path.join(projectFolderPath,"sinescan","output"))):
        rmtree(os.path.join(projectFolderPath,"sinescan","output"))
    if(path.isdir(os.path.join(projectFolderPath,"sinescan","final"))):
        rmtree(os.path.join(projectFolderPath,"sinescan","final"))
    os.mkdir(os.path.join(projectFolderPath,"sinescan","result"))
    os.mkdir(os.path.join(projectFolderPath,"sinescan","output"))
    os.mkdir(os.path.join(projectFolderPath,"sinescan","final"))
    os.system("sinescan -s 123 -g "+os.path.join(projectFolderPath,"sequence.fasta")+" -o "+os.path.join(projectFolderPath,"sinescan","output")+" -d "+os.path.join(projectFolderPath,"sinescan","result")+" -z "+os.path.join(projectFolderPath,"sinescan","final")+" "+addCommand)

def runTirVish(projectFolderPath, addCommand):
    os.system("gt suffixerator -db "+os.path.join(projectFolderPath,"sequence.fasta")+" -indexname "+os.path.join(projectFolderPath,"tirvish","sequence.index")+" -tis -suf -lcp -des -ssp -sds -dna -mirrored")
    os.system("gt tirvish -index "+os.path.join(projectFolderPath,"tirvish","sequence.index")+" > "+os.path.join(projectFolderPath,"tirvish","result.txt")+" "+addCommand)
    os.system("gt suffixerator -db "+os.path.join(projectFolderPath,"sequence_rc.fasta")+" -indexname "+os.path.join(projectFolderPath,"tirvish_rc","sequence.index")+" -tis -suf -lcp -des -ssp -sds -dna -mirrored")
    os.system("gt tirvish -index "+os.path.join(projectFolderPath,"tirvish_rc","sequence.index")+" > "+os.path.join(projectFolderPath,"tirvish_rc","result.txt")+" "+addCommand)

def runLtrHarvest(projectFolderPath, addCommand):
    os.system("gt suffixerator -db "+os.path.join(projectFolderPath,"sequence.fasta")+" -indexname "+os.path.join(projectFolderPath,"ltrHarvest","sequence.index")+" -tis -suf -lcp -des -ssp -sds -dna")
    os.system("gt ltrharvest -index "+os.path.join(projectFolderPath,"ltrHarvest","sequence.index")+" > "+os.path.join(projectFolderPath,"ltrHarvest","result.txt")+" "+addCommand)

def runRepeatModeler(projectFolderPath, addCommand):
    copyfile(os.path.join(projectFolderPath,"sequence.fasta"), os.path.join(projectFolderPath, "repeatmodel", "sequence.fasta"))
    os.system("cd "+os.path.join(projectFolderPath,"repeatmodel")+" && BuildDatabase -name sequence_index -engine ncbi sequence.fasta")
    if(addCommand==""):
        os.system("cd "+os.path.join(projectFolderPath,"repeatmodel")+" && RepeatModeler -engine ncbi -pa 10 -database sequence_index")
    else:
        os.system("cd "+os.path.join(projectFolderPath,"repeatmodel")+" && RepeatModeler -engine ncbi -database sequence_index "+addCommand)
    os.remove(os.path.join(projectFolderPath, "repeatmodel", "sequence.fasta"))

def runRepeatMasker(projectFolderPath, addCommand):
    copyfile(os.path.join(projectFolderPath,"sequence.fasta"), os.path.join(projectFolderPath, "repMasker", "sequence.fasta"))
    if(addCommand==""):
        os.system("cd "+os.path.join(projectFolderPath,"repMasker")+" && RepeatMasker -pa 10 sequence.fasta")
    else:
        os.system("cd "+os.path.join(projectFolderPath,"repMasker")+" && RepeatMasker sequence.fasta "+addCommand)
    os.remove(os.path.join(projectFolderPath, "repMasker", "sequence.fasta"))

def runTransposonPSI(projectFolderPath):
    os.mkdir(os.path.join(projectFolderPath,"transposonPSI","temp"))
    os.system("transposonPSI -fastaFile "+os.path.join(projectFolderPath,"sequence.fasta")+" -resultFolder "+os.path.join(projectFolderPath,"transposonPSI")+" -tempFolder "+os.path.join(projectFolderPath,"transposonPSI","temp")+" -mode nuc")

def runNCBICDD1000(projectFolderPath):
    os.mkdir(os.path.join(projectFolderPath,"NCBICDD1000","temp"))
    os.system("proteinNCBICDD1000 -fastaFile "+os.path.join(projectFolderPath,"sequence.fasta")+" -resultFolder "+os.path.join(projectFolderPath,"NCBICDD1000","temp"))
    
def runAnnotation(projectFolder, projectName, tool, addCommand):
    projectFolderPath = os.path.abspath(os.path.join(projectFolder,projectName))
    if(tool=="helitronScanner"):
        runHelitronScanner(projectFolderPath)
    elif(tool=="ltrHarvest"):
        runLtrHarvest(projectFolderPath, addCommand)
    elif(tool=="mitefind"):
        runMiteFinder(projectFolderPath, addCommand)
    elif(tool=="mitetracker"):
        runMiteTracker(projectFolderPath, addCommand)
    elif(tool=="must"):
        runMust(projectFolderPath, addCommand)
    elif(tool=="repeatmodel"):
        runRepeatModeler(projectFolderPath, addCommand)
    elif(tool=="repMasker"):
        runRepeatMasker(projectFolderPath, addCommand)
    elif(tool=="sinefind"):
        runSineFinder(projectFolderPath, addCommand)
    elif(tool=="sinescan"):
        runSineScan(projectFolderPath, addCommand)
    elif(tool=="tirvish"):
        runTirVish(projectFolderPath, addCommand)
    elif(tool=="transposonPSI"):
        runTransposonPSI(projectFolderPath)
    elif(tool=="NCBICDD1000"):
        runNCBICDD1000(projectFolderPath)
    elif(tool=="all"):
        runHelitronScanner(projectFolderPath)
        runLtrHarvest(projectFolderPath, "")
        runMiteFinder(projectFolderPath, "")
        runMiteTracker(projectFolderPath, "")
        runMust(projectFolderPath, "")
        runRepeatModeler(projectFolderPath, "")
        runRepeatMasker(projectFolderPath, "")
        runSineFinder(projectFolderPath, "")
        runSineScan(projectFolderPath, "")
        runTirVish(projectFolderPath, "")
        runTransposonPSI(projectFolderPath)
        runNCBICDD1000(projectFolderPath)
    else:
        print("ERROR: tool not found!")
        return "EXIT"
    
# Parameters
#projectFolder = "projects"
#projectName = "testProject"
#tool = "ltrHarvest"
#projectFolderPath = os.path.abspath(os.path.join(projectFolder,projectName))
#print(projectFolderPath)

#args = sys.argv
#projectFolder = getArgument(args, "projectFolder")
#projectName   = getArgument(args, "projectName")
#tool          = getArgument(args, "tool")
#runAnnotation(projectFolder, projectName, tool)

############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from os import path
import os.path

# Methods
def make_rc_record(record):
    return SeqRecord(seq = record.seq.reverse_complement(), id = record.id, description="")

def copySequenceClean(fromFile,projectFolderPath):
    # Copy sequence and clean heads
    f1 = open(fromFile,"r")
    f2 = open(os.path.join(projectFolderPath,"sequence.fasta"),"w+")
    f3 = open(os.path.join(projectFolderPath,"sequence_heads.txt"),"w+")
    line = f1.readline()
    counter = 0
    while line!="":
        if(line.startswith(">")):
            counter += 1
            f3.write(">seq"+str(counter)+"\t"+line)
            f2.write(">seq"+str(counter)+"\n")
        else:
            f2.write(line.upper())
        line = f1.readline()
    f1.close()
    f2.close()
    f3.close()
    # Create reverse complement Fasta file
    records = map(make_rc_record, SeqIO.parse(os.path.join(projectFolderPath,"sequence.fasta"), "fasta"))
    SeqIO.write(records, os.path.join(projectFolderPath,"sequence_rc.fasta"), "fasta")
    records = map(make_rc_record, SeqIO.parse(os.path.join(projectFolderPath,"sequence_rc.fasta"), "fasta"))
    SeqIO.write(records, os.path.join(projectFolderPath,"sequence.fasta"), "fasta")
    
def createProject(projectFolder, projectName, inputFasta):
    # Check if project folder exists
    if(not path.isdir(projectFolder)):
        os.mkdir(projectFolder)    
    # Check if given project already exits
    projectFolderPath = os.path.join(projectFolder,projectName)
    if(path.isdir(projectFolderPath)):
        print("Project already exists, process aborted")
        return "EXIT"
    os.mkdir(projectFolderPath)
    # Create folder structure for annotation softwares
    os.mkdir(os.path.join(projectFolderPath,"tirvish"))
    os.mkdir(os.path.join(projectFolderPath,"tirvish_rc"))
    os.mkdir(os.path.join(projectFolderPath,"sinescan"))
    os.mkdir(os.path.join(projectFolderPath,"sinefind"))
    os.mkdir(os.path.join(projectFolderPath,"sinefind_rc"))
    os.mkdir(os.path.join(projectFolderPath,"repMasker"))
    os.mkdir(os.path.join(projectFolderPath,"repeatmodel"))
    os.mkdir(os.path.join(projectFolderPath,"must"))
    os.mkdir(os.path.join(projectFolderPath,"mitetracker"))
    os.mkdir(os.path.join(projectFolderPath,"mitetracker_rc"))
    os.mkdir(os.path.join(projectFolderPath,"mitefind"))
    os.mkdir(os.path.join(projectFolderPath,"mitefind_rc"))
    os.mkdir(os.path.join(projectFolderPath,"ltrPred"))
    os.mkdir(os.path.join(projectFolderPath,"ltrHarvest"))
    os.mkdir(os.path.join(projectFolderPath,"helitronScanner"))
    os.mkdir(os.path.join(projectFolderPath,"helitronScanner_rc")) 
    os.mkdir(os.path.join(projectFolderPath,"transposonPSI")) 
    os.mkdir(os.path.join(projectFolderPath,"NCBICDD1000")) 
    os.mkdir(os.path.join(projectFolderPath,"parsedAnnotations")) 
    os.mkdir(os.path.join(projectFolderPath,"transposonCandA")) 
    os.mkdir(os.path.join(projectFolderPath,"transposonCandB")) 
    os.mkdir(os.path.join(projectFolderPath,"transposonCandC")) 
    os.mkdir(os.path.join(projectFolderPath,"transposonCandD")) 
    os.mkdir(os.path.join(projectFolderPath,"transposonCandE")) 
    os.mkdir(os.path.join(projectFolderPath,"transposonCandF")) 
    os.mkdir(os.path.join(projectFolderPath,"finalResults")) 
    # Copy DNA into folder
    copySequenceClean(inputFasta,projectFolderPath)

#createProject("projects", "testProject", "G:/CambridgeGenData/GenSeq/RHIZIPHAGUS_IRR/rir17contigs.fasta")

from __future__ import print_function
############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
import sys
import os
import os.path

from ProjectCreator import createProject
from AnnotationCommander import runAnnotation
from AnnotationChecker import checkAnnotations
from AnnotationParser import parseAvailableResults
from Statistics import createParsedAnnotationStatistics
from Statistics import createFinalStatistics
from DuplicateFilterA import doFiltering
from DuplicateFilterB import doFilteringB
from SequenceRenamer import renameSequences
from TransposonClustering import doTransposonClustering
from TransposonClusterBlasting import doClusterBlasting
from TransposonClusterBlastAnalysis import doAnalysis
from FinalResultsCreator import createToolAnnotation_Files, createPipelineAnnotation_Files, createFinalAnnotation_Files
from TransposonAnnotator_Help import helpExplanations

# Methods
def getArgument(args, title):
    for i in range(0, len(args)):
        if(args[i].startswith("-"+title)):
            if(i<len(args)-1):
                return args[i+1]
            else:
                return ""
    return ""

#  Main Code
args = sys.argv
mode = getArgument(args, "mode")

if("-help" in args or "-h" in args):
    helpExplanations()
elif(mode=="createProject"):
    arg1 = getArgument(args,"projectFolder")
    arg2 = getArgument(args,"projectName")
    arg3 = getArgument(args,"inputFasta")
    
    error = False
    if(not os.path.isdir(arg1) or arg1==""):
        print("ERROR: projectFolder ",arg1," doesnt exist!")
        error = True
    if(not os.path.isfile(arg3) or arg3==""):
        print("ERROR: inputFasta ",arg3,"could not be found!")
        error = True
    if(os.path.isdir(os.path.join(arg1,arg2)) or arg2==""):
        print("ERROR: project for projectName ",arg2," already exists!")
        error = True
    if(not error):
        createProject(arg1, arg2, arg3)
        print("Project ",arg2," was created successfully in ",arg1)
    
elif(mode=="annotate"):
    arg1 = getArgument(args,"projectFolder")
    arg2 = getArgument(args,"projectName")
    arg3 = getArgument(args,"tool")
    arg4 = ""
    inputLine = " ".join(args)
    addComments = inputLine.split("xxxxx")
    if(len(addComments)>1):
        arg4 = addComments[1]

    validSoftwares = ["helitronScanner", "ltrHarvest", "mitefind", "mitetracker", "must", "repeatmodel", "repMasker", "sinefind", "sinescan", "tirvish", "transposonPSI", "NCBICDD1000", "all"]
    error = False
    if(not os.path.isdir(arg1) or arg1==""):
        print("ERROR: projectFolder ",arg1," doesnt exist!")
        error = True
    if(not os.path.isdir(os.path.join(arg1,arg2))):
        print("ERROR: project for projectName ",arg2," does not exist!")
        error = True
    if(arg3 not in validSoftwares):
        print("ERROR: selected tool ",arg3," does not exist...")
        print("Choose a tool from the following list of allowed tools ")
        print(str(validSoftwares))
        error=True
    if(not error):
        runAnnotation(arg1, arg2, arg3, arg4)
        print("Annotation by software ",arg3," finished successfully...")
        
elif(mode=="checkAnnotations"):
    arg1 = getArgument(args,"projectFolder")
    arg2 = getArgument(args,"projectName")
    error = False
    if(not os.path.isdir(arg1) or arg1==""):
        print("ERROR: projectFolder ",arg1," doesnt exist!")
        error = True
    if(not os.path.isdir(os.path.join(arg1,arg2))):
        print("ERROR: project for projectName ",arg2," does not exist!")
        error = True
    if(not error):
        checkAnnotations(arg1, arg2)
    
elif(mode=="parseAnnotations"):
    arg1 = getArgument(args,"projectFolder")
    arg2 = getArgument(args,"projectName")
    error = False
    if(not os.path.isdir(arg1) or arg1==""):
        print("ERROR: projectFolder ",arg1," doesnt exist!")
        error = True
    if(not os.path.isdir(os.path.join(arg1,arg2))):
        print("ERROR: project for projectName ",arg2," does not exist!")
        error = True
    if(not error):
        projectFolderPath = os.path.join(arg1,arg2)
        parseAvailableResults(projectFolderPath)
        print("Parsing of annotations by softwares finished successfully...")
        
elif(mode=="checkParsed"):
    arg1 = getArgument(args,"projectFolder")
    arg2 = getArgument(args,"projectName")
    error = False
    if(not os.path.isdir(arg1) or arg1==""):
        print("ERROR: projectFolder ",arg1," doesnt exist!")
        error = True
    if(not os.path.isdir(os.path.join(arg1,arg2))):
        print("ERROR: project for projectName ",arg2," does not exist!")
        error = True
    if(not error):
         createParsedAnnotationStatistics(os.path.join(arg1,arg2,"parsedAnnotations"))
         
elif(mode=="statistics"):
    arg1 = getArgument(args,"projectFolder")
    arg2 = getArgument(args,"projectName")
    error = False
    if(not os.path.isdir(arg1) or arg1==""):
        print("ERROR: projectFolder ",arg1," doesnt exist!")
        error = True
    if(not os.path.isdir(os.path.join(arg1,arg2))):
        print("ERROR: project for projectName ",arg2," does not exist!")
        error = True
    if(not error):
         createFinalStatistics(os.path.join(arg1,arg2,"sequence_heads.txt"),os.path.join(arg1,arg2,"finalResults","ToolAnnotations_Transposons.gff3"),os.path.join(arg1,arg2,"finalResults","ToolAnnotations_TransposonMask.gff3"),os.path.join(arg1,arg2,"Statistics_ToolAnnotations.txt"))
         createFinalStatistics(os.path.join(arg1,arg2,"sequence_heads.txt"),os.path.join(arg1,arg2,"finalResults","FinalAnnotations_Transposons.gff3"),os.path.join(arg1,arg2,"finalResults","FinalAnnotations_TransposonMask.gff3"),os.path.join(arg1,arg2,"Statistics_FinalAnnotations.txt"))

elif(mode=="pipeline"):
    arg1 = getArgument(args,"projectFolder")
    arg2 = getArgument(args,"projectName")
    error = False
    if(not os.path.isdir(arg1) or arg1==""):
        print("ERROR: projectFolder ",arg1," doesnt exist!")
        error = True
    if(not os.path.isdir(os.path.join(arg1,arg2))):
        print("ERROR: project for projectName ",arg2," does not exist!")
        error = True
    if(not error):
        # DuplicateFilterA.py
        doFiltering(os.path.join(arg1,arg2,"parsedAnnotations"), os.path.join(arg1,arg2,"transposonCandA"))
        # DuplicateFilterB.py
        doFilteringB(os.path.join(arg1,arg2,"transposonCandA"), os.path.join(arg1,arg2,"transposonCandB"))
        # TransposonClustering.py
        doTransposonClustering(os.path.join(arg1,arg2,"transposonCandB"),os.path.join(arg1,arg2,"transposonCandC"),os.path.join(arg1,arg2,"sequence.fasta"))
        # TransposonClusterBlasting.py
        doClusterBlasting(os.path.join(arg1,arg2,"transposonCandC"), os.path.join(arg1,arg2,"transposonCandD"), os.path.join(arg1,arg2,"sequence.fasta"))
        # TransposonClusterBlastAnalysis.py
        doAnalysis(os.path.join(arg1,arg2,"sequence.fasta"),os.path.join(arg1,arg2,"transposonCandB"),os.path.join(arg1,arg2,"transposonCandC"),os.path.join(arg1,arg2,"transposonCandD"),os.path.join(arg1,arg2,"transposonCandE"))
        # FinalResultsCreator.py
        createToolAnnotation_Files(os.path.join(arg1,arg2), os.path.join(arg1,arg2,"finalResults"), os.path.join(arg1,arg2,"parsedAnnotations"), os.path.join(arg1,arg2,"transposonCandB"), os.path.join(arg1,arg2,"transposonCandF"))
        createPipelineAnnotation_Files(os.path.join(arg1,arg2), os.path.join(arg1,arg2,"finalResults"), os.path.join(arg1,arg2,"parsedAnnotations"), os.path.join(arg1,arg2,"transposonCandE"), os.path.join(arg1,arg2,"transposonCandF"))
        createFinalAnnotation_Files(os.path.join(arg1,arg2,"finalResults"), os.path.join(arg1,arg2,"parsedAnnotations"))
        # Finished
        print("Running pipeline successfully...")

elif(mode=="convertGFF3names"):
    arg1 = getArgument(args,"seqNames")
    arg2 = getArgument(args,"inputGFF")
    arg3 = getArgument(args,"outputGFF")
    renameSequences(arg1, arg2, arg3)
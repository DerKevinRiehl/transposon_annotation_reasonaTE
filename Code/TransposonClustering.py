############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import os
import shutil

# Methods
def getQuery(records, queryChrom, start, end, seq):
    return SeqFeature(FeatureLocation(start, end), type="CDS").extract(seq)

def findOptimal_FactorForClusteringBatchSize(start, end, total):
    limit = start*1.0
    lastND = 1
    lastNU = total
    lastTimes = 0
    n = lastND+(lastNU-lastND)/2
    while True:
        F = (end/start)**(1.0/n)
        times = 0
        target = 0
        counter = 0
        last = 1
        counter = 1
        for i in range(0,total):
            counter += 1
            if(counter>limit):
                times += 1
        #        print(i-last,"\t",i,"\t",limit)
                target = i-last
                last = i
                counter = 0
                limit = limit * F
                if(limit < end):
                    limit = end*1.0
    #    print(target)
        if(target-2>end):
            lastNU = n
        else:
            if(lastTimes<times):
                lastND = n
                lastTimes = times
            else:
                break
        n = int(lastND+(lastNU-lastND)/2)
#        print(n,"\t",lastND,"\t",lastNU,"\t",target,"\t",times)
        if(lastNU-lastND<=2):
            break
    return n, F

def createTransposonSequences(folderIn, folderOut, fastaFile):
    MAXLEN=100000
    # Determine latest iteration
    files = os.listdir(folderIn)
    filesN = list()
    for f in files:
        if(f.startswith("CandidatesB_it")):
            filesN.append(int(f.replace("CandidatesB_it","").replace(".gff3","")))
    filesN.sort()
    fileIn = os.path.join(folderIn,"CandidatesB_it"+str(filesN[-1])+".gff3")
    # Load Annotations, sort by size
    transposons = {}
    counter = 0
    f = open(fileIn,"r")
    line = f.readline()
    while line!="":
        if(not line.startswith("#")):
            parts = line.split("\t")
            if(parts[2]=="transposon"):
                if(parts[0] not in transposons):
                    transposons[parts[0]] = list()
                transposons[parts[0]].append([parts[0],int(parts[3]),int(parts[4]),int(parts[8]),int(parts[4])-int(parts[3])]) 
                counter+=1
        line = f.readline()
    f.close() 
    # Create Fasta File of sequences
    fileOut = os.path.join(folderOut,"TransposonSequences.fasta")
    f = open(fileOut,"w+")
    records = SeqIO.parse(fastaFile, "fasta")
    for r in records:
        for key in list(transposons.keys()):
            if(r.id==key):
                seq = r.seq
                for idx in range(0,len(transposons[key])):
                    query = transposons[key][idx]
                    queryChrom = query[0]
                    queryStart = query[1]
                    queryEnd   = query[2]
                    queryID    = query[3]
                    if(abs(queryEnd-queryStart)<MAXLEN):
                        querySequence   = getQuery(r, queryChrom, queryStart, queryEnd, seq)
                        f.write(">transposon"+str(queryID)+"\n"+str(querySequence)+"\n")
    f.close()
    return counter
    
# Parameters
def doTransposonClustering(folderIn,folderOut,fastaFile,batchSizeStart=200,batchSizeEnd=10):
#    folderIn  = "transposonCandB"
#    folderOut = "transposonCandC"
#    fastaFile = "parsedAnnotations/sequence.fasta"
#    batchSizeStart = 200
#    batchSizeEnd   = 10
#    counter = 21498
    
    # Create Transposon Sequences
    counter = createTransposonSequences(folderIn, folderOut, fastaFile)
    print("Transposon Sequences extracted...")
    
    # Sort Transposon Sequences by size
    os.system("seqkit sort -l "+os.path.join(folderOut,"TransposonSequences.fasta")+" > "+os.path.join(folderOut,"TransposonSequences_sort.fasta"))
    print("Transposon Sequences sorted...")
    
    # Create folder with parts of sequences for clustering
    if(os.path.isdir(os.path.join(folderOut,"clusterIt1"))):
        shutil.rmtree(os.path.join(folderOut,"clusterIt1"))
    os.mkdir(os.path.join(folderOut,"clusterIt1"))
    
    # Start to cluster small portions
    n,F = findOptimal_FactorForClusteringBatchSize(batchSizeStart,batchSizeEnd,counter)
    batchSize = batchSizeStart
    counterZ = 0
    fileCounter = 1
    records = SeqIO.parse(os.path.join(folderOut,"TransposonSequences_sort.fasta"), "fasta")
    selection = list()
    for r in records:
        if(counterZ>batchSize):
            print(counterZ, batchSize)
            SeqIO.write(selection, os.path.join(folderOut,"clusterIt1","SequencesBatch"+str(fileCounter)+".fasta"), "fasta")
            os.system("cd-hit -i "+os.path.join(folderOut,"clusterIt1","SequencesBatch"+str(fileCounter)+".fasta")+" -o "+os.path.join(folderOut,"clusterIt1","Results"+str(fileCounter)+".txt"))
            counterZ = 0
            selection = list()
            fileCounter += 1
            batchSize = batchSize*F
        else:
            counterZ += 1
            selection.append(r)
    print("Clustering finished successfully...")
    SeqIO.write(selection, os.path.join(folderOut,"clusterIt1","SequencesBatch"+str(fileCounter)+".fasta"), "fasta")
    os.system("cd-hit -i "+os.path.join(folderOut,"clusterIt1","SequencesBatch"+str(fileCounter)+".fasta")+" -o "+os.path.join(folderOut,"clusterIt1","Results"+str(fileCounter)+".txt"))
    fileCounter += 1
    counterZ = 0
    selection = list()
    
    # Collect cluster results 
    clusterData = {}
    clusterFileAss = {}
    collected = list()
    previousClusterNo = -1
    previousClusterB = 0
    for i in range(1,fileCounter):
        previousClusterNo = previousClusterB
        f = open(os.path.join(folderOut,"clusterIt1","Results"+str(i)+".txt.clstr"),"r")
        line = f.readline()
        while line!="":
            if(line.startswith(">")):
                if(len(collected)==0):
                    pass
                else:
                    clusterData[previousClusterNo] = collected
                    clusterFileAss[previousClusterNo] = i
                    collected = list()
                previousClusterNo = int(line.replace(">Cluster ","").replace("\n","").replace("",""))+previousClusterB
            else:
                collected.append(int(line.split(">transposon")[1].split(".")[0]))
            line = f.readline()
        f.close()
        clusterData[previousClusterNo] = collected
        clusterFileAss[previousClusterNo] = i
        collected = list()
        previousClusterB = previousClusterNo+1
    print("Collecting clustering results finished...")
    
    # Sort Transposon Sequences by size
    os.system("seqkit sort -n "+os.path.join(folderOut,"TransposonSequences.fasta")+" > "+os.path.join(folderOut,"TransposonSequences_sort2.fasta"))
    print("Transposon Sequences sorted...")
    
    # Create a fasta file and save cluster sequences
    transpIDs = list()
    for key in list(clusterData.keys()):
        cluster = clusterData[key]
        transpIDs.append(["transposon"+str(cluster[0]),key])
    transpIDs.sort()
    fileOut = os.path.join(folderOut,"ClusterSequences.fasta")
    records = SeqIO.parse(os.path.join(folderOut,"TransposonSequences_sort2.fasta"), "fasta")
    f = open(fileOut,"w+")
    for t in transpIDs:
        print(t)
        found = False
        for r in records:
            if(r.id==t[0]):
                seq = r.seq
                f.write(">cluster"+str(t[1])+"\n"+str(seq)+"\n")
                found = True
                break
        if(not found):
            print("ERROR",t)
    f.close()
    
    # Create a Table with cluster and transposon assignment
    fileOut = os.path.join(folderOut,"ClusterSequencesTransposons.txt")
    f = open(fileOut,"w+")
    for key in clusterData:
        f.write(">cluster"+str(key)+"\n"+str(clusterData[key])+"\n")
    f.close()
    
    ## Create a Fasta with a sequence for each cluster at begin and end of file to check if there is potential in clustering these together
    #fileCluAss = {}
    #for i in range(1,fileCounter):
    #    fileCluAss[i] = list()
    #    for ass in clusterFileAss:
    #        if(clusterFileAss[ass]==i):
    #            fileCluAss[i].append(ass)
    #transpIDs = list()
    #for i in range(1,fileCounter):
    #    if(len(fileCluAss[i])>1):
    #        firstCluster = fileCluAss[i][0]
    #        lastCluster  = fileCluAss[i][-1]
    #        transpIDs.append([firstCluster,clusterData[firstCluster][0]])
    #        transpIDs.append([lastCluster,clusterData[lastCluster][0]])
    #    else:
    #        firstCluster = fileCluAss[i][0]
    #        transpIDs.append([firstCluster,clusterData[firstCluster][0]])
    #transpIDs.sort()
    #fileOut = os.path.join(folderOut,"TransposonSequences_ClusterSeqs.fasta")
    #f = open(fileOut,"w+")
    #for t in transpIDs:
    #    print(t)
    #    records = SeqIO.parse(os.path.join(folderOut,"TransposonSequences.fasta"), "fasta")
    #    for r in records:
    #        if(r.id=="transposon"+str(t[1])):
    #            seq = r.seq
    #            f.write(">cluster"+str(t[0])+"\n"+str(seq)+"\n")
    #f.close()
    
    # Cluster again on them
    #os.system("cd-hit -i "+os.path.join(folderOut,"TransposonSequences_ClusterSeqs.fasta")+" -o "+os.path.join(folderOut,"Results_TransposonSequences_ClusterSeqs.txt"))
            
    ## Collect cluster results  2
    #clusterData2 = list()
    #collected = list()
    #f = open(os.path.join(folderOut,"Results_TransposonSequences_ClusterSeqs.txt.clstr"),"r")
    #line = f.readline()
    #while line!="":
    #    if(line.startswith(">")):
    #        if(len(collected)==0):
    #            pass
    #        else:
    #            clusterData2.append(collected)
    #            collected = list()
    #    else:
    #        collected.append(int(line.split(">cluster")[1].split(".")[0]))
    #    line = f.readline()
    #f.close()
    #clusterData2.append(collected)
    #collected = list()
    #print("Collecting clustering results finished...")
    

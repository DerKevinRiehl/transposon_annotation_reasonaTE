from __future__ import print_function
############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

def renameSequences(seqHeadsTxt, gffIn, gffOut):
    # Load seqnames
    f = open(seqHeadsTxt, "r")
    lines = f.readlines()
    f.close()
    keys = []
    vals = []
    for l in lines:
        parts = l.replace("\n","").replace(" ","\t").replace("|","\t").split("\t")
        keys.append(parts[0])
        vals.append(parts[1])     
    n_keys = len(list(set(keys)))
    n_vals = len(list(set(vals)))
    # Convert files
    if(n_keys != n_vals):
        print("ERROR: there are ",n_keys," unique sequences, but the corresponding second column (tabulator, space and | separated) with the original sequence names contains only",n_vals," values! The renaming could not be done. OPERATION CANCELLED!\n\n")
    else:
        dic = {}
        for i in range(0,n_keys):
            dic[keys[i].replace(">","")] = vals[i].replace(">","")
        
        f1 = open(gffIn, "r")
        f2 = open(gffOut, "w+")
        
        line = f1.readline()
        while line!="":
            parts = line.split("\t")
            if(not line.startswith("#")):
                newLine = dic[parts[0]]+"\t"
                newLine += "\t".join(parts[1:])
                f2.write(newLine)
            else:
                f2.write(newLine)
            line = f1.readline()
        f1.close()
        f2.close()
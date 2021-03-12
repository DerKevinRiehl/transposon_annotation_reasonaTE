############################################################################
##### Transposon Annotator reasonaTE - part of Transposon Ultimate #########
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

def findOptimal_FactorForClusteringBatchSize(start, end, total):
    limit = start
    lastND = 1
    lastNU = total
    lastTimes = 0
    n = lastND+(lastNU-lastND)/2
    while True:
        F = (end/start)**(1/n)
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
                    limit = end
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
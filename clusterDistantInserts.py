import sys

isizediffcutoff, numDistInsertCutoff, delLenCutoff = [int(x) for x in sys.argv[1:]]
diffsumcutoff=2*isizediffcutoff

def overlap(cluster1,cluster2):
    global diffsumcutoff
    badcomps = 0
    allcomps = 0
    for ls1,le1,rs1,re1,isize1,read1 in cluster1:
        for ls2,le2,rs2,re2,isize2,read2 in cluster2:
            diffsum = abs(le1-le2) + abs(rs1-rs2)
            if not (diffsum < diffsumcutoff and abs(abs(isize1)-abs(isize2)) < isizediffcutoff):
                badcomps += 1
            allcomps += 1
    return badcomps/float(allcomps) < 0.25 #completely arbitrary cutoff
    
def merge(cluster1,cluster2):
    newclusterh = {}
    for item in cluster1:
        newclusterh[item] = 1
    for item in cluster2:
        newclusterh[item] = 1
    return newclusterh.keys()

def guessBreakpoints(allReads):
    allReads.sort(lambda x,y: cmp(x[1], y[1]))
    med = len(allReads)/2
    sGuess = allReads[med][1]
    eGuess = allReads[med][2]
    i1, i2 = med+1, med-1
    while 1:
        if i1 < len(allReads):
            if allReads[i1][2] > sGuess and allReads[i1][1] < eGuess:
                if allReads[i1][1] > sGuess:
                    sGuess = allReads[i1][1]
                if allReads[i1][2] < eGuess:
                    eGuess = allReads[i1][2]
        if i2 >= 0:
            if allReads[i2][2] > sGuess and allReads[i2][1] < eGuess:
                if allReads[i2][1] > sGuess:
                    sGuess = allReads[i2][1]
                if allReads[i2][2] < eGuess:
                    eGuess = allReads[i2][2]
        i1+=1
        i2-=1
        if i1 >= len(allReads) and i2 < 0:
            break
    return sGuess, eGuess

def pruneCluster(inserts):
    inserts.sort(lambda x,y: cmp(x[1]+x[2],y[1]+y[2]))
    medInsert = inserts[len(inserts)/2]
    le,rs = medInsert[1:3]
    le = int(le)
    rs = int(rs)
    isize = int(medInsert[4])
    newinserts = []
    readCoords = []
    for currinsert in inserts:
        cls,cle,crs,cre,cisize,cread = currinsert
        diffsum = abs(cle-le) + abs(crs-rs)
        if diffsum < diffsumcutoff and abs(abs(int(cisize))-abs(isize)) < isizediffcutoff:
            newinserts.append(currinsert)
            readCoords.append((cls, cle, crs, cre))
    le, rs = guessBreakpoints(readCoords)
    return le,rs,newinserts

clusterh = {}
reads = {}
line = sys.stdin.readline()
while line:
    line = line.strip().split("\t")
    read,c,ls,le,strand1,rs,re,strand2,svtype,edsum,avg_quals,isize = line
    if not clusterh.has_key(c):
        clusterh[c] = []
    ls,le,rs,re,isize = [int(x) for x in [ls,le,rs,re,isize]]
    if not (ls < le and le < rs and rs < re):
        sys.tderr.write("skipping %s\n" %(line))
    else:
        clusterh[c].append([(ls,le,rs,re,isize,read)])
    line = sys.stdin.readline()

keys = clusterh.keys()
keys.sort()
for c in keys:
    i = 0
    maxi = 0
    startover = 0
    clusters = clusterh[c]
    clusters.sort()
    sys.stderr.write("clustering %s\n" %(c))
    while 1:
        change = 0
        while i < len(clusters):
            #sys.stderr.write("examining cluster %d of %d\n" %(i, len(clusters)))
            ichange = 0
            currcluster = clusters[i]
            overlappers = []
            for j in range(i+1,len(clusters)):
                if overlap(currcluster,clusters[j]):
                    change = 1
                    ichange = 1
                    currcluster = merge(currcluster,clusters[j])
                    overlappers.append(j)
            if not ichange:
                i += 1
            else:
                clusters[i] = currcluster
                for j in overlappers[::-1]:
                    clusters.pop(j)
        if not change:
            break
    clusters.sort()
    for cluster in clusters:
        delS, delE, prunedCluster = pruneCluster(cluster)
        if len(prunedCluster) >= numDistInsertCutoff and rs-le >= delLenCutoff:
            outinserts = []
            for ls,le,rs,re,isize,read in prunedCluster:
                outinserts.append("%s,%s,%s,%s,%s,%s" %(ls,le,rs,re,isize,read))
            print "%s,%s,%s\t%s\t%s" %(c, delS, delE, len(outinserts), "\t".join(outinserts))

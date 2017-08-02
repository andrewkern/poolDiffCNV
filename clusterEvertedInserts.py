import sys,os

isizediffcutoff = int(sys.argv[1])
cutoff = isizediffcutoff*2

def overlappingClusters(clust1,clust2):
    global cutoff
    for c1,left_s1,left_e1,left_strand1,left_mapqual1,right_s1,right_e1,right_strand1,right_mapqual1,read1 in clust1:
        for c2,left_s2,left_e2,left_strand2,left_mapqual2,right_s2,right_e2,right_strand2,right_mapqual2,read2 in clust2:
            if c1 == c2 and abs(left_s1-left_s2) < cutoff and abs(right_e1-right_e2) < cutoff:
                return 1

def merge(clust1,clust2):
    h = {}
    for item in clust1:
        h[item] = 1
    for item in clust2:
        h[item] = 1
    return h.keys()

clusters = []
line = sys.stdin.readline()
while line:
    #IRIS_0016:5:1:15275:1500#ACAGTG chr2L   1962446 1962481 -       37      CCTGCGACATGATAGTTAAATATTGGGTTAGGGCTT    effffafff^d^eg_gfffaffdYfee^ggggYagd    1962469 1962504 +       37      TGGGTTAGGGCTTGTTATTACCATGTGTAAGGGATA      ffffchghh_hghhghgaf^cf[ffWffffffff_f
    line = line.strip()
    line = line.split("\t")
    read,c,left_s,left_e,left_strand,left_mapqual,left_readseq,left_qualseq,right_s,right_e,right_strand,right_mapqual,right_readseq,right_qualseq = line
    left_s,left_e,right_s,right_e = int(left_s),int(left_e),int(right_s),int(right_e)
    if left_strand+right_strand != "-+":
        raise Exception
    if int(left_mapqual) > 0 or int(right_mapqual) > 0:
        clusters.append([(c,left_s,left_e,left_strand,left_mapqual,right_s,right_e,right_strand,right_mapqual,read)])
    line = sys.stdin.readline()

maxi = 0
startover = 0
i = 0
while 1:
    change = 0
    while i < len(clusters):
        ichange = 0
        currcluster = clusters[i]
        overlappers = []
        for j in range(i+1,len(clusters)):
            if overlappingClusters(currcluster,clusters[j]):
                change = 1
                ichange = 1
                currcluster = merge(currcluster,clusters[j])
                overlappers.append(j)
                #print len(clusters)
                #raise Exception
        if not ichange:
            sys.stderr.write("on %s of %s\r" %(i,len(clusters)))
            i += 1
        else:
            clusters[i] = currcluster
            for j in overlappers[::-1]:
                clusters.pop(j)
    if not change:
        break

clusters.sort()
for cluster in clusters:
    outstr = []
    allcoords = []
    for c,left_s,left_e,left_strand,left_mapqual,right_s,right_e,right_strand,right_mapqual,read in cluster:
        allcoords.append(left_s)
        allcoords.append(left_e)
        allcoords.append(right_s)
        allcoords.append(right_e)
        outstr.append("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" %(c,left_s,left_e,left_strand,left_mapqual,right_s,right_e,right_strand,right_mapqual,read))
    s = min(allcoords)
    e = max(allcoords)
    print "%s,%s,%s\t%s" %(c,s,e,len(outstr)) + "\t" + "\t".join(outstr)

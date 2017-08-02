import sys,os

pool1FileName, pool2FileName, normconst1, normconst2, cutoff = sys.argv[1:]

cutoff = int(cutoff)
normconst1, normconst2 = float(normconst1), float(normconst2)

def overlappingClustersSameChr(coordls1,coordls2):
    goodcount = 0
    totalcount = 0
    for left_s1,left_e1,right_s1,right_e1 in coordls1:
        for left_s2,left_e2,right_s2,right_e2, in coordls2:
            d1 = abs(left_e1-left_e2)
            d2 = abs(right_s1-right_s2)
            if d1 < cutoff and d2 < cutoff:
                goodcount += 1
            totalcount += 1
    if goodcount/float(totalcount) > 0:
        return 1
    else:
        return 0

pool1File = open(pool1FileName)
lines = pool1File.readlines()
pool1File.close()

h1 = {}
h2 = {}
freqh = {}
otherfreqh = {}
for line in lines:
    origline = line.strip()
    line = origline.split("\t")
    #chr2L,149,335   5,79,339,413,-409,6169BAAXX:1:22:5566:19672#0   44,118,335,409,366,6169BAAXX:1:97:13572:10361#0 26,100,362,436,411,6169BAAXX:1:88:16481:12107#0 53,127,345,419,-367,6169BAAXX:1:19:1320:5560#0  48,122,351,425,378,6169BAAXX:1:51:12313:3019#0  63,138,344,419,-357,UNC2-RDR300275_0059_FC:4:87:15828:16067#0   68,142,366,440,-373,6169BAAXX:1:29:14968:19184#0        70,144,366,440,-371,6169BAAXX:1:79:7843:2106#0  70,144,366,440,371,6169BAAXX:1:101:14969:18345#0        74,149,365,440,-367,UNC2-RDR300275_0059_FC:4:42:8826:17304#0    71,146,391,466,-396,UNC2-RDR300275_0059_FC:4:41:6077:13833#0
    currh = {}
    c,s,e = line[0].split(",")
    s,e = int(s),int(e)
    coverage = int(line[1])
    for insert in line[2:]:
        #5,79,339,413,-409,6169BAAXX:1:22:5566:19672#0
        s1,e1,s2,e2,isize,readid = insert.split(",")
        s1,e1,s2,e2 = int(s1),int(e1),int(s2),int(e2)
        currh[(s1,e1,s2,e2)] = 1
    readcoords = currh.keys()
    readcoords.sort()
    l = e-s+1
    #if (c == "chrX" and ((s > 8000000 and s < 9000000) or (e > 8000000 and e < 9000000))) or l > 10000:
    #    continue
    if not h1.has_key(c):
        h1[c] = {}
        h2[c] = {}
    h1[c][(s,e)] = (readcoords,origline)
    freqh[(c,s,e)] = (normconst1 * len(readcoords),(s,e))
    otherfreqh[(c,s,e)] = (0,"NA")

def getDists(coordls1,coordls2):
    dls = []
    for left_s1,left_e1,right_s1,right_e1 in coordls1:
        for left_s2,left_e2,right_s2,right_e2, in coordls2:
            d1 = abs(left_e1-left_e2)
            d2 = abs(right_s1-right_s2)
            dls.append((d1,d2))
    dls.sort()
    sum1,sum2 = 0,0
    total = 0
    for item in dls:
        total += 1
        sum1 += item[0]
        sum2 += item[1]
    return sum1/float(total),sum2/float(total)

pool2File = open(pool2FileName)
lines = pool2File.readlines()
pool2File.close()

currc = "adsf"
for line in lines:
    origline = line.strip()
    line = origline.split("\t")
    #chr2L,22786,22861,-,37,22792,22867,+,37,UNC2-RDR300275_0059_FC:4:120:5058:6519#0        chr2L,22786,22860,-,37,22793,22867,+,37,6169BAAXX:1:20:17586:4431#0
    currh = {}
    c,s,e = line[0].split(",")
    sys.stderr.write("%s:%s-%s\r" %(c,s,e))
    if currc != c:
        currc=c
        sys.stderr.write("starting %s----------------\n" %(c))
    s,e = int(s),int(e)
    l = e-s+1
    coverage = int(line[1])
    for insert in line[2:]:
        s1,e1,s2,e2,isize,readid = insert.split(",")
        s1,e1,s2,e2 = int(s1),int(e1),int(s2),int(e2)
        currh[(s1,e1,s2,e2)] = 1
    readcoords = currh.keys()
    readcoords.sort()
    freq = normconst2 * len(readcoords)
    found = 0
    if not h2.has_key(c):
        h1[c] = {}
        h2[c] = {}
    for firsts,firste in h1[c].keys():
        if overlappingClustersSameChr(readcoords,h1[c][(firsts,firste)][0]):
            if otherfreqh[(c,firsts,firste)][0] != 0:
                meandist1,meandist2 = getDists(h1[c][(firsts,firste)][0],h2[c][(firsts,firste)][0])
                meandist1b,meandist2b = getDists(h1[c][(firsts,firste)][0],readcoords)
                if meandist1b+meandist2b < meandist1+meandist2:
                    h2[c][(firsts,firste)] = (readcoords,origline)
                    otherfreqh[(c,firsts,firste)] = (freq,(s,e))
            else:
                otherfreqh[(c,firsts,firste)] = (freq,(s,e))
                h2[c][(firsts,firste)] = (readcoords,origline)
                found = 1
    if not found:
        otherfreqh[(c,s,e)] = (freq,(s,e))
        freqh[(c,s,e)] = (0,"NA")
        h2[c][(s,e)] = (readcoords,origline)

def guessBreakpoints(reads1, reads2):
    allReads = reads1+reads2
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

keys = freqh.keys()
keys.sort()
for key in keys:
    c,s,e = key
    freq1,coords1 = freqh[key]
    if coords1 != "NA":
        s1,e1 = coords1
        coords1 = "%s,%s,%s" %(c,s1,e1)
    if h1[c].has_key((s,e)):
        inserts1 = "|".join(h1[c][(s,e)][1].split("\t"))
        readCoords1 = h1[c][(s,e)][0]
    else:
        inserts1 = "NA"
        readCoords1 = []
    freq2,coords2 = otherfreqh[key]
    if coords2 != "NA":
        s2,e2 = coords2
        coords2 = "%s,%s,%s" %(c,s2,e2)
    if h2[c].has_key((s,e)):
        inserts2 = "|".join(h2[c][(s,e)][1].split("\t"))
        readCoords2 = h2[c][(s,e)][0]
    else:
        inserts2 = "NA"
        readCoords2 = []
    s, e = guessBreakpoints(readCoords1, readCoords2)
    coords = "%s,%s,%s" %(c, s, e)
    #print "\t".join([coords,coords1,str(freq1),inserts1,coords2,str(freq2),inserts2])
    print "\t".join([coords,coords1,str(freq1),coords2,str(freq2)])

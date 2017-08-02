import sys,os
from readDepthCorrections import readDepthCorrections

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
    #chr2L,22786,22861,-,37,22792,22867,+,37,UNC2-RDR300275_0059_FC:4:120:5058:6519#0        chr2L,22786,22860,-,37,22793,22867,+,37,6169BAAXX:1:20:17586:4431#0
    currh = {}
    allcoords = []
    coverage = int(line[1])
    for insert in line[2:]:
        c,s1,e1,strand1,mapqual1,s2,e2,strand2,mapqual2,readid = insert.split(",")
        s1,e1,s2,e2 = int(s1),int(e1),int(s2),int(e2)
        currh[(s1,e1,s2,e2)] = 1
        allcoords.append(s1)
        allcoords.append(e1)
        allcoords.append(s2)
        allcoords.append(e2)
    s = min(allcoords)
    e = max(allcoords)
    l = e-s+1
    readcoords = currh.keys()
    readcoords.sort()
    if not h1.has_key(c):
        h1[c] = {}
        h2[c] = {}
    h1[c][(s,e)] = (readcoords,origline)
    freqh[(c,s,e)] = (normconst1 * len(readcoords),(s,e))
    otherfreqh[(c,s,e)] = (0,"NA")

def getMeanDists(coordls1,coordls2):
    dls = []
    for left_s1,left_e1,right_s1,right_e1 in coordls1:
        for left_s2,left_e2,right_s2,right_e2, in coordls2:
            d1 = abs(left_e1-left_e2)
            d2 = abs(right_s1-right_s2)
            dls.append((d1,d2))
    dsum = 0
    for item in dls:
        dsum += item[0] + item[1]
    return dsum / float(len(dls))

pool2File = open(pool2FileName)
lines = pool2File.readlines()
pool2File.close()

collisions = 0
for line in lines:
    origline = line.strip()
    line = origline.split("\t")
    #chr2L,22786,22861,-,37,22792,22867,+,37,UNC2-RDR300275_0059_FC:4:120:5058:6519#0        chr2L,22786,22860,-,37,22793,22867,+,37,6169BAAXX:1:20:17586:4431#0
    currh = {}
    allcoords = []
    coverage = int(line[1])
    for insert in line[2:]:
        c,s1,e1,strand1,mapqual1,s2,e2,strand2,mapqual2,readid = insert.split(",")
        s1,e1,s2,e2 = int(s1),int(e1),int(s2),int(e2)
        currh[(s1,e1,s2,e2)] = 1
        allcoords.append(s1)
        allcoords.append(e1)
        allcoords.append(s2)
        allcoords.append(e2)
    s = min(allcoords)
    e = max(allcoords)
    l = e-s+1
    readcoords = currh.keys()
    readcoords.sort()
    freq = normconst2 * len(readcoords)
    found = 0
    matchcoords = 0,0
    for firsts,firste in h1[c].keys():
        if overlappingClustersSameChr(readcoords,h1[c][(firsts,firste)][0]):
            #if the current cluster in population 2 already has been assigned a match
            if found:
                collisions += 1
                dists1 = getMeanDists(h1[c][(matchcoords[0],matchcoords[1])][0],readcoords)
                dists2 = getMeanDists(h1[c][(firsts,firste)][0],readcoords)
                #if the previous match is better then do nothing
                if dists2 > dists1:
                    continue
                else:
                    #otherwise, match the current cluster in pop 2 with the current cluster in pop 1
                    otherfreqh[(c,firsts,firste)] = (freq,(s,e))
                    h2[c][(firsts,firste)] = (readcoords,origline)
                    found = 1
                    matchcoords = firsts,firste

                    #record that the old match from population 1 now has no match in population 2
                    otherfreqh[(c,matchcoords[0],matchcoords[1])] = (0,"NA")
                    del h2[c][(matchcoords[0],matchcoords[1])]
                    #if the collision below happens as well, then things will not be dealt with well at all, but this should be rare
            #if the matching cluster in population 1 already has been assigned a match
            if otherfreqh[(c,firsts,firste)][0] != 0:
                collisions += 1
                dists1 = getMeanDists(h1[c][(firsts,firste)][0],h2[c][(firsts,firste)][0])
                dists2 = getMeanDists(h1[c][(firsts,firste)][0],readcoords)
                #if the new match is better than the old match for the current cluster in pop 1
                if dists2 < dists1:
                    #record that the old match from population 2 now has no match in population 1 (and never will!)
                    oldfreq,oldse = otherfreqh[(c,firsts,firste)]
                    olds,olde = oldse
                    otherfreqh[(c,olds,olde)] = (oldfreq, (olds,olde))
                    freqh[(c,olds,olde)] = (0, "NA")
                    h2[c][(olds,olde)] = (readcoords,origline)
    
                    #replace the old match with the new match from population 2
                    otherfreqh[(c,firsts,firste)] = (freq,(s,e))
                    h2[c][(firsts,firste)] = (readcoords,origline)
                    found = 1
                    matchcoords = firsts,firste
            else:
                otherfreqh[(c,firsts,firste)] = (freq,(s,e))
                h2[c][(firsts,firste)] = (readcoords,origline)
                found = 1
                matchcoords = firsts,firste
    if not found:
        otherfreqh[(c,s,e)] = (freq,(s,e))
        freqh[(c,s,e)] = (0,"NA")
        h2[c][(s,e)] = (readcoords,origline)
keys = freqh.keys()
keys.sort()
for key in keys:
    c,s,e = key
    freq1,coords1 = freqh[key]
    allCoords = []
    if coords1 != "NA":
        s1,e1 = coords1
        allCoords += [s1, e1]
        coords1 = "%s,%s,%s" %(c,s1,e1)
    if h1[c].has_key((s,e)):
        inserts1 = "|".join(h1[c][(s,e)][1].split("\t"))
    else:
        inserts1 = "NA"
    freq2,coords2 = otherfreqh[key]
    if coords2 != "NA":
        s2,e2 = coords2
        allCoords += [s2, e2]
        coords2 = "%s,%s,%s" %(c,s2,e2)
    if h2[c].has_key((s,e)):
        inserts2 = "|".join(h2[c][(s,e)][1].split("\t"))
    else:
        inserts2 = "NA"
    coords = "%s,%s,%s" %(c,min(allCoords),max(allCoords))
    #print "\t".join([coords,coords1,str(freq1),inserts1,coords2,str(freq2),inserts2])
    print "\t".join([coords,coords1,str(freq1),coords2,str(freq2)])
sys.stderr.write("%s total collisions\n" %(collisions))

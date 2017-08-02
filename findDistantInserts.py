import sys,os,math

maxcutoff = int(sys.argv[1])

def getReadEndpoint(start,mapinfo):
    currnum = ""
    totallen = 0
    for char in mapinfo:
        if char.isdigit():
            currnum += char
        else:
            currnum = int(currnum)
            if char == "D":
                totallen += currnum
            elif char == "I":
                pass
            elif char == "M":
                totallen += currnum
            elif char in ["S", "H"]:
                pass
            else:
                raise Exception
            currnum = ""
    return start+totallen-1

def bin(n):
    maxn = int(math.log(n) / math.log(2))
    newstr = ""
    for i in range(maxn,-1,-1):
        if n >= 2**i:
            newstr += "1"
            n = n - 2**i
        else:
            newstr += "0"
    return newstr

badh = {}
readh = {}
i = 0
line = sys.stdin.readline()
headers = ["@HD", "@PG", "@RG", "@SQ"]
while line:
    i += 1
    if not line[:3] in headers:
        line = line.strip().split("\t")
        read,flag,c,pos1,mapqual,mapinfo,c2,pos2,isize,reads,quals = line[:11]
        isize = int(isize)
        intFlag = int(flag)
        flag = bin(intFlag)
        #the read is paired, primary, and passes qual checks and is not a PCR/optical duplicate
        if flag[-1] == "1" and intFlag <= 512:
            #make sure both reads are mapped
            if flag[-3] == "0" and flag[-4] == "0":
                pos1 = int(pos1)
                pos2 = int(pos2)
                if c2 == "=" and pos1 < pos2:
                    if flag[-5] == "0":
                        strand = "F"
                    else:
                        strand = "R"
                    if not readh.has_key(read):
                        readh[read] = ["",""]
                    s = pos1
                    e = getReadEndpoint(s,mapinfo)
                    if readh[read][0] != "":
                        del readh[read]
                        badh[read] = 1
                    else:
                        readh[read][0] = (c,s,e,strand)
                elif c2 == "=" and pos2 < pos1:
                    if flag[-5] == "0":
                        strand = "F"
                    else:
                        strand = "R"
                    if not readh.has_key(read):
                        readh[read] = ["",""]
                    s = pos1
                    e = getReadEndpoint(s,mapinfo)
                    if readh[read][1] != "":
                        del readh[read]
                        badh[read] = 1
                    else:
                        readh[read][1] = (c,s,e,strand)
                if readh.has_key(read) and readh[read][0] != "" and readh[read][1] != "":
                    lcoords,rcoords = readh[read]
                    c1,ls,le,strand1 = lcoords
                    c2,rs,re,strand2 = rcoords
                    if c1 != c2:
                        raise Exception
                    if rs < ls:
                        print read,lcoords,rcoords
                        raise Exception
                    span = (rs - le) - 1
                    strands = strand1+strand2
                    del readh[read]
                    if strands == "FR" and abs(isize) > maxcutoff:
                        svtype = "deletion"
                        if not badh.has_key(read):
                            print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(read,c1,ls,le,strand1,rs,re,strand2,svtype,"N/A","N/A",abs(isize))
    line = sys.stdin.readline()

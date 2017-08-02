import sys,os

#IRIS_0016:5:1:8236:2009#ACAGTG  163     chrX    11337514        29      16M3I17M        =       11337654        176     CCAACCGCTCTACTACTTTTTTTTTTTTGCACTACC    hhhhhcfhhhhhfhgaffffhgfhhhcccfRcffff    XT:A:M  NM:i:4  SM:i:29 AM:i:29 XM:i:1  XO:i:1  XG:i:3  MD:Z:14C18

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

line = sys.stdin.readline()
readh = {}
headers = ["@HD", "@PG", "@RG", "@SQ"]
while line:
    if not line[:3] in headers:
        line = line.strip().split("\t")
        #IRIS_0016:5:1:16019:2055#ACAGTG 97      chr2L   4534245 37      76M     =       4535255 1086    GGATGAGTAGATTTTACTATAGTTTCAGTTGTTCTGACGCTCGGAGGAGCAGTAGTCTGTTCTGTCACTTGCCCGC      dff_fc_cfc[`ff]ffcccfffffccfff[beddfffcf]dcc[d[^`WJa`[[c^aa_c``__W``dW[fc]cZ    XT:A:U  NM:i:2  SM:i:37 AM:i:18 X0:i:1    X1:i:0  XM:i:2  XO:i:0  XG:i:0  MD:Z:47T26A1
        read,flag,c,pos1,mapqual,mapinfo,c2,pos2,isize,reads,quals = line[:11]
        mapqual = int(mapqual)
        isize = int(isize)
        intFlag = int(flag)
        flag = bin(intFlag)
        #the read is paired, primary, and passes qual checks and is not a PCR/optical duplicate
        if flag[-1] == "1" and intFlag <= 512:
            #make sure both reads are mapped
            if flag[-3] == "0" and flag[-4] == "0":
                pos1 = int(pos1)
                pos2 = int(pos2)
                if flag[-5] == "0":
                    strand = "+"
                else:
                    strand = "-"
                if flag[-6] == "0":
                    strand2 = "+"
                else:
                    strand2 = "-"
                if c2 == "=" and pos1 < pos2 and strand+strand2 == "-+":
                    if not readh.has_key(read):
                        readh[read] = ["",""]
                    s = pos1
                    e = getReadEndpoint(s,mapinfo)
                    if readh[read][0] != "":
                        print read, flag, intFlag
                        raise Exception
                    readh[read][0] = (c,s,e,strand,mapqual,reads,quals)
                elif c2 == "=" and pos2 < pos1 and strand2+strand == "-+":
                    #same chromosome and this one is on the right
                    if not readh.has_key(read):
                        readh[read] = ["",""]
                    s = pos1
                    e = getReadEndpoint(s,mapinfo)
                    if readh[read][1] != "":
                        raise Exception
                    readh[read][1] = (c,s,e,strand,mapqual,reads,quals)
                if readh.has_key(read) and readh[read][0] != "" and readh[read][1] != "":
                    lcoords,rcoords = readh[read]
                    c1,ls,le,strand1,mapqual1,reads1,quals1 = lcoords
                    c2,rs,re,strand2,mapqual2,reads2,quals2 = rcoords
                    if c1 != c2:
                        raise Exception
                    if rs < ls:
                        print read,lcoords,rcoords
                        raise Exception
                    span = (rs - le) - 1
                    strands = strand1+strand2
                    if strands == "-+":
                        print "\t".join([str(x) for x in [read,c1,ls,le,strand1,mapqual1,reads1,quals1,rs,re,strand2,mapqual2,reads2,quals2]])
                        del readh[read]
    line = sys.stdin.readline()

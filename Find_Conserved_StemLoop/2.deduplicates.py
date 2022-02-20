
import GAP, General, Structure, Visual, Colors, sys, os

####################
##    参数校验
###################

if len(sys.argv)<2:
    print(f"Usage: {sys.argv[0]} 1.initial_stemloop.txt 2.uniq_stemloop.txt")
    print(f"1.initial_stemloop.txt example: /150T/zhangqf/lipan/SMART_SHAPE/stemloop_candidates/mac_pro/1.mac_pro_initial_stemloop.txt")
    print(f"2.uniq_stemloop.txt example: /150T/zhangqf/lipan/SMART_SHAPE/stemloop_candidates/mac_pro/2.mac_pro_uniq_stemloop.txt")
    exit(-1)

mouse_gap = GAP.init("/150T/zhangqf/GenomeAnnotation/Gencode/mm10.genomeCoor.bed")

def read_data(inFile, mouse_gap, outFile):
    OUT = open(outFile, "w")
    
    DATA = []
    cur_data = ""
    cur_genome_pos = []
    for line in open(inFile):
        if line[0] == ">":
            if cur_data:
                DATA.append( (cur_genome_pos, cur_data) )
            cur_data = line
            #print(line)
            #data = line[1:-1].split('\t')
            tid, region = line[1:-1].split('\t')[:2]
            start, end = region.split('-')
            start, end = int(start), int(end)
            chrID, chrPos, Strand = mouse_gap.transCoor2genomeCoor(tid, start)
            cur_genome_pos = [ chrID+Strand, chrPos ]
        else:
            cur_data += line
    
    if line[0] == ">":
        if cur_data:
            DATA.append( (cur_genome_pos, cur_data) )
    
    DATA.sort(key=lambda x: x[0])
    
    OUT.writelines(DATA[0][1])
    i = 1
    while i < len(DATA):
        lastchrID, lastchrPos = DATA[i-1][0]
        chrID, chrPos = DATA[i][0]
        if lastchrID==chrID and abs(lastchrPos-chrPos)<5:
            pass
        else:
            OUT.writelines(DATA[i][1])
        i += 1
    
    OUT.close()

read_data(sys.argv[1], mouse_gap, sys.argv[2])



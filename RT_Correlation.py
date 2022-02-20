
importCommon()

def readRT(inFile):
    D_RT = {}
    N_RT = {}
    D_BD = {}
    N_BD = {}
    lineCount = 0
    for line in open(inFile):
        lineCount += 1
        if lineCount % 10000 == 0:
            print("\tlines: ", lineCount)
        
        data = line.strip().split()
        tid = data[0]
        length = int(data[1])
        D_RT[tid] = [0] * length
        N_RT[tid] = [0] * length
        
        D_BD[tid] = [0] * length
        N_BD[tid] = [0] * length
        for i, rtbd in enumerate(data[2:]):
            if rtbd != 'NULL':
                Array = [int(it) for it in rtbd.split(',')]
                D_RT[tid][i] = Array[0] + Array[2]
                D_BD[tid][i] = Array[1] + Array[3]
                N_RT[tid][i] = Array[4] + Array[6]
                N_BD[tid][i] = Array[5] + Array[7]
        if length - N_RT[tid].count(0) < 50 or length - D_RT[tid].count(0) < 50:
            del D_RT[tid]
            del D_BD[tid]
            del N_RT[tid]
            del N_BD[tid]
    
    return D_RT, D_BD, N_RT, N_BD

def common_RT_pair(D_RT_1, D_BD_1, N_RT_1, N_BD_1, D_RT_5, D_BD_5, N_RT_5, N_BD_5, D_RT_25, D_BD_25, N_RT_25, N_BD_25, D_RT_125, D_BD_125, N_RT_125, N_BD_125):
    import numpy as np
    
    RT_matrix = []
    
    min_cov = 200
    tCount = 0
    for tid in set(N_RT_1)&set(N_RT_5)&set(N_RT_25)&set(N_RT_125):
        tCount += 1
        if tCount % 1000 == 0:
            print("\tCount: ", tCount)
        
        for d1_bd,d5_bd,d25_bd,d125_bd,n1_bd,n5_bd,n25_bd,n125_bd, d1,d5,d25,d125,n1,n5,n25,n125 in zip(D_BD_1[tid], N_BD_1[tid], D_BD_5[tid], N_BD_5[tid], D_BD_25[tid], N_BD_25[tid], D_BD_125[tid], N_BD_125[tid], D_RT_1[tid], D_RT_5[tid], D_RT_25[tid], D_RT_125[tid], N_RT_1[tid], N_RT_5[tid], N_RT_25[tid], N_RT_125[tid]):
            if min(d1_bd, d5_bd, d25_bd, d125_bd, n1_bd, n5_bd, n25_bd, n125_bd) < min_cov:
                continue
            else:
                #RT_matrix.append( (np.log2(d1+0.1),np.log2(d5+0.1),np.log2(d25+0.1),np.log2(d125+0.1),np.log2(n1+0.1),np.log2(n5+0.1),np.log2(n25+0.1),np.log2(n125+0.1)) )
                RT_matrix.append( (d1,d5,d25,d125,n1,n5,n25,n125) )
    
    return RT_matrix

def calc_correlation(sample_df):
    Col = len(sample_df.columns)
    rect = General.init_pd_rect(Col, Col, rowNames=list(sample_df.columns), colNames=list(sample_df.columns), init_value=0)
    
    shape_columns = sample_df.columns
    for idx in range(len(shape_columns)):
        for idy in range(len(shape_columns)):
            coor = scipy.stats.pearsonr( list(sample_df.iloc[:,idx]), list(sample_df.iloc[:,idy]) )[0]
            sys.stdout.writelines("%.3f\t" % (coor, ))
            rect.iloc[idx, idy] = round(coor, 3)
        print("")
    
    return rect

####### 1. Read Files

inFile = "/150T/zhangqf/lipan/SMART_SHAPE/shape_score/STAR_hg38_genome/smart_%s_RT/small_RTBD.out"

D_RT_1, D_BD_1, N_RT_1, N_BD_1 = readRT("/150T/zhangqf/lipan/SMART_SHAPE/shape_score/STAR_hg38_genome/smart_1_RT/old/small_RTBD.out") #readRT(inFile % ('5', ))
D_RT_5, D_BD_5, N_RT_5, N_BD_5 = readRT(inFile % ('5', ))
D_RT_25, D_BD_25, N_RT_25, N_BD_25 = readRT(inFile % ('25', ))
D_RT_125, D_BD_125, N_RT_125, N_BD_125 = readRT(inFile % ('125', ))

####### 2. Get common base SHAPE

RT_matrix = common_RT_pair(D_RT_1, D_BD_1, N_RT_1, N_BD_1, D_RT_5, D_BD_5, N_RT_5, N_BD_5, D_RT_25, D_BD_25, N_RT_25, N_BD_25, D_RT_125, D_BD_125, N_RT_125, N_BD_125); print (len(RT_matrix))

####### 3. Calculate correlation

sample_matrix = random.sample(RT_matrix, 1000000)
sample_df = pd.DataFrame(sample_matrix, columns=['D1', 'D5', 'D25', 'D125', 'N1', 'N5', 'N25', 'N125'])
cor_rect = calc_correlation(sample_df)
sns.heatmap(data=cor_rect, annot=True, cmap=sns.diverging_palette(220, 10, sep=80, n=100))
plt.savefig("figs/Fig2d.pdf")
plt.close()

cor_rect.to_csv("figs/Fig2d_data.csv", sep=",")



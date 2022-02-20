
importCommon()

Parser = GAP.init("/150T/zhangqf/GenomeAnnotation/Gencode/mm10.genomeCoor.bed")
ROOT1="/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-09-27-ko/Processing/3.mapping"
ROOT2="/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-10-21-ko/Processing/3.mapping"
ROOT3="/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/3.mapping"

WT_0h_1=ROOT1 + "/wt-0h-1/ReadsPerGene.out.tab"
WT_0h_2=ROOT1 + "/wt-0h-2/ReadsPerGene.out.tab"
WT_0h_3=ROOT2 + "/wt-0h-3/ReadsPerGene.out.tab"
KO_0h_1=ROOT3 + "/ko2-13_0h_rep1/ReadsPerGene.out.tab"
KO_0h_2=ROOT3 + "/ko2-13_0h_rep2/ReadsPerGene.out.tab"
KO_0h_3=ROOT3 + "/ko2-13_0h_rep3/ReadsPerGene.out.tab"

WT_3h_1=ROOT1 + "/wt-3h-1/ReadsPerGene.out.tab"
WT_3h_2=ROOT2 + "/wt-3h-2/ReadsPerGene.out.tab"
WT_3h_3=ROOT2 + "/wt-3h-3/ReadsPerGene.out.tab"
KO_3h_1=ROOT3 + "/ko2-13_3h_rep1/ReadsPerGene.out.tab"
KO_3h_2=ROOT3 + "/ko2-13_3h_rep2/ReadsPerGene.out.tab"
KO_3h_3=ROOT3 + "/ko2-13_3h_rep3/ReadsPerGene.out.tab"


def read_genecount(inFn):
    table = pd.read_csv(inFn, skiprows=4,sep="\t", header=None)
    table.index = table.iloc[:,0]
    table = table.iloc[:,[1,2,3]]
    table.columns = ['unstranded', '1st-strand', '2nd-strand']
    return table

KO_3h_1 = read_genecount(KO_3h_1); print(KO_3h_1.shape)
KO_3h_2 = read_genecount(KO_3h_2); print(KO_3h_2.shape)
KO_3h_3 = read_genecount(KO_3h_3); print(KO_3h_3.shape)
WT_3h_1 = read_genecount(WT_3h_1); print(WT_3h_1.shape)
WT_3h_2 = read_genecount(WT_3h_2); print(WT_3h_2.shape)
WT_3h_3 = read_genecount(WT_3h_3); print(WT_3h_3.shape)

count_table = pd.concat([KO_3h_1['2nd-strand'], KO_3h_2['2nd-strand'], KO_3h_3['2nd-strand'],
    WT_3h_1['2nd-strand'], WT_3h_2['2nd-strand'], WT_3h_3['2nd-strand']],axis=1)
count_table.columns = ['KO_1','KO_2','KO_3','WT_1','WT_2','WT_3']
count_table = count_table[count_table.sum(axis=1)>0]

outfn = '/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/5.cufflinks/genecount_3h.csv'
count_table.to_csv(outfn)



KO_0h_1 = read_genecount(KO_0h_1); print(KO_0h_1.shape)
KO_0h_2 = read_genecount(KO_0h_2); print(KO_0h_2.shape)
KO_0h_3 = read_genecount(KO_0h_3); print(KO_0h_3.shape)
WT_0h_1 = read_genecount(WT_0h_1); print(WT_0h_1.shape)
WT_0h_2 = read_genecount(WT_0h_2); print(WT_0h_2.shape)
WT_0h_3 = read_genecount(WT_0h_3); print(WT_0h_3.shape)

count_table = pd.concat([KO_0h_1['2nd-strand'], KO_0h_2['2nd-strand'], KO_0h_3['2nd-strand'],
    WT_0h_1['2nd-strand'], WT_0h_2['2nd-strand'], WT_0h_3['2nd-strand']],axis=1)
count_table.columns = ['KO_1','KO_2','KO_3','WT_1','WT_2','WT_3']
count_table = count_table[count_table.sum(axis=1)>0]

outfn = '/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/5.cufflinks/genecount_0h.csv'
count_table.to_csv(outfn)



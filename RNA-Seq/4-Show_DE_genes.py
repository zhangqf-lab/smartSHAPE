importCommon()
Parser = GAP.init("/150T/zhangqf/GenomeAnnotation/Gencode/mm10.genomeCoor.bed")

DeSeq2_fn0 = "/150T/zhangqf/lipan/SMART_SHAPE/RNA-Seq/5.DESeq2/DESeq2_0h.csv"
DeSeq2_fn3 = "/150T/zhangqf/lipan/SMART_SHAPE/RNA-Seq/5.DESeq2/DESeq2_3h.csv"

table0 = pd.read_csv(DeSeq2_fn0,sep=",",index_col=0);table0.head()
table3 = pd.read_csv(DeSeq2_fn3,sep=",",index_col=0);table3.head()

def rename_index(table):
    geneNames = []
    for geneID in table.index:
        geneName = Parser.getGeneParser()[geneID]['gene_name']
        geneNames.append(geneName)
    table.index = geneNames

rename_index(table0)
rename_index(table3)

### log2FoldChange是log(WT/KO)!!!!!!!!!
sorted_table0 = table0.sort_values(by='log2FoldChange')
sorted_table3 = table3.sort_values(by='log2FoldChange')

def print_DE_genes(sorted_table):
    i = 0
    count = 0
    for line in sorted_table.values:
        g = sorted_table.index[i]
        i += 1
        baseMean,log2FoldChange,lfcSE,stat,pvalue,padj = line
        if np.isnan(padj) or padj>0.05: continue
        count += 1
        bc = 'default'
        if g.startswith('Il') or g.startswith('Cd') or g.startswith('Ptgs2'):
            bc = 'yellow'
        if log2FoldChange>0.5: # log(WT/KO)>0.5 也就是下调
            print(Colors.f(g,fc='red',bc=bc),end=" ")
        elif log2FoldChange<-0.5: # log(WT/KO)<-0.5 也就是上调
            print(Colors.f(g,fc='blue',bc=bc),end=" ")

# 红色的基因下调，蓝色的基因上调

print_DE_genes(sorted_table0)
print_DE_genes(sorted_table3)


target_genes = ['Itgav','Shc1','Insig1','Ascc2','Plxna1','Plec',
    'Man1a','Sema4b','Nfkbid','Msl2','Id2','Ier5',
    'Ctdsp2','Atp6v0b','Arpc5','Cd47','Rapgef2',
    'Grk2','Furin','Zeb2','Ctnnb1','Kmt2d','Plekhm2',
    'Brd4','Ccl4','Lfng','Spag9','Etf1','Dusp2','Kdm7a',
    'Pcyt1a','Rel','Pfkfb3','Mafb','Git1','Rnf19b']

def show_gene_DE(table, geneName):
    baseMean,log2FoldChange,lfcSE,stat,pvalue,padj = table.loc[geneName]
    if pvalue<0.05 and log2FoldChange<0: # log(WT/KO)<0 显著上调
        geneNameColor = Colors.f(geneName,bc='green',fc='red')
    else:
        geneNameColor = Colors.f(geneName,bc='default',fc='default')
    print( geneNameColor, "log2FC="+str(log2FoldChange)[:5], "pvalue={:.7f}".format(pvalue), "adj-p={:.7f}".format(padj), sep="\t" )

for geneName in target_genes:
    show_gene_DE(sorted_table3, geneName)




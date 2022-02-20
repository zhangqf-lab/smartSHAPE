

def add_label(x,y,text,pos='up',text_c='#000000',line_c='#000000',x_range=[-1,1],y_up_range=[0,1],y_down_range=[-1,0]):
    assert pos in ('up', 'down')
    n_x = x + np.random.uniform(x_range[0], x_range[1])
    if pos=='up':
        n_y = y + np.random.uniform(y_up_range[0], y_up_range[1])
    else:
        n_y = y + np.random.uniform(y_down_range[0], y_down_range[1])
    plt.plot([x, n_x], [y, n_y], c=line_c, lw=0.3)
    if pos=='down':
        n_y -= 0.1
    plt.text(n_x, n_y, text, fontdict={'color':text_c})


def rename_index(table):
    geneNames = []
    for geneID in table.index:
        geneName = Parser.getGeneParser()[geneID]['gene_name']
        geneNames.append(geneName)
    table.index = geneNames

##########################
### 画散点图
##########################

importCommon()
Parser = GAP.init("/150T/zhangqf/GenomeAnnotation/Gencode/mm10.genomeCoor.bed")
DeSeq2_fn3 = "/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/5.cufflinks/DESeq2_3h.csv"
table3 = pd.read_csv(DeSeq2_fn3,sep=",",index_col=0);table3.head()
table3['log2FoldChange'] = -table3['log2FoldChange']

rename_index(table3)

table3 = table3.sort_values(by=['padj'], ascending=False)

target_genes = ['Itgav','Shc1','Insig1','Ascc2','Plxna1','Plec',
    'Man1a','Sema4b','Nfkbid','Msl2','Id2','Ier5',
    'Ctdsp2','Atp6v0b','Arpc5','Cd47','Rapgef2',
    'Grk2','Furin','Zeb2','Ctnnb1','Kmt2d','Plekhm2',
    'Brd4','Ccl4','Lfng','Spag9','Etf1','Dusp2','Kdm7a', 'Il6',
    'Pcyt1a','Rel','Pfkfb3','Mafb','Git1','Rnf19b','Ptgs2', 'Cd274', 
    'Nfkbiz', 'Tnf', 'Mafk', 'Id1', 'Il2', 'Ox40', 'Il12b', 'Cxcl1', 'Cxcl2', 'Cxcl3',
    'Tm2d3', 'Stat3']

colors = np.array( [Colors.RGB['gray']]*table3.shape[0] )
sizes = np.array( [1]*table3.shape[0] )
colors[table3.padj<0.05] = Colors.RGB['green']
for i in range(table3.shape[0]):
    if table3.index[i] in target_genes:
        #colors[i] = Colors.RGB['orange']
        sizes[i] = 3

plt.figure(figsize=(10,10))
plt.axhline(y=0, xmin=-1, xmax=1, linewidth=0.5, linestyle='--', color=Colors.RGB['blue'])
plt.axhline(y=-1, xmin=-1, xmax=1, linewidth=0.5, linestyle='--', color=Colors.RGB['blue'])
plt.axhline(y=1, xmin=-1, xmax=1, linewidth=0.5, linestyle='--', color=Colors.RGB['blue'])
plt.scatter(x=np.log10(table3.baseMean), y=table3.log2FoldChange, s=sizes, c=colors)
x_range=[-2,1]
y_up_range=[0.5,3]
y_down_range=[-4,-0.5]
for i in range(table3.shape[0]):
    if table3.index[i] in target_genes:
        baseMean,log2FoldChange,lfcSE,stat,pvalue,padj = table3.iloc[i]
        if pvalue>0.01 and padj>=0.05: continue
        if log2FoldChange>0:
            pos = 'up'
        else:
            pos = 'down'
        if padj<0.05:
            text_c = Colors.RGB['red']
            line_c = Colors.RGB['red']
        elif pvalue<0.01:
            text_c = Colors.RGB['green']
            line_c = Colors.RGB['green']
        else:
            text_c='#000000'
            line_c='#000000'
        add_label(np.log10(baseMean), log2FoldChange,table3.index[i],pos=pos,text_c=text_c,line_c=line_c,
            x_range=x_range,y_up_range=y_up_range,y_down_range=y_down_range)

plt.xlabel("Log10 of mean of normalized count")
plt.ylabel("Log2 of FoldChange")
plt.savefig("figs/my_figure.pdf")
plt.close()

 

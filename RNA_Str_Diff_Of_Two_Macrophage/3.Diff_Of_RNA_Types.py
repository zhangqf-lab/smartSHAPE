
########################
### C:\Users\hnsfy\Seafile\美玲的分析\20200425PAS分析,更多结构,结构差异\结构差异
### Scan and find all windows >= defined cutoff (quantile 0.95)
########################

importCommon()
random.seed(1234)
matplotlib.use('Agg')

shape_m = General.load_shape("/150T/zhangqf/lipan/SMART_SHAPE/macrophage/7_calcScore/N_m_smartSHAPE.out")
shape_pro = General.load_shape("/150T/zhangqf/lipan/SMART_SHAPE/macrophage/7_calcScore/N_pro_smartSHAPE.out")
gaper = GAP.init("/150T/zhangqf/GenomeAnnotation/Gencode/mm10.genomeCoor.bed")
sequence = General.load_fasta("/150T/zhangqf/GenomeAnnotation/Gencode/mm10_transcriptome.fa")

###############
### Read difference
###############

class Window:
    def __init__(self, tid, start, end, diff, is_dynamic, genename, geneid, genetype):
        self.tid = tid
        self.start = int(start)
        self.end = int(end)
        self.diff = float(diff)
        self.is_dynamic = True if is_dynamic=='+' else False
        self.genename = genename
        self.geneid = geneid
        self.genetype = genetype
    def __repr__(self):
        return f"{self.tid}:{self.start}-{self.end} {self.diff} {self.is_dynamic}"

transDynamic = {}
dynamic_fn = '/150T/zhangqf/lipan/SMART_SHAPE/macrophage/dynamics/m-pro-diff_windows.txt'
IN = open(dynamic_fn)
for line in IN:
    data = line.strip().split()
    if data[0] not in transDynamic:
        transDynamic[data[0]] = []
    transDynamic[data[0]].append( Window(*data) )


###############
### Check the number of total >= dynamics
### trans_windowNum_list record the number of windows of total/dynamics for each transcript
###############

trans_count_dict = {}
trans_windowNum_list = []
for tid in transDynamic:
    trans_count_dict[tid] = [ len(transDynamic[tid]), len([ d for d in transDynamic[tid] if d.is_dynamic ]) ]
    trans_windowNum_list.append( len(transDynamic[tid]) )

print( np.median(trans_windowNum_list) ) # 48.0
print( np.quantile(trans_windowNum_list, 0.05) ) # 4.0
print( np.quantile(trans_windowNum_list, 0.95) ) # 457.0


###############
### Violinplot: 各类RNA的大差异窗口占比
###############

RNA_type_diff_ratio = {}
RNA_type_count = {}
for tid in transDynamic:
    if trans_count_dict[tid][0] <= 10: continue
    ft = gaper.getTransFeature(tid)
    if ft['gene_type']=='protein_coding':
        total = {}
        dynamics = {}
        for window in transDynamic[tid]:
            if window.start<ft['cds_start']:
                total['5utr'] = total.get('5utr',0)+1
            elif window.start<ft['cds_end']:
                total['cds'] = total.get('cds',0)+1
            else:
                total['3utr'] = total.get('3utr',0)+1
        for window in transDynamic.get(tid):
            if not window.is_dynamic:
                continue
            if window.start<ft['cds_start']:
                dynamics['5utr'] = dynamics.get('5utr',0)+1
            elif window.start<ft['cds_end']:
                dynamics['cds'] = dynamics.get('cds',0)+1
            else:
                dynamics['3utr'] = dynamics.get('3utr',0)+1
        for mtype in total:
            if total[mtype] > 10:
                overal_ratio = dynamics.get(mtype, 0) / total[mtype]
                RNA_type_diff_ratio[mtype] = RNA_type_diff_ratio.get(mtype,[]) + [overal_ratio]
                RNA_type_count[mtype] = RNA_type_count.get(mtype,0) + 1
    gt = Seq.format_gene_type(ft['gene_type'])
    overal_ratio = trans_count_dict[tid][1] / trans_count_dict[tid][0]
    RNA_type_diff_ratio[gt] = RNA_type_diff_ratio.get(gt,[]) + [overal_ratio]
    RNA_type_count[gt] = RNA_type_count.get(gt,0) + 1

## Show meidan
for mtype in RNA_type_diff_ratio:
    print(mtype, np.median(RNA_type_diff_ratio[mtype]))

## Show count
for mtype in RNA_type_count:
    print(mtype, RNA_type_count[mtype])

## Show P-values
rna_sorting = ['mRNA','5utr','cds','3utr','pseudogene','lncRNA','miRNA','snoRNA']
p_values = General.init_pd_rect(len(rna_sorting),len(rna_sorting),rna_sorting,rna_sorting)
for mtype1 in rna_sorting:
    for mtype2 in rna_sorting:
        s, p = scipy.stats.ks_2samp(RNA_type_diff_ratio[mtype1], RNA_type_diff_ratio[mtype2], 'two-sided')
        p_values.loc[mtype1, mtype2] = p

print(p_values)


## Save violinplot figure
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(10, 10))

axs[0].set_title('Big Diff Window Ratio')
axs[0].set_ylabel('Ratio')
axs[0].set_ylim(-0.05, 1)
data = [ RNA_type_diff_ratio[mtype] for mtype in rna_sorting ]
colors = [Colors.RGB['pale_red']] * len(data)
Figures.violinPlot(axs[0], data, rna_sorting,colors=colors)

axs[1].set_title('Data size')
axs[1].set_ylabel('count')
data = [ RNA_type_count[item] for item in rna_sorting ]
axs[1].bar(range(len(data)), data)

fig.tight_layout()
fig.savefig(os.environ['HOME']+"/figs/DynamicWindowRatioDistribution.pdf")
plt.close()

###############
### Violinplot: 各类RNA直接做差值的分布
###############

def average_diff(shape_list1, shape_list2, min_shape_count=20):
    diff = [ abs(float(d1)-float(d2)) for d1,d2 in zip(shape_list1, shape_list2) if d1!='NULL' and d2!='NULL' ]
    if len(diff)>=min_shape_count:
        return np.mean(diff)
    else:
        return -1

RNA_type_diff_ratio = {}
RNA_type_count = {}
for tid in set(shape_m)&set(shape_pro):
    ft = gaper.getTransFeature(tid)
    m = shape_m[tid]
    pro = shape_pro[tid]
    if ft['gene_type']=='protein_coding':
        d1,d2 = ft['cds_start'], ft['cds_end']
        utr5_diff = average_diff( m[:d1], pro[:d1] )
        cds_diff = average_diff( m[d1:d2], pro[d1:d2] )
        utr3_diff = average_diff( m[d2:], pro[d2:] )
        if utr5_diff != -1:
            RNA_type_diff_ratio['5utr'] = RNA_type_diff_ratio.get('5utr', []) + [utr5_diff]
            RNA_type_count['5utr'] = RNA_type_count.get('5utr',0) + 1
        if cds_diff != -1:
            RNA_type_diff_ratio['cds'] = RNA_type_diff_ratio.get('cds', []) + [cds_diff]
            RNA_type_count['cds'] = RNA_type_count.get('cds',0) + 1
        if utr3_diff != -1:
            RNA_type_diff_ratio['3utr'] = RNA_type_diff_ratio.get('3utr', []) + [utr3_diff]
            RNA_type_count['3utr'] = RNA_type_count.get('3utr',0) + 1
    overall_diff = average_diff( m, pro )
    if overall_diff != -1:
        gt = Seq.format_gene_type(ft['gene_type'])
        RNA_type_diff_ratio[gt] = RNA_type_diff_ratio.get(gt, []) + [overall_diff]
        RNA_type_count[gt] = RNA_type_count.get(gt,0) + 1


for mtype in RNA_type_diff_ratio:
    print(mtype, np.median(RNA_type_diff_ratio[mtype]))

for mtype in RNA_type_count:
    print(mtype, RNA_type_count[mtype])


rna_sorting = ['mRNA','5utr','cds','3utr','pseudogene','lncRNA','miRNA','snoRNA']
p_values = General.init_pd_rect(len(rna_sorting),len(rna_sorting),rna_sorting,rna_sorting)
for mtype1 in rna_sorting:
    for mtype2 in rna_sorting:
        s, p = scipy.stats.ks_2samp(RNA_type_diff_ratio[mtype1], RNA_type_diff_ratio[mtype2], 'two-sided')
        p_values.loc[mtype1, mtype2] = p

print(p_values)

fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(10, 10))

axs[0].set_title('Big Diff Window Ratio')
axs[0].set_ylabel('Ratio')
data = [ RNA_type_diff_ratio[mtype] for mtype in rna_sorting ]
colors = [Colors.RGB['pale_red']] * len(data)
Figures.violinPlot(axs[0], data, rna_sorting,colors=colors)

axs[1].set_title('Data size')
axs[1].set_ylabel('count')
data = [ RNA_type_count[item] for item in rna_sorting ]
axs[1].bar(range(len(data)), data)

fig.tight_layout()
fig.savefig(os.environ['HOME']+"/figs/DirectBaseDifferenceDistribution.pdf")
plt.close()





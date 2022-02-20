
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

### Get Dynamics file

def output_strdiff_regions(output_fn, cutoff=0.3168):
    windows_record = []
    
    bar = tqdm(total=len(set(shape_m)&set(shape_pro)), leave=True)
    for tid in set(shape_m)&set(shape_pro):
        bar.update(1)
        ft = gaper.getTransFeature(tid)
        diff = [ abs(float(m)-float(pro)) if m!='NULL' and pro!='NULL' else 'NULL' for m,pro in zip(shape_m[tid],shape_pro[tid]) ]
        ws = 5
        i = 0
        while i<len(diff)-ws:
            if 'NULL' not in diff[i:i+ws]:
                wdiff = np.mean(diff[i:i+ws])
                if wdiff > cutoff:
                    tmp = [tid, i+1, i+ws, round(wdiff,3), '+', ft['gene_name'], ft['gene_id'], ft['gene_type']]
                else:
                    tmp = [tid, i+1, i+ws, round(wdiff,3), '-', ft['gene_name'], ft['gene_id'], ft['gene_type']]
                windows_record.append(tmp)
            i += ws    
    bar.close()
    
    OUT = open(output_fn, 'w')
    for record in windows_record:
        print(record[0], record[1], record[2], record[3], record[4], record[5], record[6], record[7], sep="\t", file=OUT)
    OUT.close()

output_fn = '/150T/zhangqf/lipan/SMART_SHAPE/macrophage/dynamics/m-pro-diff_windows.txt'
output_strdiff_regions(output_fn, cutoff=0.3168)

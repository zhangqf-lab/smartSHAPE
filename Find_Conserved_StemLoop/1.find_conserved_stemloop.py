#-*- coding:utf-8 -*-


####################
##    中间的loop区域时 嘧啶-嘌呤-嘧啶 的stemloop
###################

from library import *
from multiprocessing import Pool, TimeoutError

####################
##    参数校验
###################

if len(sys.argv)<3:
    print(f"Usage: {sys.argv[0]} shape_file outputfile")
    print(f"shape_file example: /150T/zhangqf/lipan/SMART_SHAPE/macrophage/7_calcScore/N_m_smartSHAPE.out")
    print(f"outputfile example: /150T/zhangqf/lipan/SMART_SHAPE/stemloop_candidates/1.initial_stemloop.txt")
    exit(-1)

#### Load SHAPE scores

# /150T/zhangqf/lipan/SMART_SHAPE/macrophage/7_calcScore/N_m_smartSHAPE.out
mouse_smartSHAPE = General.load_shape(sys.argv[1])

if os.path.exists(sys.argv[2]):
    print(f"Error: {sys.argv[2]} exists, please delete it at first")
    exit(-1)

# /150T/zhangqf/lipan/SMART_SHAPE/stemloop_candidates/1.initial_stemloop.txt
OUT = open(sys.argv[2], 'w')
ERR = open(sys.argv[2]+".err", 'w')
LOCK = threading.Lock()

######################
## 准备要查找的转录本列表
######################

tids_to_run = []
for tid in mouse_smartSHAPE:
    ft = Homo_Parser['mouse'].getTransFeature(tid)
    if ft['gene_type'] in ('mRNA', 'protein_coding'):
        cds_start = ft['cds_start']
        cds_end = ft['cds_end']
        trans_len = ft['trans_len']
        valid_shape_num_3 = trans_len - mouse_smartSHAPE[tid].count('NULL')
        valid_shape_num_5 = mouse_smartSHAPE[tid][:cds_start].count('NULL')
        if valid_shape_num_5 > 100 or valid_shape_num_3>150:
            tids_to_run.append(tid)

min_structure_score=0.65
min_identity=0.7
ext=3
find_region=['5', '3']




params = []
for i, tid in enumerate(tids_to_run):
    args = (tids_to_run[i], Homo_Parser, geneDescription, mouse_smartSHAPE, homogenes, 'mouse', 
                find_region, OUT, ERR, min_structure_score, min_identity, LOCK, ext)
    params.append( args )

cores = 7
pool = Pool(processes=cores)
pool.map( find_stable_conserved_stemloop_for_trans,  params)




importCommon()

SHAPE_Base = "/150T/zhangqf/lipan/SMART_SHAPE/shape_score/human_rRNA_RT/"
SHAPE_1 = SHAPE_Base+"1/high_shape.out"
SHAPE_5 = SHAPE_Base+"5/high_shape.out"
SHAPE_25 = SHAPE_Base+"25/high_shape.out"
SHAPE_125 = SHAPE_Base+"125/high_shape.out"
SHAPE_wc = '/150T/zhangqf/lipan/SMART_SHAPE/shape_score/human_rRNA_standard/paris_wc/high_shape.out'

shape = {}
shape['1ng'] = General.load_shape(SHAPE_1)
shape['5ng'] = General.load_shape(SHAPE_5)
shape['25ng'] = General.load_shape(SHAPE_25)
shape['125ng'] = General.load_shape(SHAPE_125)
shape['wc'] = General.load_shape(SHAPE_wc)

#sl_cy_18S = "/150T/zhangqf/lipan/SMART_SHAPE/rRNA_accessibility/sunlei_rRNA/sunlei_18S.shape"
#sl_cy_28S = "/150T/zhangqf/lipan/SMART_SHAPE/rRNA_accessibility/sunlei_rRNA/sunlei_28S.shape"
#shape['wc'] = { '18S':General.load_shape(sl_cy_18S)['human_SSU'], '28S':General.load_shape(sl_cy_28S)['human_LSU'] }


Seq = General.load_fasta("/150T/zhangqf/lipan/SMART_SHAPE/ref_sequence/human_rRNA.fa")

dot_18S = General.load_dot("/150T/zhangqf/lipan/SMART_SHAPE/standard_ss/human_18S.dot")['human_18S'][1]
dot_28S = General.load_dot("/150T/zhangqf/lipan/SMART_SHAPE/standard_ss/human_28S.dot")['human_28S'][1]

root = "/150T/zhangqf/lipan/SMART_SHAPE/rRNA_accessibility/"
access_18S = [ line.strip().split()[2] for line in open(root+"access_18S.txt") ]
access_seq_18S = "".join([line.strip().split()[1] for line in open(root+"access_18S.txt")])
access_28S = [ line.strip().split()[2] for line in open(root+"access_28S.txt") ]
access_seq_28S = "".join([line.strip().split()[1] for line in open(root+"access_28S.txt")])

assert access_seq_18S.replace('T','U') == Seq['18S'].replace('T','U')
assert access_seq_28S.replace('T','U') == Seq['28S'].replace('T','U')

assert len(access_seq_18S)==len(shape['wc']['18S'])==len(shape['1ng']['18S'])
assert len(access_seq_28S)==len(shape['wc']['28S'])==len(shape['1ng']['28S'])

def filter_shape(shape_list, dot_str, access_list, min_access=5):
    myshape = []
    mydot = ""
    for s,d,ac in zip(shape_list, dot_str, access_list):
        if ac!='NULL' and float(ac)<min_access:
            continue
        myshape.append(s)
        mydot += d
    return mydot, myshape

########################################
#### 各个不同剂起始量的rRNA的AUC
########################################

## 18S

plt.figure(figsize=(5.5,5))

mydot, myshape = filter_shape(shape['1ng']['18S'], dot_18S, access_18S, min_access=3)
ROC = General.calc_shape_structure_ROC(mydot, myshape);
auc = General.calc_AUC_v2(mydot, myshape)
plt.plot([i[0] for i in ROC], [i[1] for i in ROC], '-', color=Colors.RGB['green'], label=f"1ng AUC={auc:.3}")

mydot, myshape = filter_shape(shape['5ng']['18S'], dot_18S, access_18S, min_access=3)
ROC = General.calc_shape_structure_ROC(mydot, myshape)
auc = General.calc_AUC_v2(mydot, myshape)
plt.plot([i[0] for i in ROC], [i[1] for i in ROC], '-', color=Colors.RGB['blue'], label=f"5ng AUC={auc:.3}")

mydot, myshape = filter_shape(shape['25ng']['18S'], dot_18S, access_18S, min_access=3)
ROC = General.calc_shape_structure_ROC(mydot, myshape)
auc = General.calc_AUC_v2(mydot, myshape)
plt.plot([i[0] for i in ROC], [i[1] for i in ROC], '-', color=Colors.RGB['red'], label=f"25ng AUC={auc:.3}")

mydot, myshape = filter_shape(shape['125ng']['18S'], dot_18S, access_18S, min_access=3)
ROC = General.calc_shape_structure_ROC(mydot, myshape)
auc = General.calc_AUC_v2(mydot, myshape)
plt.plot([i[0] for i in ROC], [i[1] for i in ROC], '-', color=Colors.RGB['yellow'], label=f"125ng AUC={auc:.3}")

mydot, myshape = filter_shape(shape['wc']['18S'], dot_18S, access_18S, min_access=3)
ROC = General.calc_shape_structure_ROC(mydot, myshape)
auc = General.calc_AUC_v2(mydot, myshape)
plt.plot([i[0] for i in ROC], [i[1] for i in ROC], '-', color=Colors.RGB['gray'], label=f"icSHAPE AUC={auc:.3}")

plt.plot([0,1], [0, 1], '--', color='black')
plt.ylim(0,1)
plt.xlim(0,1)
plt.legend()
plt.title(f"18S Number={len(myshape)}")
plt.savefig( os.path.join(os.environ['HOME'], "figs/18S_ROC.pdf") )
plt.show()

## 28S

plt.figure(figsize=(5.5,5))

mydot, myshape = filter_shape(shape['1ng']['28S'], dot_28S, access_28S, min_access=3)
ROC = General.calc_shape_structure_ROC(mydot, myshape)
auc = General.calc_AUC_v2(mydot, myshape)
plt.plot([i[0] for i in ROC], [i[1] for i in ROC], '-', color=Colors.RGB['green'], label=f"1ng AUC={auc:.3}")

mydot, myshape = filter_shape(shape['5ng']['28S'], dot_28S, access_28S, min_access=3)
ROC = General.calc_shape_structure_ROC(mydot, myshape)
auc = General.calc_AUC_v2(mydot, myshape)
plt.plot([i[0] for i in ROC], [i[1] for i in ROC], '-', color=Colors.RGB['blue'], label=f"5ng AUC={auc:.3}")

mydot, myshape = filter_shape(shape['25ng']['28S'], dot_28S, access_28S, min_access=3)
ROC = General.calc_shape_structure_ROC(mydot, myshape)
auc = General.calc_AUC_v2(mydot, myshape)
plt.plot([i[0] for i in ROC], [i[1] for i in ROC], '-', color=Colors.RGB['red'], label=f"25ng AUC={auc:.3}")

mydot, myshape = filter_shape(shape['125ng']['28S'], dot_28S, access_28S, min_access=3)
ROC = General.calc_shape_structure_ROC(mydot, myshape)
auc = General.calc_AUC_v2(mydot, myshape)
plt.plot([i[0] for i in ROC], [i[1] for i in ROC], '-', color=Colors.RGB['yellow'], label=f"125ng AUC={auc:.3}")

mydot, myshape = filter_shape(shape['wc']['28S'], dot_28S, access_28S, min_access=3)
ROC = General.calc_shape_structure_ROC(mydot, myshape)
auc = General.calc_AUC_v2(mydot, myshape)
plt.plot([i[0] for i in ROC], [i[1] for i in ROC], '-', color=Colors.RGB['gray'], label=f"icSHAPE AUC={auc:.3}")

plt.plot([0,1], [0, 1], '--', color='black')
plt.ylim(0,1)
plt.xlim(0,1)
plt.legend()
plt.title(f"28S Number={len(myshape)}")
plt.savefig( os.path.join(os.environ['HOME'], "figs/28S_ROC.pdf") )
plt.show()



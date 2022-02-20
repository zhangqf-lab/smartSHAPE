########################
### C:\Users\hnsfy\Seafile\美玲的分析\20200425PAS分析,更多结构,结构差异\结构差异
### Collect Site difference and window difference
### Print the windows difference quantile and site difference quantile
########################

importCommon()
random.seed(1234)
matplotlib.use('Agg')

shape_m = General.load_shape("/150T/zhangqf/lipan/SMART_SHAPE/macrophage/7_calcScore/N_m_smartSHAPE.out")
shape_pro = General.load_shape("/150T/zhangqf/lipan/SMART_SHAPE/macrophage/7_calcScore/N_pro_smartSHAPE.out")

ws = 5
site_diff_overall = []
window_diff_overall = []
bar = tqdm(total=len(set(shape_m)&set(shape_pro)), leave=False)
for tid in set(shape_m)&set(shape_pro):
    bar.update(1)
    m = [ 'NULL' if d=='NULL' else float(d) for d in shape_m[tid] ]
    pro = [ 'NULL' if d=='NULL' else float(d) for d in shape_pro[tid] ]
    diff = [ abs(pro_-m_) if 'NULL' not in (pro_,m_) else 'NULL' for pro_,m_ in zip(pro, m) ]
    site_diff_overall += [ d for d in diff if 'NULL'!=d ]
    i = 0
    while i<len(diff)-ws:
        if 'NULL' not in diff[i:i+ws]:
            window_diff_overall.append( np.mean(diff[i:i+ws]) )
        i += ws

bar.close()

print( f"site mean {np.mean(site_diff_overall)}; site median {np.median(site_diff_overall)}; window mean {np.mean(window_diff_overall)}; window median {np.median(window_diff_overall)}"  )
# site mean 0.15360594486075202; site median 0.10799999999999998; window mean 0.15190220317092026; window median 0.138
# site mean 0.13956885537252736; site median 0.098; window mean 0.1379846041037647; window median 0.1256

print( "We will define a cutoff as=", np.quantile(window_diff_overall, 0.95) )
# We will define a cutoff as= 0.31679999999999997
# We will define a cutoff as= 0.2836

############
## Plot the distribution of site difference and window difference
############

sns.distplot(random.sample(site_diff_overall,1000000), bins=100, color=Colors.RGB['blue'], label='site_difference')
sns.distplot(random.sample(window_diff_overall,1000000), bins=100, color=Colors.RGB['red'], label='window_difference')
plt.grid(True)
plt.legend()
plt.savefig(os.environ['HOME']+"/figs/diff_dist.pdf")
plt.close()



importCommon()
from matplotlib.gridspec import GridSpec

def countSHAPETrans(inFile, Parser, min_ratio=0.8, min_len=50):
    count = 0
    TransTypeCount = {}
    for line in open(inFile):
        data = line.strip().split()
        nullCount = data[3:].count('NULL')
        totalCount = int(data[1])
        valid_count = totalCount - nullCount
        if 1.0*valid_count/totalCount>=min_ratio and totalCount>=min_len:
            count += 1
            try:
                gt = hg38_parser.getTransFeature(data[0])['gene_type'] #seq.anno_Methods.gene_type(hg38_parser.getTransFeature(data[0])['gene_type'])
            except KeyError:
                continue
            gt = Seq.format_gene_type(gt)
            #gt = hg38_parser.getTransFeature(data[0])['gene_type']
            TransTypeCount[gt] = TransTypeCount.get(gt, 0) + 1
    return count, TransTypeCount

def gene_type_pie(gtc):
    total = sum(gtc.values())
    mRNA = gtc['mRNA']
    lncRNA = gtc['lncRNA']
    others = total - mRNA - lncRNA
    plt.pie((1.0*mRNA/total, 1.0*lncRNA/total, 1.0*others/total), explode=(0.05, 0, 0), colors=('red', 'green', 'gray'), labels=('mRNA', 'lncRNA', 'others'), shadow=False, autopct='%.3f%%')



####### Read Parser

hg38_parser = GAP.init("/150T/zhangqf/GenomeAnnotation/Gencode/hg38.genomeCoor.bed")

####### Read Parser

ROOT = "/150T/zhangqf/lipan/SMART_SHAPE/sample_reads/%s/N_%s_smartSHAPE.out"

cover_ratio = 0.8
count_50M_1, gtc_50M_1 = countSHAPETrans(ROOT % ('50M', '1'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_50M_5, gtc_50M_5 = countSHAPETrans(ROOT % ('50M', '5'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_50M_25, gtc_50M_25 = countSHAPETrans(ROOT % ('50M', '25'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_50M_125, gtc_50M_125 = countSHAPETrans(ROOT % ('50M', '125'), hg38_parser, min_ratio=cover_ratio, min_len=50)

count_100M_1, gtc_100M_1 = countSHAPETrans(ROOT % ('100M', '1'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_100M_5, gtc_100M_5 = countSHAPETrans(ROOT % ('100M', '5'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_100M_25, gtc_100M_25 = countSHAPETrans(ROOT % ('100M', '25'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_100M_125, gtc_100M_125 = countSHAPETrans(ROOT % ('100M', '125'), hg38_parser, min_ratio=cover_ratio, min_len=50)

count_150M_1, gtc_150M_1 = countSHAPETrans(ROOT % ('150M', '1'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_150M_5, gtc_150M_5 = countSHAPETrans(ROOT % ('150M', '5'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_150M_25, gtc_150M_25 = countSHAPETrans(ROOT % ('150M', '25'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_150M_125, gtc_150M_125 = countSHAPETrans(ROOT % ('150M', '125'), hg38_parser, min_ratio=cover_ratio, min_len=50)

count_200M_1, gtc_200M_1 = countSHAPETrans(ROOT % ('200M', '1'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_200M_5, gtc_200M_5 = countSHAPETrans(ROOT % ('200M', '5'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_200M_25, gtc_200M_25 = countSHAPETrans(ROOT % ('200M', '25'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_200M_125, gtc_200M_125 = countSHAPETrans(ROOT % ('200M', '125'), hg38_parser, min_ratio=cover_ratio, min_len=50)

count_250M_1, gtc_250M_1 = countSHAPETrans(ROOT % ('250M', '1'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_250M_5, gtc_250M_5 = countSHAPETrans(ROOT % ('250M', '5'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_250M_25, gtc_250M_25 = countSHAPETrans(ROOT % ('250M', '25'), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_250M_125, gtc_250M_125 = countSHAPETrans(ROOT % ('250M', '125'), hg38_parser, min_ratio=cover_ratio, min_len=50)

inROOT = "/150T/zhangqf/lipan/SMART_SHAPE/shape_score/STAR_hg38_genome/smart_%s_RT/smartSHAPE.out"

count_1, gtc_1 = countSHAPETrans(inROOT % ('1', ), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_5, gtc_5 = countSHAPETrans(inROOT % ('5', ), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_25, gtc_25 = countSHAPETrans(inROOT % ('25', ), hg38_parser, min_ratio=cover_ratio, min_len=50)
count_125, gtc_125 = countSHAPETrans(inROOT % ('125', ), hg38_parser, min_ratio=cover_ratio, min_len=50)


data = []
data.append( (count_50M_1, '1 ng', 50) )
data.append( (count_50M_5, '5 ng', 50) )
data.append( (count_50M_25, '25 ng', 50) )
data.append( (count_50M_125, '125 ng', 50) )

data.append( (count_100M_1, '1 ng', 100) )
data.append( (count_100M_5, '5 ng', 100) )
data.append( (count_100M_25, '25 ng', 100) )
data.append( (count_100M_125, '125 ng', 100) )

data.append( (count_150M_1, '1 ng', 150) )
data.append( (count_150M_5, '5 ng', 150) )
data.append( (count_150M_25, '25 ng', 150) )
data.append( (count_150M_125, '125 ng', 150) )

data.append( (count_200M_1, '1 ng', 200) )
data.append( (count_200M_5, '5 ng', 200) )
data.append( (count_200M_25, '25 ng', 200) )
data.append( (count_200M_125, '125 ng', 200) )

data.append( (count_250M_1, '1 ng', 250) )
data.append( (count_250M_5, '5 ng', 250) )
data.append( (count_250M_25, '25 ng', 250) )
data.append( (count_250M_125, '125 ng', 250) )

data.append( (count_1, '1 ng', 663) )
data.append( (count_5, '5 ng', 390) )
data.append( (count_25, '25 ng', 302) )
data.append( (count_125, '125 ng', 277) )


######################
## Sample Line Plot
######################


x1 = (50, 100, 150, 200, 250, 663)
y1 = (count_50M_1, count_100M_1, count_150M_1, count_200M_1, count_250M_1, count_1)
plt.plot(x1, y1, 'o', color='#4C72B0', label='1 ng')
plt.plot(x1, y1, '-', color='#4C72B0')

x2 = (50, 100, 150, 200, 250, 390)
y2 = (count_50M_5, count_100M_5, count_150M_5, count_200M_5, count_250M_5, count_5)
plt.plot(x2, y2, 'o', color='#55A868', label='5 ng')
plt.plot(x2, y2, '-', color='#55A868')

x3 = (50, 100, 150, 200, 250, 302)
y3 = (count_50M_25, count_100M_25, count_150M_25, count_200M_25, count_250M_25, count_25)
plt.plot(x3, y3, 'o', color='#C44E52', label='25 ng')
plt.plot(x3, y3, '-', color='#C44E52')

x4 = (50, 100, 150, 200, 250, 277)
y4 = (count_50M_125, count_100M_125, count_150M_125, count_200M_125, count_250M_125, count_125)
plt.plot(x4, y4, 'o', color='#8172B2', label='125 ng')
plt.plot(x4, y4, '-', color='#8172B2')

plt.legend()
plt.savefig("figs/sample_reads.pdf")
plt.show()

######################
## Pie plot
######################

the_grid = GridSpec(2, 2)
plt.subplot(the_grid[0, 0], aspect=1)
plt.title("1 ng")
gene_type_pie(gtc_1)
plt.subplot(the_grid[0, 1], aspect=1)
plt.title("5 ng")
gene_type_pie(gtc_5)
plt.subplot(the_grid[1, 0], aspect=1)
plt.title("25 ng")
gene_type_pie(gtc_25)
plt.subplot(the_grid[1, 1], aspect=1)
gene_type_pie(gtc_125)
plt.title("125 ng")

plt.savefig("figs/50M_trans_ratio.pdf")
plt.show()








#??????????????????????????????????
#
#     For 1ng,5ng,25ng,125ng of human species
#
#??????????????????????????????????

source /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/load_data.sh

cd /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/1ng-5ng-25ng-125ng

mkdir -p 1.map_rRNA 2.map_smallRNA 3.map_genome 4.sam2tab 5.calc_SHAPE 6.shape 7.countRT small_RNAs
ln -s /150T/zhangqf/GenomeAnnotation/INDEX/bowtie2/human_rRNA_tRNA_mtRNA human_rRNA
falen /150T/zhangqf/GenomeAnnotation/INDEX/bowtie2/human_rRNA_tRNA_mtRNA/human_rRNA_tRNA_mtRNA.fa > human_rRNA_tRNA_mtRNA.len
fetchSmallRNA.py /150T/zhangqf/GenomeAnnotation/Gencode/hg38_transcriptome.fa small_RNAs/smallRNA.fa 250
bowtie2-build small_RNAs/smallRNA.fa small_RNAs/smallRNA
falen small_RNAs/smallRNA.fa > small_RNAs/smallRNA.len

########################################
#### 1. Map to rRNA
########################################

function map_rRNA_filtermode(){
    input_fq=$1
    output_fq=$2
    icSHAPE-pipe cleanFq -i ${input_fq} \
        -o ${output_fq} \
        -x /150T/zhangqf/GenomeAnnotation/INDEX/bowtie2/human_rRNA_tRNA_mtRNA/human_rRNA_tRNA_mtRNA \
        -p 20 \
        --mode Local   
}
function map_rRNA_alignmode(){
    input_fq=$1
    output_sam=$2
    icSHAPE-pipe cleanFq -i ${input_fq} \
        -o /dev/null \
        -x /150T/zhangqf/GenomeAnnotation/INDEX/bowtie2/human_rRNA_tRNA_mtRNA/human_rRNA_tRNA_mtRNA \
        -p 20 \
        --mode EndToEnd \
        --sam ${output_sam} \
        --bowparam "--norc"
}
export -f map_rRNA_filtermode map_rRNA_alignmode

declare a=(D N)
declare n=(1ng 5ng 25ng 125ng wc cy)
declare r=(1 2)
for name in ${a[@]};
do
    for num in ${n[@]};
    do
        for rep in ${r[@]};
        do
            declare var_name=${name}_${num}_r${rep};
            ll ${!var_name}
            bsub -q Z-ZQF -e error -o log -n 20 "map_rRNA_alignmode ${!var_name} ${var_name}.sam"
            bsub -q Z-ZQF -e error -o log -n 20 "map_rRNA_filtermode ${!var_name} ${var_name}.fq"
        done
    done
done

########################################
#### 2. Map to smallRNA
########################################

function map_smallRNAs(){
    input_fq=$1
    output_fq=$2
    output_sam=$3
    icSHAPE-pipe cleanFq -i ${input_fq} \
        -o ${output_fq} \
        -x ../small_RNAs/smallRNA \
        -p 20 \
        --mode EndToEnd \
        --sam ${output_sam}
}
export -f map_smallRNAs

cd /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/1ng-5ng-25ng-125ng/2.map_smallRNA

declare a=(D N)
declare n=(1ng 5ng 25ng 125ng wc cy)
declare r=(1 2)
for name in ${a[@]};
do
    for num in ${n[@]};
    do
        for rep in ${r[@]};
        do
            declare var_name=${name}_${num}_r${rep};
            bsub -q Z-ZQF -e error -o log -n 20 "map_smallRNAs ../1.map_rRNA/${var_name}.fq ${var_name}.fq ${var_name}.sam"
        done
    done
done

########################################
#### 3. Map to genome
########################################

function map_genome()
{
    infn=$1
    outprex=$2
    icSHAPE-pipe mapGenome -i ${infn} \
        -o ${outprex} \
        -x /150T/zhangqf/GenomeAnnotation/INDEX/STAR/hg38_Gencode \
        --maxMMap 10 \
        --maxMisMatch 2 \
        --alignMode EndToEnd \
        --maxReport 1 \
        -p 20 \
        --noWithin
}
export -f map_genome

cd /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/1ng-5ng-25ng-125ng/3.map_genome

declare a=(D N)
declare n=(1ng 5ng 25ng 125ng wc cy)
declare r=(1 2)
for name in ${a[@]};
do
    for num in ${n[@]};
    do
        for rep in ${r[@]};
        do
            declare var_name=${name}_${num}_r${rep};
            bsub -q Z-ZQF -n 20 -e error -o log "map_genome ../2.map_smallRNA/${var_name}.fq ${var_name}"
        done
    done
done

########################################
#### 4. sam2tab
########################################

cd /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/1ng-5ng-25ng-125ng/4.sam2tab 

declare a=(D N)
declare n=(1ng 5ng 25ng 125ng wc cy)
declare r=(1 2)
for name in ${a[@]};
do
    for num in ${n[@]};
    do
        for rep in ${r[@]};
        do
            declare var_name=${name}_${num}_r${rep};
            bsub -q Z-ZQF -n 5 -e error -o log "icSHAPE-pipe sam2tab -in ../1.map_rRNA/${var_name}.sam -out ${var_name}.rRNA.tab"
            bsub -q Z-ZQF -n 5 -e error -o log "icSHAPE-pipe sam2tab -in ../2.map_smallRNA/${var_name}.sam -out ${var_name}.small.tab"
            bsub -q Z-ZQF -n 5 -e error -o log "icSHAPE-pipe sam2tab -in ../3.map_genome/${var_name}.sorted.bam -out ${var_name}.tab"
        done
    done
done

########################################
#### 5. calcSHAPE
########################################

cd /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/1ng-5ng-25ng-125ng/5.calc_SHAPE

function calcSHAPEScore(){
    rep1=$1
    rep2=$2
    out_gTab=$3
    icSHAPE-pipe calcSHAPENoCont \
        -N ${rep1},${rep2} \
        -size /150T/zhangqf/GenomeAnnotation/INDEX/STAR/hg38_Gencode/chrNameLength.txt \
        -ijf /150T/zhangqf/GenomeAnnotation/INDEX/STAR/hg38_Gencode/sjdbList.fromGTF.out.tab \
        -out ${out_gTab} \
        -wsize 200 \
        -wstep 5 \
        -genome /150T/zhangqf/GenomeAnnotation/genome/hg38.fa \
        -bases A,C,T,G
}
function calcSHAPEScore_rRNA(){
    rep1=$1
    rep2=$2
    out_gTab=$3
    icSHAPE-pipe calcSHAPENoCont \
        -N ${rep1},${rep2} \
        -size ../human_rRNA_tRNA_mtRNA.len \
        -out ${out_gTab} \
        -genome /150T/zhangqf/GenomeAnnotation/INDEX/bowtie2/human_rRNA_tRNA_mtRNA/human_rRNA_tRNA_mtRNA.fa \
        -bases A,C,T,G \
        -non-sliding \
        -noparam
}
function calcSHAPEScore_small(){
    rep1=$1
    rep2=$2
    out_gTab=$3
    icSHAPE-pipe calcSHAPENoCont \
        -N ${rep1},${rep2} \
        -size ../small_RNAs/smallRNA.len \
        -out ${out_gTab} \
        -genome ../small_RNAs/smallRNA.fa \
        -bases A,C,T,G \
        -non-sliding \
        -noparam
}
export -f calcSHAPEScore calcSHAPEScore_rRNA calcSHAPEScore_small

bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore ../4.sam2tab/N_1ng_r1.tab ../4.sam2tab/N_1ng_r2.tab 1ng.gTab"
bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore ../4.sam2tab/N_5ng_r1.tab ../4.sam2tab/N_5ng_r2.tab 5ng.gTab"
bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore ../4.sam2tab/N_25ng_r1.tab ../4.sam2tab/N_25ng_r2.tab 25ng.gTab"
bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore ../4.sam2tab/N_125ng_r1.tab ../4.sam2tab/N_125ng_r2.tab 125ng.gTab"

bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore_rRNA ../4.sam2tab/N_1ng_r1.rRNA.tab ../4.sam2tab/N_1ng_r2.rRNA.tab 1ng.rRNA.gTab"
bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore_rRNA ../4.sam2tab/N_5ng_r1.rRNA.tab ../4.sam2tab/N_5ng_r2.rRNA.tab 5ng.rRNA.gTab"
bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore_rRNA ../4.sam2tab/N_25ng_r1.rRNA.tab ../4.sam2tab/N_25ng_r2.rRNA.tab 25ng.rRNA.gTab"
bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore_rRNA ../4.sam2tab/N_125ng_r1.rRNA.tab ../4.sam2tab/N_125ng_r2.rRNA.tab 125ng.rRNA.gTab"

bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore_small ../4.sam2tab/N_1ng_r1.small.tab ../4.sam2tab/N_1ng_r2.small.tab 1ng.small.gTab"
bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore_small ../4.sam2tab/N_5ng_r1.small.tab ../4.sam2tab/N_5ng_r2.small.tab 5ng.small.gTab"
bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore_small ../4.sam2tab/N_25ng_r1.small.tab ../4.sam2tab/N_25ng_r2.small.tab 25ng.small.gTab"
bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore_small ../4.sam2tab/N_125ng_r1.small.tab ../4.sam2tab/N_125ng_r2.small.tab 125ng.small.gTab"

########################################
#### 6. genSHAPE2TransSHAPE
########################################

function genSHAPE2TransSHAPE(){
    in_gTab=$1
    out_shape=$2
    icSHAPE-pipe genSHAPEToTransSHAPE \
        -i $in_gTab \
        -g /150T/zhangqf/GenomeAnnotation/Gencode/hg38.genomeCoor.bed \
        -p 20 \
        -c 100 \
        -T 0.5 \
        -n 30 \
        -m 0.2 \
        -o $out_shape
}
function genSHAPE2TransSHAPE_rRNA(){
    in_gTab=$1
    out_shape=$2
    icSHAPE-pipe genSHAPEToTransSHAPE \
        -i $in_gTab \
        -s ../human_rRNA_tRNA_mtRNA.len \
        -p 1 \
        -c 100 \
        -T 2 \
        -n 30 \
        -m 0.2 \
        -o $out_shape \
        --app
}
function genSHAPE2TransSHAPE_small(){
    in_gTab=$1
    out_shape=$2
    icSHAPE-pipe genSHAPEToTransSHAPE \
        -i $in_gTab \
        -s ../small_RNAs/smallRNA.len \
        -p 1 \
        -c 100 \
        -T 2 \
        -n 30 \
        -m 0.6 \
        -o $out_shape \
        --app
}
export -f genSHAPE2TransSHAPE genSHAPE2TransSHAPE_rRNA genSHAPE2TransSHAPE_small

cd /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/1ng-5ng-25ng-125ng/6.shape

bsub -q Z-ZQF -n 20 -e error -o log "genSHAPE2TransSHAPE ../5.calc_SHAPE/1ng.gTab 1ng.shape"
bsub -q Z-ZQF -n 20 -e error -o log "genSHAPE2TransSHAPE ../5.calc_SHAPE/5ng.gTab 5ng.shape"
bsub -q Z-ZQF -n 20 -e error -o log "genSHAPE2TransSHAPE ../5.calc_SHAPE/25ng.gTab 25ng.shape"
bsub -q Z-ZQF -n 20 -e error -o log "genSHAPE2TransSHAPE ../5.calc_SHAPE/125ng.gTab 125ng.shape"

bsub -q Z-ZQF -e error -o log "genSHAPE2TransSHAPE_rRNA ../5.calc_SHAPE/1ng.rRNA.gTab 1ng.shape"
bsub -q Z-ZQF -e error -o log "genSHAPE2TransSHAPE_rRNA ../5.calc_SHAPE/5ng.rRNA.gTab 5ng.shape"
bsub -q Z-ZQF -e error -o log "genSHAPE2TransSHAPE_rRNA ../5.calc_SHAPE/25ng.rRNA.gTab 25ng.shape"
bsub -q Z-ZQF -e error -o log "genSHAPE2TransSHAPE_rRNA ../5.calc_SHAPE/125ng.rRNA.gTab 125ng.shape"

bsub -q Z-ZQF -e error -o log "genSHAPE2TransSHAPE_small ../5.calc_SHAPE/1ng.small.gTab 1ng.shape"
bsub -q Z-ZQF -e error -o log "genSHAPE2TransSHAPE_small ../5.calc_SHAPE/5ng.small.gTab 5ng.shape"
bsub -q Z-ZQF -e error -o log "genSHAPE2TransSHAPE_small ../5.calc_SHAPE/25ng.small.gTab 25ng.shape"
bsub -q Z-ZQF -e error -o log "genSHAPE2TransSHAPE_small ../5.calc_SHAPE/125ng.small.gTab 125ng.shape"

########################################
#### 7. CountRT
########################################

cd /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/1ng-5ng-25ng-125ng/7.countRT

bsub -q Z-ZQF -e error -o log \
    "icSHAPE-pipe countRT \
        -in ../4.sam2tab/N_1ng_r1.tab,../4.sam2tab/N_1ng_r2.tab,../4.sam2tab/N_5ng_r1.tab,../4.sam2tab/N_5ng_r2.tab,../4.sam2tab/N_25ng_r1.tab,../4.sam2tab/N_25ng_r2.tab,../4.sam2tab/N_125ng_r1.tab,../4.sam2tab/N_125ng_r2.tab \
        -genome /150T/zhangqf/GenomeAnnotation/genome/hg38.fa \
        -size /150T/zhangqf/GenomeAnnotation/INDEX/STAR/hg38_Gencode/chrNameLength.txt \
        -ijf /150T/zhangqf/GenomeAnnotation/INDEX/STAR/hg38_Gencode/sjdbList.fromGTF.out.tab \
        -omc 50 \
        -out countRT.gTab"

bsub -q Z-ZQF -e error -o log \
    "icSHAPE-pipe genRTBDToTransRTBD \
        -i countRT.gTab \
        -g /150T/zhangqf/GenomeAnnotation/Gencode/hg38.genomeCoor.bed \
        -p 20 \
        --col 5-20 \
        -o countRT.txt"







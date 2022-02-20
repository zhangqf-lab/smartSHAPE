
#??????????????????????????????????
#
#     Calculate smartSHAPE score for macrophage
#
#??????????????????????????????????

source /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/load_data.sh

cd /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/mac_m-mac_pro

mkdir -p 1.map_rRNA 2.map_smallRNA 3.map_genome 4.sam2tab 5.calc_SHAPE 6.shape 7.countRT small_RNAs
ln -s /150T/zhangqf/GenomeAnnotation/INDEX/bowtie2/mouse_rRNA_tRNA_mtRNA mouse_rRNA
falen /150T/zhangqf/GenomeAnnotation/INDEX/bowtie2/mouse_rRNA_tRNA_mtRNA/mouse_rRNA_tRNA_mtRNA.fa > mouse_rRNA_tRNA_mtRNA.len
fetchSmallRNA.py /150T/zhangqf/GenomeAnnotation/Gencode/mm10_transcriptome.fa small_RNAs/smallRNA.fa 250
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
        -x /150T/zhangqf/GenomeAnnotation/INDEX/bowtie2/mouse_rRNA_tRNA_mtRNA/mouse_rRNA_tRNA_mtRNA \
        -p 20 \
        --mode Local   
}
function map_rRNA_alignmode(){
    input_fq=$1
    output_sam=$2
    icSHAPE-pipe cleanFq -i ${input_fq} \
        -o /dev/null \
        -x /150T/zhangqf/GenomeAnnotation/INDEX/bowtie2/mouse_rRNA_tRNA_mtRNA/mouse_rRNA_tRNA_mtRNA \
        -p 20 \
        --mode EndToEnd \
        --sam ${output_sam} \
        --bowparam "--norc"
}
export -f map_rRNA_filtermode map_rRNA_alignmode

bsub -q Z-ZQF -e error -o log -n 20 "map_rRNA_filtermode ${mac_m_r1} mac_m_r1.fq"
bsub -q Z-ZQF -e error -o log -n 20 "map_rRNA_alignmode ${mac_m_r1} mac_m_r1.sam"

bsub -q Z-ZQF -e error -o log -n 20 "map_rRNA_filtermode ${mac_m_r2} mac_m_r2.fq"
bsub -q Z-ZQF -e error -o log -n 20 "map_rRNA_alignmode ${mac_m_r2} mac_m_r2.sam"

bsub -q Z-ZQF -e error -o log -n 20 "map_rRNA_filtermode ${mac_pro_r1} mac_pro_r1.fq"
bsub -q Z-ZQF -e error -o log -n 20 "map_rRNA_alignmode ${mac_pro_r1} mac_pro_r1.sam"

bsub -q Z-ZQF -e error -o log -n 20 "map_rRNA_filtermode ${mac_pro_r2} mac_pro_r2.fq"
bsub -q Z-ZQF -e error -o log -n 20 "map_rRNA_alignmode ${mac_pro_r2} mac_pro_r2.sam"

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

cd /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/mac_m-mac_pro/2.map_smallRNA

bsub -q Z-ZQF -e error -o log -n 20 "map_smallRNAs ../1.map_rRNA/mac_m_r1.fq mac_m_r1.fq mac_m_r1.sam"
bsub -q Z-ZQF -e error -o log -n 20 "map_smallRNAs ../1.map_rRNA/mac_m_r2.fq mac_m_r2.fq mac_m_r2.sam"
bsub -q Z-ZQF -e error -o log -n 20 "map_smallRNAs ../1.map_rRNA/mac_pro_r1.fq mac_pro_r1.fq mac_pro_r1.sam"
bsub -q Z-ZQF -e error -o log -n 20 "map_smallRNAs ../1.map_rRNA/mac_pro_r2.fq mac_pro_r2.fq mac_pro_r2.sam"

########################################
#### 3. Map to genome
########################################

function map_genome()
{
    infn=$1
    outprex=$2
    icSHAPE-pipe mapGenome -i ${infn} \
        -o ${outprex} \
        -x /150T/zhangqf/GenomeAnnotation/INDEX/STAR/mm10_Gencode \
        --maxMMap 10 \
        --maxMisMatch 2 \
        --alignMode EndToEnd \
        --maxReport 1 \
        -p 20 \
        --noWithin
}
export -f map_genome

cd /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/mac_m-mac_pro/3.map_genome

bsub -q Z-ZQF -n 20 -e error -o log "map_genome ../2.map_smallRNA/mac_m_r1.fq mac_m_r1"
bsub -q Z-ZQF -n 20 -e error -o log "map_genome ../2.map_smallRNA/mac_m_r2.fq mac_m_r2"
bsub -q Z-ZQF -n 20 -e error -o log "map_genome ../2.map_smallRNA/mac_pro_r1.fq mac_pro_r1"
bsub -q Z-ZQF -n 20 -e error -o log "map_genome ../2.map_smallRNA/mac_pro_r2.fq mac_pro_r2"

########################################
#### 4. sam2tab
########################################

cd /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/mac_m-mac_pro/4.sam2tab

samples=(mac_m_r1 mac_m_r2 mac_pro_r1 mac_pro_r2)
for sample in ${samples[@]}; do
    bsub -q Z-ZQF -n 5 -e error -o log "icSHAPE-pipe sam2tab -in ../1.map_rRNA/${sample}.sam -out ${sample}.rRNA.tab"
    bsub -q Z-ZQF -n 5 -e error -o log "icSHAPE-pipe sam2tab -in ../2.map_smallRNA/${sample}.sam -out ${sample}.small.tab"
    bsub -q Z-ZQF -n 5 -e error -o log "icSHAPE-pipe sam2tab -in ../3.map_genome/${sample}.sorted.bam -out ${sample}.tab"
done

########################################
#### 5. calcSHAPE
########################################

cd /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/mac_m-mac_pro/5.calc_SHAPE

function calcSHAPEScore(){
    rep1=$1
    rep2=$2
    out_gTab=$3
    icSHAPE-pipe calcSHAPENoCont \
        -N ${rep1},${rep2} \
        -size /150T/zhangqf/GenomeAnnotation/INDEX/STAR/mm10_Gencode/chrNameLength.txt \
        -ijf /150T/zhangqf/GenomeAnnotation/INDEX/STAR/mm10_Gencode/sjdbList.fromGTF.out.tab \
        -out ${out_gTab} \
        -wsize 200 \
        -wstep 5 \
        -genome /150T/zhangqf/GenomeAnnotation/genome/mm10.fa \
        -bases A,C,T,G
}
function calcSHAPEScore_rRNA(){
    rep1=$1
    rep2=$2
    out_gTab=$3
    icSHAPE-pipe calcSHAPENoCont \
        -N ${rep1},${rep2} \
        -size ../mouse_rRNA_tRNA_mtRNA.len \
        -out ${out_gTab} \
        -genome /150T/zhangqf/GenomeAnnotation/INDEX/bowtie2/mouse_rRNA_tRNA_mtRNA/mouse_rRNA_tRNA_mtRNA.fa \
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


bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore ../4.sam2tab/mac_m_r1.tab ../4.sam2tab/mac_m_r2.tab mac_m.gTab"
bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore ../4.sam2tab/mac_pro_r1.tab ../4.sam2tab/mac_pro_r2.tab mac_pro.gTab"

bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore_rRNA ../4.sam2tab/mac_m_r1.rRNA.tab ../4.sam2tab/mac_m_r2.rRNA.tab mac_m.rRNA.gTab"
bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore_rRNA ../4.sam2tab/mac_pro_r1.rRNA.tab ../4.sam2tab/mac_pro_r2.rRNA.tab mac_pro.rRNA.gTab"

bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore_small ../4.sam2tab/mac_m_r1.small.tab ../4.sam2tab/mac_m_r2.small.tab mac_m.small.gTab"
bsub -q Z-ZQF -n 5 -e error -o log "calcSHAPEScore_small ../4.sam2tab/mac_pro_r1.small.tab ../4.sam2tab/mac_pro_r2.small.tab mac_pro.small.gTab"

########################################
#### 6. genSHAPE2TransSHAPE
########################################

function genSHAPE2TransSHAPE(){
    in_gTab=$1
    out_shape=$2
    icSHAPE-pipe genSHAPEToTransSHAPE \
        -i $in_gTab \
        -g /150T/zhangqf/GenomeAnnotation/Gencode/mm10.genomeCoor.bed \
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
        -s ../mouse_rRNA_tRNA_mtRNA.len \
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

cd /150T/zhangqf/lipan/SMART_SHAPE/2020-ReAnalyze/mac_m-mac_pro/6.shape

bsub -q Z-ZQF -n 20 -e error -o log "genSHAPE2TransSHAPE ../5.calc_SHAPE/mac_m.gTab mac_m.shape"
bsub -q Z-ZQF -n 20 -e error -o log "genSHAPE2TransSHAPE ../5.calc_SHAPE/mac_pro.gTab mac_pro.shape"

bsub -q Z-ZQF -n 20 -e error -o log "genSHAPE2TransSHAPE_rRNA ../5.calc_SHAPE/mac_m.rRNA.gTab mac_m.shape"
bsub -q Z-ZQF -n 20 -e error -o log "genSHAPE2TransSHAPE_rRNA ../5.calc_SHAPE/mac_pro.rRNA.gTab mac_pro.shape"

bsub -q Z-ZQF -n 20 -e error -o log "genSHAPE2TransSHAPE_small ../5.calc_SHAPE/mac_m.small.gTab mac_m.shape"
bsub -q Z-ZQF -n 20 -e error -o log "genSHAPE2TransSHAPE_small ../5.calc_SHAPE/mac_pro.small.gTab mac_pro.shape"


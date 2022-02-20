
IN=/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko
OUT=/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing
adaptor=/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/adaptor.fa

mkdir 1.trim 2.removeAbound 3.mapping 4.cuffdiff 5.cufflinks

#######  ko2-13_0h_rep1

bsub -q Z-ZQF -n 20 -e error -o log \
    "java -jar ~/usr/icSHAPE-pipe/bin/Functions/trimmomatic-0.38.jar PE \
        -threads 20 \
        ${IN}/ko2-13_0h_rep1_R1.fastq \
        ${IN}/ko2-13_0h_rep1_R2.fastq \
        ${OUT}/1.trim/ko2-13_0h_rep1_forward_paired.fastq.gz \
        ${OUT}/1.trim/ko2-13_0h_rep1_forward_unpaired.fastq.gz \
        ${OUT}/1.trim/ko2-13_0h_rep1_reverse_paired.fastq.gz \
        ${OUT}/1.trim/ko2-13_0h_rep1_reverse_unpaired.fastq.gz \
        ILLUMINACLIP:${adaptor}:2:30:10:2:keepBothReads \
        LEADING:3 \
        TRAILING:3 \
        MINLEN:36"

#######  ko2-13_0h_rep2

bsub -q Z-ZQF -n 20 -e error -o log \
    "java -jar ~/usr/icSHAPE-pipe/bin/Functions/trimmomatic-0.38.jar PE \
        -threads 20 \
        ${IN}/ko2-13_0h_rep2_R1.fastq \
        ${IN}/ko2-13_0h_rep2_R2.fastq \
        ${OUT}/1.trim/ko2-13_0h_rep2_forward_paired.fastq.gz \
        ${OUT}/1.trim/ko2-13_0h_rep2_forward_unpaired.fastq.gz \
        ${OUT}/1.trim/ko2-13_0h_rep2_reverse_paired.fastq.gz \
        ${OUT}/1.trim/ko2-13_0h_rep2_reverse_unpaired.fastq.gz \
        ILLUMINACLIP:${adaptor}:2:30:10:2:keepBothReads \
        LEADING:3 \
        TRAILING:3 \
        MINLEN:36"

#######  ko2-13_0h_rep3

bsub -q Z-ZQF -n 20 -e error -o log \
    "java -jar ~/usr/icSHAPE-pipe/bin/Functions/trimmomatic-0.38.jar PE \
        -threads 20 \
        ${IN}/ko2-13_0h_rep3_R1.fastq \
        ${IN}/ko2-13_0h_rep3_R2.fastq \
        ${OUT}/1.trim/ko2-13_0h_rep3_forward_paired.fastq.gz \
        ${OUT}/1.trim/ko2-13_0h_rep3_forward_unpaired.fastq.gz \
        ${OUT}/1.trim/ko2-13_0h_rep3_reverse_paired.fastq.gz \
        ${OUT}/1.trim/ko2-13_0h_rep3_reverse_unpaired.fastq.gz \
        ILLUMINACLIP:${adaptor}:2:30:10:2:keepBothReads \
        LEADING:3 \
        TRAILING:3 \
        MINLEN:36"

#######  ko2-13_3h_rep1

bsub -q Z-ZQF -n 20 -e error -o log \
    "java -jar ~/usr/icSHAPE-pipe/bin/Functions/trimmomatic-0.38.jar PE \
        -threads 20 \
        ${IN}/ko2-13_3h_rep1_R1.fastq \
        ${IN}/ko2-13_3h_rep1_R2.fastq \
        ${OUT}/1.trim/ko2-13_3h_rep1_forward_paired.fastq.gz \
        ${OUT}/1.trim/ko2-13_3h_rep1_forward_unpaired.fastq.gz \
        ${OUT}/1.trim/ko2-13_3h_rep1_reverse_paired.fastq.gz \
        ${OUT}/1.trim/ko2-13_3h_rep1_reverse_unpaired.fastq.gz \
        ILLUMINACLIP:${adaptor}:2:30:10:2:keepBothReads \
        LEADING:3 \
        TRAILING:3 \
        MINLEN:36"

#######  ko2-13_3h_rep2

bsub -q Z-ZQF -n 20 -e error -o log \
    "java -jar ~/usr/icSHAPE-pipe/bin/Functions/trimmomatic-0.38.jar PE \
        -threads 20 \
        ${IN}/ko2-13_3h_rep2_R1.fastq \
        ${IN}/ko2-13_3h_rep2_R2.fastq \
        ${OUT}/1.trim/ko2-13_3h_rep2_forward_paired.fastq.gz \
        ${OUT}/1.trim/ko2-13_3h_rep2_forward_unpaired.fastq.gz \
        ${OUT}/1.trim/ko2-13_3h_rep2_reverse_paired.fastq.gz \
        ${OUT}/1.trim/ko2-13_3h_rep2_reverse_unpaired.fastq.gz \
        ILLUMINACLIP:${adaptor}:2:30:10:2:keepBothReads \
        LEADING:3 \
        TRAILING:3 \
        MINLEN:36"

#######  ko2-13_3h_rep3

bsub -q Z-ZQF -n 20 -e error -o log \
    "java -jar ~/usr/icSHAPE-pipe/bin/Functions/trimmomatic-0.38.jar PE \
        -threads 20 \
        ${IN}/ko2-13_3h_rep3_R1.fastq \
        ${IN}/ko2-13_3h_rep3_R2.fastq \
        ${OUT}/1.trim/ko2-13_3h_rep3_forward_paired.fastq.gz \
        ${OUT}/1.trim/ko2-13_3h_rep3_forward_unpaired.fastq.gz \
        ${OUT}/1.trim/ko2-13_3h_rep3_reverse_paired.fastq.gz \
        ${OUT}/1.trim/ko2-13_3h_rep3_reverse_unpaired.fastq.gz \
        ILLUMINACLIP:${adaptor}:2:30:10:2:keepBothReads \
        LEADING:3 \
        TRAILING:3 \
        MINLEN:36"


########################
### Remove Abundant RNAs
########################


function removeAbundant()
{
    if [[ $# -ne 3 ]]; then
        echo "Usage: removeAbundant forward.fq.gz reversed.fq.gz output"
        return
    fi
    forward_reads=$1
    reversed_reads=$2
    OUTPUT=$3
    STAR --runThreadN 20 \
        --genomeDir /150T/zhangqf/GenomeAnnotation/INDEX/STAR/mouse_rRNA_tRNA_mtRNA \
        --readFilesIn ${forward_reads} ${reversed_reads} \
        --readFilesCommand \"gzip -dc\" \
        --outFilterType Normal \
        --outFilterMultimapNmax 999 \
        --alignSJoverhangMin 8 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 21 \
        --alignIntronMax 0 \
        --alignMatesGapMax 0 \
        --outFileNamePrefix ${OUTPUT} \
        --outSAMmultNmax -1 \
        --outSAMattributes All \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM SortedByCoordinate \
        --chimOutType WithinBAM \
        --alignEndsType Local \
        --outReadsUnmapped Fastx
}
export -f removeAbundant

IN=/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/1.trim
OUT=/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/2.removeAbound

bsub -q Z-ZQF -n 20 -e error \
    "removeAbundant ${IN}/ko2-13_0h_rep1_forward_paired.fastq.gz ${IN}/ko2-13_0h_rep1_reverse_paired.fastq.gz ${OUT}/ko2-13_0h_rep1/"

bsub -q Z-ZQF -n 20 -e error \
    "removeAbundant ${IN}/ko2-13_0h_rep2_forward_paired.fastq.gz ${IN}/ko2-13_0h_rep2_reverse_paired.fastq.gz ${OUT}/ko2-13_0h_rep2/"

bsub -q Z-ZQF -n 20 -e error \
    "removeAbundant ${IN}/ko2-13_0h_rep3_forward_paired.fastq.gz ${IN}/ko2-13_0h_rep3_reverse_paired.fastq.gz ${OUT}/ko2-13_0h_rep3/"

bsub -q Z-ZQF -n 20 -e error \
    "removeAbundant ${IN}/ko2-13_3h_rep1_forward_paired.fastq.gz ${IN}/ko2-13_3h_rep1_reverse_paired.fastq.gz ${OUT}/ko2-13_3h_rep1/"

bsub -q Z-ZQF -n 20 -e error \
    "removeAbundant ${IN}/ko2-13_3h_rep2_forward_paired.fastq.gz ${IN}/ko2-13_3h_rep2_reverse_paired.fastq.gz ${OUT}/ko2-13_3h_rep2/"

bsub -q Z-ZQF -n 20 -e error \
    "removeAbundant ${IN}/ko2-13_3h_rep3_forward_paired.fastq.gz ${IN}/ko2-13_3h_rep3_reverse_paired.fastq.gz ${OUT}/ko2-13_3h_rep3/"

########################
### Pair-end mapping
########################

function mapping()
{
    if [[ $# -ne 3 ]]; then
        echo "Usage: mapping forward.fq.gz reversed.fq.gz output"
        return
    fi
    forward_reads=$1
    reversed_reads=$2
    OUTPUT=$3
    STAR --runThreadN 20 \
        --genomeDir /150T/zhangqf/GenomeAnnotation/INDEX/STAR/mm10_Gencode \
        --readFilesIn ${forward_reads} ${reversed_reads} \
        --outFilterType Normal \
        --outFilterMultimapNmax 10 \
        --alignSJoverhangMin 8 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 21 \
        --alignIntronMax 0 \
        --alignMatesGapMax 0 \
        --outFileNamePrefix ${OUTPUT} \
        --outSAMattrIHstart 0 \
        --outSAMmultNmax -1 \
        --outSAMattributes All \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM Unsorted \
        --chimOutType WithinBAM \
        --alignEndsType Local \
        --quantMode TranscriptomeSAM GeneCounts
}

export -f mapping

IN=/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/2.removeAbound
OUT=/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/3.mapping

mkdir -p ${OUT}/ko2-13_0h_rep1 ${OUT}/ko2-13_0h_rep2 ${OUT}/ko2-13_0h_rep3 ${OUT}/ko2-13_3h_rep1 ${OUT}/ko2-13_3h_rep2 ${OUT}/ko2-13_3h_rep3
bsub -q Z-ZQF -n 20 -e error "mapping ${IN}/ko2-13_0h_rep1/Unmapped.out.mate1 ${IN}/ko2-13_0h_rep1/Unmapped.out.mate2 ${OUT}/ko2-13_0h_rep1/"
bsub -q Z-ZQF -n 20 -e error "mapping ${IN}/ko2-13_0h_rep2/Unmapped.out.mate1 ${IN}/ko2-13_0h_rep2/Unmapped.out.mate2 ${OUT}/ko2-13_0h_rep2/"
bsub -q Z-ZQF -n 20 -e error "mapping ${IN}/ko2-13_0h_rep3/Unmapped.out.mate1 ${IN}/ko2-13_0h_rep3/Unmapped.out.mate2 ${OUT}/ko2-13_0h_rep3/"
bsub -q Z-ZQF -n 20 -e error "mapping ${IN}/ko2-13_3h_rep1/Unmapped.out.mate1 ${IN}/ko2-13_3h_rep1/Unmapped.out.mate2 ${OUT}/ko2-13_3h_rep1/"
bsub -q Z-ZQF -n 20 -e error "mapping ${IN}/ko2-13_3h_rep2/Unmapped.out.mate1 ${IN}/ko2-13_3h_rep2/Unmapped.out.mate2 ${OUT}/ko2-13_3h_rep2/"
bsub -q Z-ZQF -n 20 -e error "mapping ${IN}/ko2-13_3h_rep3/Unmapped.out.mate1 ${IN}/ko2-13_3h_rep3/Unmapped.out.mate2 ${OUT}/ko2-13_3h_rep3/"


########################
### sort bam
########################

DIR=/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/3.mapping

bsub -q Z-ZQF -n 20 -e error "samtools sort --threads 20 -m 100G ${DIR}/ko2-13_0h_rep1/Aligned.out.bam > ${DIR}/ko2-13_0h_rep1/Aligned.SortedByCoordinate.bam"
bsub -q Z-ZQF -n 20 -e error "samtools sort --threads 20 -m 100G ${DIR}/ko2-13_0h_rep2/Aligned.out.bam > ${DIR}/ko2-13_0h_rep2/Aligned.SortedByCoordinate.bam"
bsub -q Z-ZQF -n 20 -e error "samtools sort --threads 20 -m 100G ${DIR}/ko2-13_0h_rep3/Aligned.out.bam > ${DIR}/ko2-13_0h_rep3/Aligned.SortedByCoordinate.bam"
bsub -q Z-ZQF -n 20 -e error "samtools sort --threads 20 -m 100G ${DIR}/ko2-13_3h_rep1/Aligned.out.bam > ${DIR}/ko2-13_3h_rep1/Aligned.SortedByCoordinate.bam"
bsub -q Z-ZQF -n 20 -e error "samtools sort --threads 20 -m 100G ${DIR}/ko2-13_3h_rep2/Aligned.out.bam > ${DIR}/ko2-13_3h_rep2/Aligned.SortedByCoordinate.bam"
bsub -q Z-ZQF -n 20 -e error "samtools sort --threads 20 -m 100G ${DIR}/ko2-13_3h_rep3/Aligned.out.bam > ${DIR}/ko2-13_3h_rep3/Aligned.SortedByCoordinate.bam"

########################
### DE analysis
########################

GTF=/150T/zhangqf/GenomeAnnotation/Gencode/mm10.gtf
ROOT1=/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-09-27-ko/Processing/3.mapping
ROOT2=/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-10-21-ko/Processing/3.mapping
ROOT3=/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/3.mapping
OUT=/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/4.cuffdiff

WT_0h_rep1=${ROOT1}/wt-0h-1/Aligned.SortedByCoordinate.bam
WT_0h_rep2=${ROOT1}/wt-0h-2/Aligned.SortedByCoordinate.bam
WT_0h_rep3=${ROOT2}/wt-0h-3/Aligned.SortedByCoordinate.bam
KO_0h_rep1=${ROOT3}/ko2-13_0h_rep1/Aligned.SortedByCoordinate.bam
KO_0h_rep2=${ROOT3}/ko2-13_0h_rep2/Aligned.SortedByCoordinate.bam
KO_0h_rep3=${ROOT3}/ko2-13_0h_rep3/Aligned.SortedByCoordinate.bam

WT_3h_rep1=${ROOT1}/wt-3h-1/Aligned.SortedByCoordinate.bam
WT_3h_rep2=${ROOT2}/wt-3h-2/Aligned.SortedByCoordinate.bam
WT_3h_rep3=${ROOT2}/wt-3h-3/Aligned.SortedByCoordinate.bam
KO_3h_rep1=${ROOT3}/ko2-13_3h_rep1/Aligned.SortedByCoordinate.bam
KO_3h_rep2=${ROOT3}/ko2-13_3h_rep2/Aligned.SortedByCoordinate.bam
KO_3h_rep3=${ROOT3}/ko2-13_3h_rep3/Aligned.SortedByCoordinate.bam

bsub -q Z-ZQF -e error -n 20 \
    "cuffdiff -o ${OUT}/0h -L wt,ko -p 20 \
        -min-alignment-count 10 \
        -FDR 0.05 \
        -library-type fr-unstranded \
        -library-norm-method classic-fpkm \
        -dispersion-method pooled \
        --no-update-check \
        ${GTF} \
        ${WT_0h_rep1},${WT_0h_rep2},${WT_0h_rep3} \
        ${KO_0h_rep1},${KO_0h_rep2},${KO_0h_rep3}"

bsub -q Z-ZQF -e error -n 20 \
    "cuffdiff -o ${OUT}/3h -L wt,ko -p 20 \
        -min-alignment-count 10 \
        -FDR 0.05 \
        -library-type fr-unstranded \
        -library-norm-method classic-fpkm \
        -dispersion-method pooled \
        --no-update-check \
        ${GTF} \
        ${WT_3h_rep1},${WT_3h_rep2},${WT_3h_rep3} \
        ${KO_3h_rep1},${KO_3h_rep2},${KO_3h_rep3}"

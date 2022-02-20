
function count_fq()
{
    infq=$1

    if [[ ! -f "$infq" ]]; then
        echo -e "${infq} does not exists"
        return -1
    fi

    declare cmd=""
    if [[ $(endswith "${infq}" ".gz") == "true" ]]; then
        cmd=${cmd}"gzip -dc "${infq}
    else
        cmd=${cmd}"cat "${infq}
    fi

    cmd=${cmd}" | wc -l"

    #echo $cmd
    count=$(eval ${cmd})
    echo $count
    return 0
}


function reads_statistics()
{
    # root-file-name
    OUT_DIR=$1

    # input_file
    FQ_FILE=$2

    # output prefix
    OUT_PREFIX=$3

    # ratio
    NUMBER=$4

    # adaptor
    ADAPTOR=$5

    # Check true parameters
    if [[ $# < 5 ]]; then
        echo "Usage: $0 out_dir fq_file out_prefix number adaptor [QUEUE] [THREADS]";
        return;
    fi

    if [[ "$(endswith "${FQ_FILE}" ".gz")" == "true" ]]; then
        echo -e "input file endswith .gz";
        total=$(gzip -dc ${FQ_FILE} | wc -l)
    else
        total=$(cat ${FQ_FILE} | wc -l)
    fi
    total=$(awk -v total="$total" 'BEGIN{print int(total/4);}')
    sample_rate=$(awk -v NUMBER="$NUMBER" -v total="$total" 'BEGIN{print NUMBER/total;}')

    if ((`bc <<< "${sample_rate}>1"`)); then
        sample_rate=1.0
    fi

    echo -e "Total ${total} reads. Want to sample ${NUMBER} reads. Sample rate=${sample_rate}"

    mkdir -p $OUT_DIR/0_sampleFQ
    mkdir -p $OUT_DIR/1_cutadapt
    mkdir -p $OUT_DIR/2_trimmomatic
    mkdir -p $OUT_DIR/3_readCollapse
    mkdir -p $OUT_DIR/4_cutadapt
    mkdir -p $OUT_DIR/5_rRNA
    mkdir -p $OUT_DIR/6_mapGenome

    SAMPLE_FQ=$OUT_DIR/0_sampleFQ/${OUT_PREFIX}.fastq.gz
    CUTADAP_1=$OUT_DIR/1_cutadapt/${OUT_PREFIX}.fastq.gz
    TRIMMO_2=$OUT_DIR/2_trimmomatic/${OUT_PREFIX}.fastq
    COLLAPSE_3=$OUT_DIR/3_readCollapse/${OUT_PREFIX}.fastq
    CUTADAP_4=$OUT_DIR/4_cutadapt/${OUT_PREFIX}.fastq.gz
    RRNA_SAM=$OUT_DIR/5_rRNA/${OUT_PREFIX}.sam
    RRNA_UNMAP=$OUT_DIR/5_rRNA/${OUT_PREFIX}.unmap.fastq

    # Cluster parameters

    THREADS=1
    QUEUE=Z-BNODE
    if (($# >= 6)); then QUEUE=$6; fi
    if (($# >= 7)); then THREADS=$7; fi

    PAR_THREADS=16

    job_name=$OUT_PREFIX
    Error=${OUT_PREFIX}_error
    Log=${OUT_PREFIX}_log
    rRNA_INDEX=/150T/zhangqf/lipan/SMART_SHAPE/ref_sequence/human_rRNA/human_rRNA
    genome_INDEX=/150T/zhangqf/GenomeAnnotation/INDEX/STAR/hg38_Gencode

    # Group id
    GROUP_ID="/reads_statistics/"$RANDOM
    echo -e "Current Group ID: "${GROUP_ID}

    ## Run

    BIN=/Share/home/zhangqf7/jinsong_zhang/icSHAPE/icSHAPE-master/scripts

    BSUB_0=$(bsub -n $THREADS -J ${job_name}_sample_fq -q $QUEUE -e $OUT_DIR/$Error -o $OUT_DIR/$Log -g $GROUP_ID \
        "seqtk seq -f $sample_rate $FQ_FILE | gzip - > $SAMPLE_FQ")
        #"sample_fq.py -i $FQ_FILE -o $SAMPLE_FQ -s $RATIO ")
    #BSUB_0=$(bsub -n $THREADS -J ${job_name}_cp_gzip -q $QUEUE -e $OUT_DIR/$Error -o $OUT_DIR/$Log -g $GROUP_ID \
    #    "cat $RATIO $FQ_FILE | gzip - > $SAMPLE_FQ")


    BSUB_1=$(bsub -n $THREADS -J ${job_name}_cutadapt_1 -w "done($(job_id $BSUB_0))" -q $QUEUE -e $OUT_DIR/$Error -o $OUT_DIR/$Log -g $GROUP_ID \
        "cutadapt -a \"${ADAPTOR}\" -m 25 $SAMPLE_FQ -o $CUTADAP_1")

    BSUB_2=$(bsub -n $THREADS -J ${job_name}_trimmomatic -w "done($(job_id $BSUB_1))" -q $QUEUE -e $OUT_DIR/$Error -o $OUT_DIR/$Log \
        "java -Xmx30g -jar /Share/home/zhangqf/shaodi/app/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads $PAR_THREADS -phred33 $CUTADAP_1 $TRIMMO_2 TRAILING:20 MINLEN:25")

    BSUB_3=$(bsub -n $THREADS -J ${job_name}_collapse -w "done($(job_id $BSUB_2))" -q $QUEUE -e $OUT_DIR/$Error -o $OUT_DIR/$Log \
        "perl $BIN/readCollapse.pl -U $TRIMMO_2 -o $COLLAPSE_3 -f /dev/null ")

    BSUB_4=$(bsub -n $THREADS -J ${job_name}_cutadapt_4 -w "done($(job_id $BSUB_3))" -q $QUEUE -e $OUT_DIR/$Error -o $OUT_DIR/$Log \
        "cutadapt -u 10 -m 25 $COLLAPSE_3 -o $CUTADAP_4")

    BSUB_5=$(bsub -n $PAR_THREADS -J ${job_name}_map_rRNA -w "done($(job_id $BSUB_4))" -q $QUEUE -e $OUT_DIR/$Error -o $OUT_DIR/$Log \
        "bowtie2 -U $CUTADAP_4 -x $rRNA_INDEX --sensitive -p $PAR_THREADS | \
            awk '{
                    if(substr(\$0,1,1)==\"@\"){ print \$0 > \"$RRNA_SAM\"; }
                    else if(\$2==0||\$2==16){ print \$0 > \"$RRNA_SAM\"; }
                    else if(\$2==4){ printf \"@%s\\n%s\\n+\\n%s\\n\",\$1,\$10,\$11 > \"$RRNA_UNMAP\"; }
                    else { print \"Impossible\!\!\"; }
                }'")

    BSUB_6=$(bsub -n $PAR_THREADS -J ${job_name}_map_genome -w "done($(job_id $BSUB_5))" -q $QUEUE -e $OUT_DIR/$Error -o $OUT_DIR/$Log \
        "STAR --outFileNamePrefix $OUT_DIR/6_mapGenome/${OUT_PREFIX}.  \
                --readFilesIn $RRNA_UNMAP \
                --genomeDir $genome_INDEX \
                --runThreadN $PAR_THREADS \
                --genomeLoad LoadAndKeep \
                --runMode alignReads \
                --outReadsUnmapped Fastx \
                --outSAMtype SAM \
                --outBAMsortingThreadN 6 \
                --outSAMmultNmax 1 \
                --outFilterMultimapNmax 10 \
                --outFilterMismatchNmax 2 \
                --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
                --outSAMstrandField intronMotif \
                --outSJfilterOverhangMin 30 12 12 12 \
                --alignEndsType EndToEnd \
                --outSAMattributes All \
                --alignIntronMin 20 \
                --alignIntronMax 1000000 \
                --alignMatesGapMax 1000000 \
                --alignSJDBoverhangMin 1")
}

function sliding_icSHAPE()
{
    # root-file-name
    OUT_DIR=$1

    # sam-file-name
    SAM_N1=$2
    SAM_N2=$3

    JUNCTION=$4
    SIZE=$5

    # /150T/zhangqf/GenomeAnnotation/Gencode/hg38.genomeCoor.bed
    ANNOTATION=$6

    OUT_PREFIX=$7

    # Check true parameters
    if [[ $# < 7 ]]; then
        echo "Usage: sliding_icSHAPE out_dir NAI_1 NAI_2 JUNCTION SIZE GenomeCoor_Annotation OUT_PREFIX [QUEUE] [THREADS]";
        return;
    fi
    
    # remove all files in $OUT_DIR
    mkdir -p $OUT_DIR
    #rm -rfi $OUT_DIR/* || true

    # Cluster parameters

    #JUNCFILE=/150T/zhangqf/lipan/SMART_SHAPE/ref_sequence/STAR_hg38_genome/sjdbList.fromGTF.out.tab
    #SIZEFILE=/150T/zhangqf/lipan/SMART_SHAPE/ref_sequence/STAR_hg38_genome/chrNameLength.txt

    QUEUE=Z-BNODE
    THREADS=1
    if (($# >= 8)); then QUEUE=$8; fi
    if (($# >= 9)); then THREADS=$9; fi

    job_name=$OUT_DIR/$OUT_PREFIX
    error_file=$OUT_DIR/${OUT_PREFIX}_error
    log_file=$OUT_DIR/${OUT_PREFIX}_log

    ## 1. sam2tab

    BSUB_1=$(bsub -n $THREADS -J ${job_name}_sam2tab_N1 -q $QUEUE -e $error_file -o $log_file \
    "sam2tab -in $SAM_N1 -out $OUT_DIR/${OUT_PREFIX}_NAI_1.tab")
    
    BSUB_2=$(bsub -n $THREADS -J ${job_name}_sam2tab_N2 -q $QUEUE -e $error_file -o $log_file \
    "sam2tab -in $SAM_N2 -out $OUT_DIR/${OUT_PREFIX}_NAI_2.tab")

    ## 2. calc_sliding_shape smartSHAPE

    BSUB_3=$(bsub -n $THREADS -J ${job_name}_sliding_smart -w "done($(job_id $BSUB_1))&&done($(job_id $BSUB_2))" -q $QUEUE -e $error_file -o $log_file \
    "calc_sliding_shape smartSHAPE -N $OUT_DIR/${OUT_PREFIX}_NAI_1.tab,$OUT_DIR/${OUT_PREFIX}_NAI_2.tab -ijf $JUNCTION -size $SIZE -out $OUT_DIR/${OUT_PREFIX}_smartSHAPE.tab")

    ## 3. CountRT
    BSUB_4=$(bsub -n $THREADS -J ${job_name}_countRT -w "done($(job_id $BSUB_1))&&done($(job_id $BSUB_2))" -q $QUEUE -e $error_file -o $log_file \
    "countRT -in $OUT_DIR/${OUT_PREFIX}_NAI_1.tab,$OUT_DIR/${OUT_PREFIX}_NAI_2.tab -size $SIZE -ijf $JUNCTION -out $OUT_DIR/${OUT_PREFIX}_RTBD.tab")

    ## 4. Genome Coordination to Transcript Coordination smartSHAPE
    BSUB_5=$(bsub -n $THREADS -J ${job_name}_genSHAPEToTransSHAPE -w "done($(job_id $BSUB_3))" -q $QUEUE -e $error_file -o $log_file \
    "genSHAPEToTransSHAPE.py -i $OUT_DIR/${OUT_PREFIX}_smartSHAPE.tab -t smartSHAPE -g $ANNOTATION -o $OUT_DIR/${OUT_PREFIX}_smartSHAPE.out -c 100")

    ## 5. Genome RTBD to Transcriptome RTBD
    BSUB_6=$(bsub -n $THREADS -J ${job_name}_genRTBDToTransRTBD -w "done($(job_id $BSUB_4))" -q $QUEUE -e $error_file -o $log_file \
    "genRTBDToTransRTBD.py -i $OUT_DIR/${OUT_PREFIX}_RTBD.tab -g $ANNOTATION -o $OUT_DIR/${OUT_PREFIX}_RTBD.out --col 4,6")
}

#############################
#### Step 1: Sample and Map
#############################

# Ratio for per 50M, we'll sample 50M, 100M, 150M, 200M, 250M
ratio_N_1=0.07539179157
#ratio_N_5=0.128083291
#ratio_N_25=0.165329018
#ratio_N_125=0.179975117



#N_1_r1=/150T/zhangqf/lipan/SMART_SHAPE/raw_data/combined_NAIN3_1/NAIN3-1-r1.fastq
#N_1_r2=/150T/zhangqf/lipan/SMART_SHAPE/raw_data/combined_NAIN3_1/NAIN3-1-r2.fastq
#N_5_r1=/150T/zhangqf7/meiling_sequencing_rawdata/low_input_icshape/smartshape-180522/NAIN3-5-r1.fastq
#N_5_r2=/150T/zhangqf7/meiling_sequencing_rawdata/low_input_icshape/smartshape-180522/NAIN3-5-r2.fastq

#N_25_r1=/150T/zhangqf7/meiling_sequencing_rawdata/low_input_icshape/smartshape-180517/NAIN3-25-r1.fastq.gz
#N_25_r2=/150T/zhangqf7/meiling_sequencing_rawdata/low_input_icshape/smartshape-180517/NAIN3-25-r2.fastq.gz
#N_125_r1=/150T/zhangqf7/meiling_sequencing_rawdata/low_input_icshape/smartshape-180517/NAIN3-125-r1.fastq.gz
#N_125_r2=/150T/zhangqf7/meiling_sequencing_rawdata/low_input_icshape/smartshape-180517/NAIN3-125-r2.fastq.gz


OUTROOT=/150T/zhangqf/lipan/SMART_SHAPE/raw_reads_all

adaptor=GGAAGAGCACACGTCTG # for smartSHAPE
#reads_statistics $OUTROOT $N_1_r1 N_1_r1 $(echo $ratio_N_1*1|bc)
#reads_statistics $OUTROOT $N_1_r2 N_1_r2 $(echo $ratio_N_1*1|bc)
#reads_statistics $OUTROOT $N_5_r1 N_5_r1 $(echo $ratio_N_5*5|bc)
#reads_statistics $OUTROOT $N_5_r2 N_5_r2 $(echo $ratio_N_5*5|bc)
#reads_statistics $OUTROOT $N_25_r1 N_25_r1 $(echo $ratio_N_25*5|bc)
#reads_statistics $OUTROOT $N_25_r2 N_25_r2 $(echo $ratio_N_25*5|bc)
#reads_statistics $OUTROOT $N_125_r1 N_125_r1 $(echo $ratio_N_125*5|bc)
#reads_statistics $OUTROOT $N_125_r2 N_125_r2 $(echo $ratio_N_125*5|bc)


icSHAPE_D_r1=/150T/zhangqf7/zxlin/RNAProbDB/raw_data/human/GSE74353/SRR3194440.fastq.gz
icSHAPE_D_r2=/150T/zhangqf7/zxlin/RNAProbDB/raw_data/human/GSE74353/SRR3194441.fastq.gz
icSHAPE_N_r1=/150T/zhangqf7/zxlin/RNAProbDB/raw_data/human/GSE74353/SRR3194444.fastq.gz
icSHAPE_N_r2=/150T/zhangqf7/zxlin/RNAProbDB/raw_data/human/GSE74353/SRR3194445.fastq.gz
OUTROOT=/150T/zhangqf/lipan/SMART_SHAPE/raw_reads_all

adaptor=AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG # for icSHAPE

# Ratio for per 50M, we'll sample 50M, 100M, 150M, 200M, 250M
number_list=(50000000 100000000 150000000 200000000 250000000)
for number in ${number_list[@]};
do
    half_num=$(awk -v number="${number}" 'BEGIN{print int(number/2)}')
    echo "half_num:"${half_num}
    reads_statistics $OUTROOT $icSHAPE_D_r1 icSHAPE_D_r1_$number $half_num ${adaptor} Z-ZQF 20
    reads_statistics $OUTROOT $icSHAPE_D_r2 icSHAPE_D_r2_$number $half_num ${adaptor} Z-ZQF 20
    reads_statistics $OUTROOT $icSHAPE_N_r1 icSHAPE_N_r1_$number $half_num ${adaptor} Z-ZQF 20
    reads_statistics $OUTROOT $icSHAPE_N_r2 icSHAPE_N_r2_$number $half_num ${adaptor} Z-ZQF 20
done



#############################
#### Step 2: Calulate SHAPE score
#############################

ROOT_IN=6_mapGenome
ANNO=/150T/zhangqf/GenomeAnnotation/Gencode/hg38.genomeCoor.bed
JUNC=/150T/zhangqf/lipan/SMART_SHAPE/ref_sequence/STAR_hg38_genome/sjdbList.fromGTF.out.tab
SIZE=/150T/zhangqf/lipan/SMART_SHAPE/ref_sequence/STAR_hg38_genome/chrNameLength.txt

sliding_icSHAPE 7_calcScore $ROOT_IN/N_1_r1.Aligned.out.sam $ROOT_IN/N_1_r2.Aligned.out.sam $JUNC $SIZE $ANNO N_1
sliding_icSHAPE 7_calcScore $ROOT_IN/N_5_r1.Aligned.out.sam $ROOT_IN/N_5_r2.Aligned.out.sam $JUNC $SIZE $ANNO N_5
sliding_icSHAPE 7_calcScore $ROOT_IN/N_25_r1.Aligned.out.sam $ROOT_IN/N_25_r2.Aligned.out.sam $JUNC $SIZE $ANNO N_25
sliding_icSHAPE 7_calcScore $ROOT_IN/N_125_r1.Aligned.out.sam $ROOT_IN/N_125_r2.Aligned.out.sam $JUNC $SIZE $ANNO N_125




sliding_icSHAPE 7_calcScore $ROOT_IN/N_1_r1.Aligned.out.sam $ROOT_IN/N_1_r2.Aligned.out.sam $JUNC $SIZE $ANNO N_1







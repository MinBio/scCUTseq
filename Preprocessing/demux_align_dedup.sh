#!/usr/bin/env bash

# AUTHOR: Luuk Harbers

export LC_ALL=C


# Parse options

# Help string
helps="
 usage: ./demux_align_dedup.sh [-h] -i inFastq -o outDir -s samplelist [-c cells] [-t threads] [-b barcodepattern]
 	[-r ref] [-m mismatches] [-k known] [-d tmpdir]
 Description:
  Demultiplex, run SamToFastq, Align reads using bwa mem, MergeBamAlignment recalibrate and deduplicate.
 Mandatory arguments:
  -i  inFastq multiplexed fastq
  -o  outDir	Output directory. Created if not found.
  -c  cells Number of cells in experiment
  -s  samplelist list of samples with barcode. 1 barcode per row.
 Optional arguments:
  -h  Show this help page.
  -t  threads	Number of threads for parallelization.
  -r  index	Path to BWA index.
  		Default: '/mnt/AchTeraD/Documents/references/hg19/hg19.fa'
  -b  barcode-pattern  Pattern of barcode (cutsite included in bc) and umi in sample. C=barcode, N=UMI
      Default: NNNNNNNNCCCCCCCCCCCC
  -m  mismatches Number of missmatches allowed in barcode and cutsite combined
      Default: 1
  -k  known  known-sites for base score recalibration
      Default: '/mnt/AchTeraD/Documents/references/vcf-files/dbsnp-b151_GRCh37p13-all.vcf.gz'
  -d  tmpdir
      Default: '/mnt/AchTeraD/tmp/'
"

# Default values
threads=1
ref="/mnt/AchTeraD/Documents/references/hg19/hg19.fa"
barcodepattern="NNNNNNNNCCCCCCCCCCCC"
mismatches=1
known='/mnt/AchTeraD/Documents/references/vcf-files/dbsnp-b151_GRCh37p13-all.vcf.gz'
tmpdir="/mnt/AchTeraD/tmp/"

while getopts ht:i:o:s:r:c:b:m:k:d: opt; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 1
		;;
		t)
			if [ 0 -ge "$OPTARG" ]; then
				echo -e "Enforcing a minimum of 1 thread.\n"
			else
				threads=$OPTARG
			fi
		;;
		i)
			in_fastq=$OPTARG
			if [ ! -e "$OPTARG" ]; then
			  echo "input fastq not found"
			  exit 1
			fi
		;;
		o)
			out_dir=$OPTARG
			if [ ! -d "$OPTARG" ]; then
				msg="Output folder not found, creating it."
        mkdir -p "${out_dir}/extracting/"
        mkdir -p "${out_dir}/demultiplexed/"
        mkdir -p "${out_dir}/bamfiles/"
			fi
		;;
    s)
      if [ -e "$OPTARG" ]; then
        samplelist=$OPTARG
			else
				msg="Invalid -s option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
    fi
    ;;
		r)
			if [ -e "$OPTARG" ]; then
				ref=$OPTARG
			else
				msg="Invalid -r option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
			fi
		;;
    c)
			if [ -e "$OPTARG" ]; then
				cells=$OPTARG
			fi
		;;
    b)
			if [ ! -e "$OPTARG" ]; then
				barcodepattern=$OPTARG
			fi
		;;
    m)
			if [ -e "$OPTARG" ]; then
				mismatches=$OPTARG
			fi
		;;
    k)
      if [ -e "$OPTARG" ]; then
        known=$OPTARG
			else
				msg="Invalid -k option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
    fi
    ;;
		d)
			if [ ! -d "$OPTARG" ]; then
				tmpdir=$OPTARG
			fi
		;;
    *) echo "usage: $0 [-v] [-r]" >&2
       exit 1
    ;;
	esac
done

# Check mandatory options
if [ -z "$out_dir" ]; then
	echo -e "$helps\n!!! Missing mandatory -o option.\n"
	exit 1
fi
if [ -z "$in_fastq" ]; then
  echo -e "$helps\n!!! Missing Mandatory -i option.\n"
  exit 1
fi

# Additional checks
if [ -z "$ref" ]; then
	ref="/mnt/AchTeraD/Documents/references/hg19/hg19.fa"
fi


# RUN ==========================================================================

set -o pipefail

#extra check on output subdirectories
if [ ! -d "${out_dir}/extracting/" ]; then
  mkdir -p "${out_dir}/extracting/"
fi
if [ ! -d "${out_dir}/demultiplexed/" ]; then
  mkdir -p "${out_dir}/demultiplexed/"
fi
if [ ! -d "${out_dir}/bamfiles/" ]; then
  mkdir -p "${out_dir}/bamfiles/"
fi

name=$(echo "${in_fastq}" | sed 's/.fastq.gz//g' | sed 's@.*/@@')

#run whitelist of barcode file
#echo "Getting whitelist of all barcodes"
#umi_tools whitelist -I "${in_fastq}" \
#                    -S "${out_dir}/extracting/${name}_whitelist.txt" \
#                    --bc-pattern "${barcodepattern}" \
#                    --set-cell-number "${cells}" \
#                    --error-correct-threshold "${mismatches}" &&

#extracting barcodes/umi's from read to header
echo "extracting barcodes and umis from read and appending to header"
umi_tools extract -I "${in_fastq}" \
                  -S "${out_dir}/extracting/${name}_extracted.fastq.gz" \
                  -L "${out_dir}/extracting/${name}_extract.log" \
                  --extract-method=regex \
                  --bc-pattern='(?P<umi_1>.{8})(?P<cell_1>.{8})CATG{s<=1}' \
                  --filter-cell-barcode \
                  --whitelist "${samplelist}" &&


#demultiplex fastq
echo "demultiplexing fastq"

/home/luukharbers/BBMap_38.76/bbmap/demuxbyname.sh  \
                  in="${out_dir}/extracting/${name}_extracted.fastq.gz" \
                  out="${out_dir}/demultiplexed/%.fastq.gz" \
                  outu="${out_dir}/demultiplexed/undetermined.fast.gz" \
                  delimiter=_ \
                  column=2 &&

#alignment

cd "${out_dir}" || exit 1
samples=$(find "demultiplexed/" -type f -name "*.fastq.gz" | sed 's/.fastq.gz//g' | sed 's/demultiplexed\///g')

#loop through files and go through alignment and deduplication
for sample in ${samples}; do

  java -Xmx40G -jar /home/luukharbers/picard.jar FastqToSam \
    FASTQ="demultiplexed/${sample}.fastq.gz" \
    OUTPUT="bamfiles/${sample}.sam" \
    READ_GROUP_NAME="$sample" \
    SAMPLE_NAME="${sample}" \
    LIBRARY_NAME="${sample}" \
    PLATFORM="Illumina" \
    TMP_DIR="$tmpdir" &&
  java -Xmx24G -jar /home/luukharbers/picard.jar MarkIlluminaAdapters \
    I="bamfiles/${sample}.sam" \
    O="bamfiles/${sample}.markedAdapters.sam" \
    M="bamfiles/${sample}.markedAdapters.metrics" \
    TMP_DIR="$tmpdir" &&

  #run samtofastq+bwa+mergebamalignment pipe
  java -Xmx40G -jar /home/luukharbers/picard.jar SamToFastq \
    I="bamfiles/${sample}.markedAdapters.sam" \
    F=/dev/stdout \
    CLIPPING_ATTRIBUTE=XT \
    CLIPPING_ACTION=2 \
    TMP_DIR="${tmpdir}" | \
  #align with bwa mem
  bwa mem -M -t "$threads" -p "$ref" /dev/stdin | \
  #merge bam alignment with uBAM
  java -Xmx40G -jar /home/luukharbers/picard.jar MergeBamAlignment \
    ALIGNED_BAM=/dev/stdin \
    UNMAPPED_BAM="bamfiles/${sample}.sam" \
    OUTPUT="bamfiles/${sample}.piped.bam" \
    REFERENCE_SEQUENCE="$ref" \
    CREATE_INDEX=true \
    ADD_MATE_CIGAR=true \
    CLIP_ADAPTERS=false \
    CLIP_OVERLAPPING_READS=true \
    INCLUDE_SECONDARY_ALIGNMENTS=true \
    MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
    ATTRIBUTES_TO_RETAIN=XS \
    TMP_DIR="$tmpdir" &&

  #remove intermediate bams
  rm "bamfiles/${sample}.sam";
  rm "bamfiles/${sample}.markedAdapters.sam";
  rm "bamfiles/${sample}.markedAdapters.metrics";

  #run alfred?
  /home/luukharbers/alfred/bin/alfred qc \
  -r ${ref} \
  -o "bamfiles/${sample}.all.tsv.gz" \
  "bamfiles/${sample}.piped.bam" &&

  #deduplication using umi_tools dedup and retaining reads with >29 mapping quality
  echo "deduplication"
  umi_tools dedup -I "bamfiles/${sample}.piped.bam" \
                  -S "bamfiles/${sample}.dedup_q30.bam" \
                  -L "bamfiles/${sample}.dedup_q30.log" \
                  --mapping-quality 30 &&


  #run BQSR+recalibration pipeline
#  gatk BQSRPipelineSpark \
#        -R "$ref" \
#        -I "$(pwd)/bamfiles/${sample}.dedup_q30.bam" \
#        --known-sites "$known" \
#        -O "$(pwd)/bamfiles/${sample}.final_q30.bam" \
#        --spark-master local["$threads"] &&

  /home/luukharbers/alfred/bin/alfred qc \
  -r ${ref} \
  -o "bamfiles/${sample}.dedup.tsv.gz" \
  "bamfiles/${sample}.dedup_q30.bam" &&

  echo "moving to next sample"
done

echo "All done."



#!/usr/bin/env bash

# AUTHOR: Luuk Harbers

export LC_ALL=C


# Parse options

# Help string
helps="
 usage: ./demux_align_dedup.sh [-h] -i inFastq -o outDir -s samplelist [-t threads] [-b barcodepattern]
 	[-r ref]
 Description:
  Demultiplex, run SamToFastq, Align reads using bwa mem, MergeBamAlignment recalibrate and deduplicate.
 Mandatory arguments:
  -i  inFastq multiplexed fastq
  -o  outDir	Output directory. Created if not found.
  -s  samplelist list of samples with barcode. 1 barcode per row.
 Optional arguments:
  -h  Show this help page.
  -t  threads	Number of threads for parallelization.
  -r  index	Path to BWA index.
  		Default: '/mnt/AchTeraD/Documents/references/hg19/hg19.fa'
  -b  barcode-pattern  Pattern of barcode (cutsite included in bc) and umi in sample. Must be in regex string,
      check umi_tools manual for instructions
      Default: '(?P<umi_1>.{8})(?P<cell_1>.{8})CATG{s<=1}'

"

# Default values
threads=1
ref="/mnt/AchTeraD/Documents/references/hg19/hg19.fa"
barcodepattern='(?P<umi_1>.{8})(?P<cell_1>.{8})CATG{s<=1}'


while getopts ht:i:o:s:r:c:b: opt; do
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
    b)
			if [ ! -e "$OPTARG" ]; then
				barcodepattern=$OPTARG
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

#extracting barcodes/umi's from read to header
echo "extracting barcodes and umis from read and appending to header"
umi_tools extract -I "${in_fastq}" \
                  -S "${out_dir}/extracting/${name}_extracted.fastq.gz" \
                  -L "${out_dir}/extracting/${name}_extract.log" \
                  --extract-method=regex \
                  --bc-pattern="${barcodepattern}" \
                  --filter-cell-barcode \
                  --whitelist "${samplelist}" &&


#demultiplex fastq
echo "demultiplexing fastq"

/home/luukharbers/BBMap_38.76/bbmap/demuxbyname.sh  \
                  in="${out_dir}/extracting/${name}_extracted.fastq.gz" \
                  out="${out_dir}/demultiplexed/%.fastq.gz" \
                  delimiter=_ \
                  column=2 \
                  Xmx=20G &&

#alignment

cd "${out_dir}" || exit 1
samples=$(find "demultiplexed/" -type f -name "*.fastq.gz" | sed 's/.fastq.gz//g' | sed 's/demultiplexed\///g')

#loop through files and go through alignment and deduplication
for sample in ${samples}; do

  #align with bwa mem and sort output
  bwa mem -M -t "$threads" -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
  -p "$ref" "demultiplexed/${sample}.fastq.gz" | \
  samtools sort -o "bamfiles/${sample}.piped.sorted.bam" &&

  #index sorted bamfile
  samtools index "bamfiles/${sample}.piped.sorted.bam" &&

  #deduplication using umi_tools dedup and retaining reads with >29 mapping quality
  echo "deduplication"
  umi_tools dedup -I "bamfiles/${sample}.piped.sorted.bam" \
                  -S "bamfiles/${sample}.dedup_q30.bam" \
                  -L "bamfiles/${sample}.dedup_q30.log" \
                  --mapping-quality 30 &&

  samtools index "bamfiles/${sample}.dedup_q30.bam" &&

  echo "moving to next sample"
done

echo "Alignment and deduplication finished. Running quality control.."

#run alfred on pre-deduplicated bam files
find . -type f -name '*.piped.sorted.bam' | parallel -I% --max-args 1 \
"/home/luukharbers/alfred/bin/alfred qc -r ${ref} -o %.all.tsv.gz %" &&

#run alfred on deduplicated bam files
find . -type f -name '*.dedup_q30.bam' | parallel -I% --max-args 1 \
"/home/luukharbers/alfred/bin/alfred qc -r ${ref} -o %.dedup.tsv.gz %" &&

#formatting final qc files
zgrep ^ME bamfiles/*all.tsv.gz | cut -f 2- | sed -n '1p;0~2p' > all.tsv
zgrep ^ME bamfiles/*dedup.tsv.gz | cut -f 2- | sed -n '1p;0~2p' > dedup.tsv
echo "All done! :)"


mkdir ../downsampling
mkdir ../downsampling/run
mkdir ../downsampling/bamfiles
mkdir ../downsampling/coverage

for INBAM in *.bam
do
    INFILEBASE=`basename $INBAM`
    INSAMPLEBASE=${INFILEBASE%.bam}
    COVHISTFILE=${INSAMPLEBASE}.covhist.txt
    OUTBAM=$INSAMPLEBASE.downsampled.bam
	
    # Downsampling
    # Save headers, which will not be included when downsampling
    printf "%s\n" "samtools view ${INSAMPLEBASE}.bam -H > ${INSAMPLEBASE}.sam" >> ../downsampling/run/headers_commands.sh
    # Downsample by randomly selecting a subset of rows besides the header, and
    # appending them to the header
    # Use the file itself as the seed for the random number generator
    printf "%s\n" "samtools view $INSAMPLEBASE.bam | shuf -n 800000 --random-source=$INSAMPLEBASE.bam >> $INSAMPLEBASE.sam" >> ../downsampling/run/reads_commands.sh
    # Sort the downsampled SAM file
    printf "%s\n" "samtools sort -o ../downsampling/bamfiles/$OUTBAM $INSAMPLEBASE.sam" >> ../downsampling/run/sort_commands.sh
    # Delete the unsorted SAM files
    printf "%s\n" "rm $INSAMPLEBASE.sam" >> ../downsampling/run/delete_commands.sh
    # Calculate coverage maximum read size of 50 to all reads
    printf "genomeCoverageBed -ibam ../downsampling/bamfiles/$OUTBAM -fs 50 -max 1 > ../downsampling/coverage/$COVHISTFILE\n" >> ../downsampling/run/genomecoveragebed_commands.sh
done

# Run jobs
parallel --jobs 30 < ../downsampling/run/headers_commands.sh
parallel --jobs 30 < ../downsampling/run/reads_commands.sh
parallel --jobs 30 < ../downsampling/run/sort_commands.sh
parallel --jobs 30 < ../downsampling/run/delete_commands.sh
parallel --jobs 30 < ../downsampling/run/genomecoveragebed_commands.sh
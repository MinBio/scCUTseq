import scyn

# create SCYN object
scyn_operator = scyn.SCYN()

# call cnv
# bam_dir is the input bam directory and output_dir is the output directory
scyn_operator.call("/mnt/AchTeraD/data/BICRO231/NZ86_new/bamfiles/",
"/mnt/AchTeraD/Documents/Projects/scCUTseq/SCYN/BICRO231/NZ86/"
)

# store cnv matrix to a csv file
#scyn_operator.cnv.to_csv('your file name')
### AneuFinder analysis script ###

# set the path to bam-file directory
inputdir <- "/mnt/AchTeraD/runextra/TCCGAGAT.dedup.bam"

# set path to analysis output
outputdir <- "/mnt/AchTeraD/runextra/"
dir.create(outputdir, recursive = T)

# set path configuration file: NULL, mouse (GRCm38), or human (GRCh38))
config <- "/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/Aneufinder/Aneufinder-dnacopy.config"

### RUN FROM HERE ###

# analyse samples
Aneufinder(inputfolder = inputdir,
           outputfolder = outputdir, 
           configfile = config)

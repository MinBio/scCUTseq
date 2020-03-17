### AneuFinder analysis script ###

# load AneuFinder
require(AneuFinder)

# set the path to bam-file directory
inputdir <- "/mnt/AchTeraD/data/BICRO213/bamfiles_final/"

# set path to analysis output
outputdir <- "/mnt/AchTeraD/Documents/Projects/scCUTseq/data/BICRO213/Aneufinder/MS21"
dir.create(outputdir, recursive = T)
# set path configuration file: NULL, mouse (GRCm38), or human (GRCh38))
config <- "/mnt/AchTeraD/Documents/Aneufinder/config-files/AneuFinder-dnacopy.config"

### RUN FROM HERE ###

# analyse samples
Aneufinder(inputfolder = inputdir,
           outputfolder = outputdir, 
           configfile = config)
  
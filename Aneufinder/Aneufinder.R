### AneuFinder analysis script ###
library(AneuFinder)
library(GenomeInfoDb)

# set the path to bam-file directory
inputdir <- "/mnt/AchTeraD/data/BICRO243/MS67/bamfiles/"

# set path to analysis output
outputdir <- "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/Aneufinder/BICRO243/MS67/"
dir.create(outputdir, recursive = T)

# set path configuration file: NULL, mouse (GRCm38), or human (GRCh38))
config <- "/mnt/AchTeraD/data/Aneufinder_refs/AneuFinder.config"

### RUN FROM HERE ### 

# analyse samples
Aneufinder(inputfolder = inputdir,
           outputfolder = outputdir,
           configfile = config, 
           numCPU = 20,
           hotspot.pval = NULL)
            
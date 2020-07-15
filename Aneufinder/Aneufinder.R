### AneuFinder analysis script ###
require(AneuFinder)

# set the path to bam-file directory
inputdir <- "/mnt/AchTeraD/data/BICRO229/NEB/"

# set path to analysis output
outputdir <- "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/Aneufinder/BICRO229/NZ74_NEB/"
dir.create(outputdir, recursive = T)

# set path configuration file: NULL, mouse (GRCm38), or human (GRCh38))
config <- "/mnt/AchTeraD/data/Aneufinder_refs/AneuFinder.config"

### RUN FROM HERE ### 

# analyse samples
Aneufinder(inputfolder = inputdir,
           outputfolder = outputdir,
           configfile = config, 
           numCPU = 15,
           hotspot.pval = NULL)
            
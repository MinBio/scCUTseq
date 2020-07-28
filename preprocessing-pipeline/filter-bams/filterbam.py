#!/usr/bin/python3

# Author: Luuk Harbers
# Date: 2020-07-27
# Email: l.j.w.harbers@gmail.com


import filterbam
import pandas as pd

def main():
	cutsites = pd.read_csv(
		"/mnt/AchTeraD/Documents/Projects/scCUTseq/cutsite-distribution/hg19-cutsites_fixed.bed", sep = "\t",
		names = ['chrom', 'start', 'end'],
		dtype = {'chrom': str, 'start': int, 'end': int})

	# Get chromosomes
	chroms = []
	for i in range(1, 23):
		chroms.append(str(i))

	chroms.append("X")
	chroms.append("Y")


	# Loop over chromosomes of bamfile and cutsites bedfile
	filterbam.filterBam(
		inbam = "/mnt/AchTeraD/data/BICRO232/NZ84_small_newpipeline/bamfiles/TCACATCA.dedup_q30.bam",
		outbam = "/mnt/AchTeraD/data/BICRO232/NZ84_small_newpipeline/test.bam",
		chroms = chroms,
		max_dist = 20,
		cutsites = cutsites,
		update = 100)

if __name__ == '__main__':
	main()
#!/usr/bin/python3

# Author: Luuk Harbers
# Date: 2020-07-18
# Email: l.j.w.harbers@gmail.com

cimport cython
import gzip
from pysam import FastxFile

# Functions
# Hamming
def cHamming(str s0, str s1):
	cdef:
		int N = len(s0)
		int i, count = 0
	for i in range(N):
		count += (s0[i] != s1[i])
	return count

def hammingDistanceLoop(str barcode, barcodelist, int mismatches):
	for bc in barcodelist:
		if cHamming(bc, barcode) <= mismatches:
			return bc
	return "unassigned"

# Open file
def openfile(filename, mode = "r"):
	if filename.endswith(".gz"):
		return gzip.open(filename, mode + "t")
	else:
		return open(filename, mode)

# Do record processing
def process_fastq(entry, barcodes, mismatches):
	umi = entry.sequence[0:8]
	barcode = entry.sequence[8:16]
	match = hammingDistanceLoop(barcode, barcodes.iloc[:,0], mismatches)
	return (match, f"@{entry.name}_{umi}_{match} {entry.comment}\n{entry.sequence[20:]}\n+\n{entry.quality[20:]}\n")

def iterateFastq(str fastq, barcodes, int mismatches, int update, str outdir):
	# Initialize counts
	cdef:
		int match_count = 0
		int unassigned_count = 0

	with FastxFile(fastq) as input_handle:
		for entry in input_handle:
			result = process_fastq(entry, barcodes, mismatches)

			# Update counts
			if result[0] == "unassigned":
				unassigned_count += 1
			else:
				match_count += 1

			# Progress updates
			if (unassigned_count + match_count) % update == 0:
				print(f"Reads processed: {unassigned_count + match_count}")

			# Write output to file
			with open(outdir + result[0] + ".fastq", "a+") as outfile:
				outfile.write(result[1])
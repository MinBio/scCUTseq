#!/usr/bin/python3

# Author: Luuk Harbers
# Date: 2020-07-18
# Email: l.j.w.harbers@gmail.com

import pandas as pd
import gzip, argparse

# ARGPARSER time
parser = argparse.ArgumentParser(
	description="""Moves umi and barcode to header
	and demultiplexes fastq file into separate fastq files""")
parser.add_argument("-f", "--fastq", type=str, help="Path to fastq file.")
parser.add_argument(
	"-o", "--outdir", type=str,
	help="Path to output directory.")
parser.add_argument("-l", "--log", type=str, help="Path to logfile.")
parser.add_argument("-b", "--barcodes", type=str, help="Path to barcode file.")
parser.add_argument("-m", "--mismatches", type=int, default=1,
	help="Number of mismatches allowed for barcode matching, default = 1")
parser.add_argument(
	"-u", "--update", type=int, default=1e6,
	help="Update to stdout every n lines, default = 1.000.000")

args = parser.parse_args()

# Functions
# Get hamming distance
def hammingDistance(barcode, barcodelist, mismatches):
	for bc in barcodelist:
		if sum(i != j for i, j in zip(barcode, bc)) <= mismatches:
			return bc
	return "unassigned"

# Open file
def openfile(filename, mode = "r"):
	if filename.endswith(".gz"):
		return gzip.open(filename, mode + "t")
	else:
		return open(filename, mode)

# Make lines into fastq record
def makeRecord(lines):
	ks = ['name', 'sequence', 'optional', 'quality']
	return {k: v for k, v in zip(ks, lines)}

# Do record processing
def process_fastq(record, barcodes, mismatches):
	umi = record['sequence'][0:8]
	barcode = record['sequence'][8:16]

	match = hammingDistance(barcode, barcodes.iloc[:,0], mismatches)

	# Get new read title
	space_index = record['name'].index(" ")
	new_name = record['name'][:space_index] + "_" + umi + "_" + match + record['name'][space_index:]

	# trim reads
	new_sequence = record['sequence'][20:]
	new_quality = record['quality'][20:]

	return [match, f"@{new_name}\n{new_sequence}\n+\n{new_quality}\n"]

# Script

# Read in barcode file
barcodes = pd.read_csv(args.barcodes, delimiter = ",")

# Initialize counts
match_count = 0
unassigned_count = 0

n = 4
with openfile(args.fastq, "r") as input_handle:
	# Get generator
	lines = []
	for line in input_handle:
		lines.append(line.rstrip())
		if len(lines) == n:
			record = makeRecord(lines)

			result = process_fastq(record, barcodes, args.mismatches)

			# Update counts
			if result[0] == "unassigned":
				unassigned_count += 1
			else:
				match_count += 1

			# Progress updates
			if (unassigned_count + match_count) % args.update == 0:
				print(f"Reads processed: {unassigned_count + match_count}")

			# Write output to file
			with open(args.outdir + result[0] + ".fastq", "a+") as outfile:
				outfile.write(result[1])
			lines = []

# Write counts
with open(args.log, "w") as logfile:
	logfile.write(
		f"args used :{args}\n\n"
		f"total reads: {match_count + unassigned_count}\n"
		f"reads with barcode: {match_count}\n"
		f"read unassigned: {unassigned_count}")
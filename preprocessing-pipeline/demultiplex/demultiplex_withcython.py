#!/usr/bin/python3

# Author: Luuk Harbers
# Date: 2020-07-18
# Email: l.j.w.harbers@gmail.com

import pandas as pd
import gzip, argparse
from pysam import FastxFile
import demultiplex

# ARGPARSER
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

# Script

# Read in barcode file
barcodes = pd.read_csv(args.barcodes, delimiter = ",")

demultiplex.iterateFastq(args.fastq, barcodes, args.mismatches, args.update, args.outdir)


# Write counts
with open(args.log, "w") as logfile:
	logfile.write(
		f"args used :{args}\n\n"
		f"total reads: {match_count + unassigned_count}\n"
		f"reads with barcode: {match_count}\n"
		f"read unassigned: {unassigned_count}")
#!/usr/bin/python3

# Author: Luuk Harbers
# Date: 2020-07-27
# Email: l.j.w.harbers@gmail.com

import pysam
from bisect import bisect_left
import math
import numpy as np
cimport cython
cimport numpy as np
# Functions
'''
def take_closest(myList, int myNumber):
	"""
	Assumes myList is sorted. Returns closest value to myNumber.

	If two numbers are equally close, return the smallest number.
	"""
	cdef int pos = bisect_left(myList, myNumber)
	if pos == 0:
		return myList[0]
	if pos == len(myList):
		return myList[-1]
	cdef int before = myList[pos - 1]
	cdef int after = myList[pos]
	if after - myNumber < myNumber - before:
		return after
	else:
		return before
'''

def take_closest(np.ndarray myList, int myNumber):
	cdef int idx = np.searchsorted(myList, myNumber, side="left")
	if idx > 0 and (idx == len(myList) or math.fabs(myNumber - myList[idx-1]) < math.fabs(myNumber - myList[idx])):
		return myList[idx-1]
	else:
		return myList[idx]

def filterBam(
	str inbam, str outbam, chroms,
	int max_dist, cutsites, int update):
	cdef:
		int inrange, outrange, closest, distance = 0

	# Open input and output bamfile
	bam_in = pysam.AlignmentFile(inbam, "rb")
	bam_out = pysam.AlignmentFile(outbam, "wb", template = bam_in)

	# Loop through chromosomes and write reads that fall within
	# Max_dist of cutsite
	for chrom in chroms:
		for read in bam_in.fetch(chrom):
			closest = take_closest(cutsites[cutsites.chrom == chrom]["start"].values,
								   			read.pos)
			distance = abs(closest - read.pos)
			# Output reads within max_dist distance of cutsite
			if read.is_reverse:
				if (abs(distance - 57)) <= max_dist:
					bam_out.write(read)
					inrange += 1
				else:
					outrange += 1
			else:
				if (abs(distance - 3)) <= max_dist:
					bam_out.write(read)
					inrange += 1
				else:
					outrange += 1

			# Print progress
			if (inrange + outrange) % update == 0:
				print(f"processed {inrange} reads")
	bam_in.close()
	bam_out.close()

#!/usr/bin/env python
"""
========================================================
Extract info on reads that align to a given region
in draft genome assembly.
========================================================
"""
from __future__ import print_function

try:
	from Bio import SeqIO
	import pysam
	import argparse
	import subprocess
	import tarfile
	import gzip
	import sys,os
except ImportError:
	print('Missing package(s)')
	quit()

verbose = False
log = list()

def main():
	# --------------------------------------------------------
	# PART 0: 	Parse input
	# --------------------------------------------------------
	parser = argparse.ArgumentParser(description='Extract and package reads within region')
	parser.add_argument('-v', '--verbose', action="store_true", default=False, required=False, dest="verbose", help="Use for verbose output with info on progress.")
	parser.add_argument('-b', '--bam', action="store", required=True, dest="bam", help="Sorted bam file created by aligning reads to the draft genome (refer to reads.sorted.bam in Nanopolish README).")
	parser.add_argument('-r', '--reads', action="store", dest="fa_filename", help="Fasta, fastq, fasta.gz, or fastq.gz file (refer to reads.fa in Nanopolish README)")
	parser.add_argument('-g', '--genome',  action="store", required=True, dest="draft_ga", help="Draft genome assembly (refer to draft.fa in Nanopolish README).")
	parser.add_argument('-w', '--window', action="store", required=True, dest="draft_ga_coords", help="Draft genome assembly coordinates wrapped in quotes ex. \"tig000001:10000-20000\".")
	parser.add_argument('-o', '--output_prefix', action="store", required=False, default="reads_subset", dest="output_prefix", help="Output prefix for tar.gz file and log file.")
	args = parser.parse_args()
	
	# Check to see if user used verbose option
	global verbose
	if args.verbose:
		verbose = True

	# Infer readdb file from fasta/q file
	readdb = args.fa_filename + ".index.readdb"

	custom_print( "===================================================" )
	custom_print( "Extract reads that align to given region" )
	custom_print( "Package all necessary files to reproduce error" )
	custom_print( "===================================================" )

	# --------------------------------------------------------
	# PART 1: 	Validate input
	# --------------------------------------------------------
	custom_print( "[ Input ]" )
	custom_print( "[+] Extracting from draft genome assembly coords: " + args.draft_ga_coords )
	custom_print( "[+] BAM file (reads.fa aligned to draft.fa): " + args.bam )
	custom_print( "[+] Readdb file: " + readdb )
	custom_print( "[+] Draft genome assembly (draft.fa): " + args.draft_ga )
	custom_print( "[+] FASTA/Q file (reads.fa): " + args.fa_filename )
	custom_print( "[+] Output prefix: " + args.output_prefix ) 

	custom_print( "[ Input check ]" )
	files = list()
	files.append(args.bam)
	files.append(readdb)
	files.append(args.fa_filename)
	files.append(args.draft_ga)
	draft_ga_fai = args.draft_ga + ".fai"
	files.append(draft_ga_fai)

	for i in files:
		if not os.path.exists(i) or not os.path.getsize(i) > 0 or not os.access(i, os.R_OK):
			print( "Expecting " + i + ". But does not exist, is empty or is not readable." )
			sys.exit(1)

	custom_print( "[ Validated input ] All input files exist, are not-empty, and are readable." )

	# --------------------------------------------------------
	# PART 2: 	Reassign input argument values	
	# --------------------------------------------------------
	# o = old/original, ga = genome assembly, fa = fasta/q file
	# coords = coordinates, op = output
	o_bam = args.bam
	o_readdb = readdb
	o_fa = args.fa_filename
	op = args.output_prefix
	draft_ga_coords = args.draft_ga_coords

	# --------------------------------------------------------
	# PART 3: 	With user input ref coords, extract all 
	#		aligned reads within these coordinates, 
	#		store read_ids, and fast5 files.
	# --------------------------------------------------------
	custom_print( "[ Extracting info on reads aligned to region ] \t" + draft_ga_coords )
	samfile = pysam.AlignmentFile(o_bam, "rb")
	region_read_ids = list()
	region_num_reads = 0

	# get all read ids of reads that are aligned to region in draft assembly
	for read in samfile.fetch(region=draft_ga_coords):
		id = read.query_name
		# add to list if not already in list
		if not id in region_read_ids:
			# store read id in list
			region_read_ids.append(id)
			# count number of reads that were aligned to the given region
			region_num_reads+=1

	# --------------------------------------------------------
	# PART 4:   Parse readdb file and find path to fast5 files
	# 		associated with each read that aligned to region
	# --------------------------------------------------------
	# readdb file has 2 columns: one indicating read_id and another indicating the fast5 file the read came from
	# each row represents a read
	custom_print( "[ Reading readdb file ]" )
	region_fast5_files = dict()
	with open (o_readdb, "r") as file:
		for line in file:
			l = line.split("\t")
			read_id = l.pop(0)
			if read_id in region_read_ids:
				fast5_file = l.pop(0)
				region_fast5_files[str(read_id)] = fast5_file.rstrip()

	# --------------------------------------------------------
	# PART 5:   Make a region BAM and BAI file
	# --------------------------------------------------------
	new_bam = "reads.bam"
	custom_print( "[ Writing to a new BAM file ] \t" + new_bam )
	region_reads = pysam.view("-b", o_bam, draft_ga_coords, "-o", new_bam, catch_stdout=False)
	
	new_bam_index = new_bam + ".bai"
	custom_print( "[ Writing to a new BAI file ] \t" + new_bam_index )
	pysam.index(new_bam, new_bam_index)

	# --------------------------------------------------------
	# PART 6: 	With user input ref coords, extract all 
	#		aligned	reads within these coordinates 
	#		and make new FASTA file
	# --------------------------------------------------------
	# detect type of sequences file then handle accordingly
	file_type = detect_fa_filetype(o_fa)
	new_fa = "reads.fasta"
	custom_print( "[ Writing to a new fasta file ]\t" +  new_fa )
	with open (new_fa, "w") as fout:
		if ".gz" in file_type:
			with gzip.open(o_fa, "rt") as fin:
				if "fasta.gz" in file_type:
					for record in SeqIO.parse(fin, "fasta"):
						if record.id in region_read_ids:
							fout.write(">" + record.id + "\n")
							fout.write(str(record.seq) + "\n")
				elif "fastq.gz" in file_type:
					for record in SeqIO.parse(fin, "fastq"):
						if record.id in region_read_ids:
							fout.write(">" + record.id + "\n")
							fout.write(str(record.seq) + "\n")
		else:
			with open(o_fa, "rt") as fin:
				if "fasta" in file_type:
					for record in SeqIO.parse(fin, "fasta"):
						if record.id in region_read_ids:
							fout.write(">" + record.id + "\n")
							fout.write(str(record.seq) + "\n")
				elif "fastq" in file_type:
					for record in SeqIO.parse(fin, "fastq"):
						if record.id in region_read_ids:
							fout.write(">" + record.id + "\n")
							fout.write(str(record.seq) + "\n")

	# --------------------------------------------------------
	# PART 7: 	Let's get to tarring
	# --------------------------------------------------------
	# While tarring, we need to fix the directory structure
	# such that the original path to files are not saved.
	# For each fast5 file we need to extract the basename,
	# and save it in tar such that we save only the basename,
	# and not the whole path from the original source.
	tar_filename = op + ".tar.gz"
	archive = tarfile.open(tar_filename, "w:gz")
	custom_print( "[ Creating a tar.gz file ] \t" + tar_filename )
	custom_print( "[+] FAST5 files: " + op + "/fast5_files/<FAST5 file(s)>" )
	# track missing fast5 files
	bad_f5_found = False # true if missing fast5 file
	bad_read_id = ""
	bad_f5_path = ""
	num_bad_cases = 0
	for r in region_fast5_files.keys():
		read_id = r
		f5 = region_fast5_files[r]

		# get basename of fast5 file
		f5_basename = extract_basename(f5)
		an = op + "/fast5_files/" + f5_basename
		try:
			archive.add(f5, arcname=an)
		except:
			bad_f5_found = True
			bad_read_id = read_id
			bad_f5_path = f5
			num_bad_cases += 1
	
	# handle missing fast5 files
	if bad_f5_found:
		print("\nERROR: For read " + read_id + ", could not add " + str(f5) + ".")
		print("This path is inferred from the readdb file.")
		print("Please check that this is the correct path in readdb file for this read.")
		if num_bad_cases > 1:
			print("There are " + str(num_bad_cases) + " other reads with this problem (out of " + str(len(region_fast5_files)) + ").")
		print("\n")
		sys.exit(1)

	# --------------------------------------------------------
	# PART 8:	Add new files to tar
	# 			new fasta, new bam, and new bai with reads 
	#			in the region given only
	# --------------------------------------------------------
	an = op + "/" + new_fa
	archive.add(new_fa, arcname=an)
	custom_print( "[+] New FASTA: " + an )
	
	an_new_bam = op + "/" + new_bam
	archive.add(new_bam, arcname=an_new_bam)
	custom_print( "[+] New BAM: " + an_new_bam )

	an_new_bam_index = op + "/" + new_bam_index
	archive.add(new_bam_index, arcname=an_new_bam_index)
	custom_print( "[+] New BAI: " + an_new_bam_index )

	# --------------------------------------------------------
	# PART 9:	Add original draft genome assembly file
	#			and the index file
	# --------------------------------------------------------
	an_draft_ga = op + "/draft.fa"
	archive.add(args.draft_ga, arcname=an_draft_ga)
	custom_print( "[+] Original draft ga: " + an_draft_ga )

	an_draft_ga_fai = op + "/draft.fa.fai"
	archive.add(i, arcname=an_draft_ga_fai)
	custom_print( "[+] Original draft ga index: " + an_draft_ga_fai )

	# --------------------------------------------------------
	# PART 10: 	Check the number of reads in all new files
	# --------------------------------------------------------
	custom_print( "[ Output check ] " )
	# check the length of bam file
	num_reads_bam = region_num_reads
	num_reads_fasta = int(float(file_length(new_fa))/2.0)
	num_fast5_files = len(region_fast5_files)
	values = list()
	values.append(num_reads_bam)
	values.append(num_reads_fasta)
	custom_print( "[+] Num reads in new BAM: \t" + str(num_reads_bam) )
	custom_print( "[+] Num reads in new FASTA: \t" + str(num_reads_fasta) )
	custom_print( "[+] Num files in fast5_files/: \t" + str(num_fast5_files))
	if not all( v == num_fast5_files for v in values ):
		print( "[!] WARNING: The number of reads in the new bam, new fasta, and num of fast5 files tarred are not equal..." )
	else:
		custom_print( "[ Validated output ] Number of reads in the new bam, new fasta, and num of fast5 files tarred are equal!" )

	# --------------------------------------------------------
	# FINAL: 	Output log if verbose flag not used
	# --------------------------------------------------------
	global log
	logfile = op + ".log"
	with open (logfile, "w") as lfile:
		for s in log:
			lfile.write(s + "\n")
	an_logfile = op + "/" + logfile
	custom_print( "[ Log file ] " +  an_logfile )
	custom_print( "[ Tar file ] " + str(tar_filename) )
	custom_print( "[ Finished ] " )
	archive.add(logfile, arcname=an_logfile)
	archive.close()

def file_length(filename):
	# ========================================================
	# Returns number of lines in a file
	# --------------------------------------------------------
	# Input: 	Filename
	# Output: 	Number of lines in the file ...
	# ========================================================
	with open(filename) as f:
		for i, l in enumerate(f):
			pass
	return int(i) + 1

def extract_basename(filename):
	# ========================================================
	# Returns base filename
	# --------------------------------------------------------
	# Input: 	Filenames with paths
	# Output: 	Base filename
	# ========================================================
	# remove backslashes at the end of the file names that could return empty basenames..
	a = filename.rstrip("\\")
	a = a.rstrip("//")
	b = os.path.basename(a)
	return str(b)

def detect_fa_filetype(fa_filename):
	# ========================================================
	# Detects filetype of sequences input
	# --------------------------------------------------------
	# Input: 	FASTA/Q filename
	# Output: 	Either ['fa.gz', 'fastq.gz', 'fasta.gz', 
	# 			'fastq', 'fasta']
	# ========================================================
	path = fa_filename
	if path.endswith('fa.gz'):
		print("Possibly using the reads file generated by nanopolish index? Use original reads file...")	
	for ext in ['fastq.gz', 'fasta.gz', 'fastq', 'fasta']:
		if path.endswith(ext):
			return ext
	print("Must be either fasta, fastq, fasta.gz, fastq.gz")
	sys.exit(1)

def custom_print(s):
	# ========================================================
	# Depending on verbose flag, will save all prints to 
	# log list, or will print to stdout
	# --------------------------------------------------------
	# Input: 	string to print
	# ========================================================
	global verbose
	global log
	if verbose:
		print(s)
	log.append(s)

if __name__ == "__main__":
	main()

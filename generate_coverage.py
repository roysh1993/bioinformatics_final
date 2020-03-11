import os
import sys
import random
import subprocess
import subprocess
import assembly_stats
OUTPUT_FRAG_1_FILE =  "data/output_frag_1.fastq"
OUTPUT_FRAG_2_FILE =  "data/output_frag_2.fastq"
SPADES_EXE_LOCATION = "/data/roy/bio_informatics/SPAdes-3.12.0-Linux/bin/spades.py"
OUTPUT_DIR = "./output_dir"



def count_fragments(fastq_files):
	frag_count = 0
	
	with open(fastq_files[0]) as f:
		i = 0
		for line in f:
			# [1] Start of fragment
			if i == 0 and line.rstrip().startswith('@'):
				i = 1
			else:
				# [2] sequence part of fragment
				if i == 1 and line.rstrip().isupper():
					i = 2
				else:
					# [3] + part of the fragment
					if i == 2 and line.rstrip() is '+':
						i = 3
					# [4] quality part of fragment
					else: 
						if i == 3:
							i = 0
							frag_count += 1
						else:
							i = 0

	return frag_count;



def calculate_new_fragment_count(total_frag_count,new_converage):
	return round(total_frag_count*new_converage)


def get_random_fragments(fastq_files,total_frag_count,new_frag_count):
	frag_random_indexes = sorted(random.sample(range(1, total_frag_count), new_frag_count))
	print ("random indices:", frag_random_indexes)
	frag_identifiers = get_fragments_by_index(fastq_files[0],frag_random_indexes)
	print("##########")
	print (frag_identifiers)
	print("##########")

	get_fragments_by_identifier(fastq_files[1],frag_identifiers)

	

def get_fragments_by_index(fastq_file,frag_indexes, output_filename= OUTPUT_FRAG_1_FILE):
	frag_count = 0
	frag_identifiers = []

	frag_identifier = ""
	frag_data = ""

	with open(fastq_file) as f, open(output_filename,"w") as out:
		i = 0
		for line in f:
			# [1] Start of fragment
			if i == 0 and line.rstrip().startswith('@'):
				i = 1
				frag_identifier = line.rstrip()
				frag_data = line.rstrip()
			else:
				# [2] sequence part of fragment
				if i == 1 and line.rstrip().isupper():
					i = 2
					frag_data += "\n" + line.rstrip()
				else:
					# [3] + part of the fragment
					if i == 2 and line.rstrip() is '+':
						i = 3
						frag_data += "\n" + line.rstrip()
					# [4] quality part of fragment
					else: 
						if i == 3:
							frag_data += "\n" + line.rstrip()
							i = 0
							frag_count += 1

							if frag_indexes:
								if frag_count == frag_indexes[0]:
									print(frag_indexes.pop(0))
									out.write(frag_data)
									out.write("\n")
									frag_identifiers.append(frag_identifier[:-1])

						else:
							i = 0
							frag_data = ""

	return frag_identifiers



def get_fragments_by_identifier(fastq_file,frag_identifiers,output_filename =OUTPUT_FRAG_2_FILE ):
	frag_count = 0

	frag_data = ""

	with open(fastq_file) as f,open(output_filename,"w") as out:
		print("HI")
		i = 0
		for line in f:
			# [1] Start of fragment
			if i == 0 and line.rstrip().startswith('@'):
				i = 1
				frag_identifier = line.rstrip()
				frag_data = line.rstrip()
			else:
				# [2] sequence part of fragment
				if i == 1 and line.rstrip().isupper():
					i = 2
					frag_data += "\n" + line.rstrip()
				else:
					# [3] + part of the fragment
					if i == 2 and line.rstrip() is '+':
						i = 3
						frag_data += "\n" + line.rstrip()
					# [4] quality part of fragment
					else: 
						if i == 3:
							frag_data += "\n" + line.rstrip()
							i = 0
							frag_count += 1

							if frag_identifiers:
								if frag_identifier == frag_identifiers[0]+"2":
									print(frag_identifiers.pop(0) + "2")
									out.write(frag_data)
									out.write("\n")

						else:
							i = 0
							frag_data = ""


def single_ineration_per_corr(fastq_files, coverage_ratio):
	total_reads_count = count_fragments(fastq_files)
	new_reads_count = calculate_new_fragment_count(total_reads_count,coverage_ratio)
	# get samples
	get_random_fragments(fastq_files,total_reads_count,new_reads_count)
	#run spades on sampled data
	subprocess.call(["python3",SPADES_EXE_LOCATION,"-1", fastq_files[0], -2, fastq_files[1],"-o" ,OUTPUT_DIR])
# 	get stats of per coverage
	stats = assembly_stats.calc_stats(os.path.join(OUTPUT_DIR, "scaffolds.fasta"))
# 	TODO delete later
	print(stats)
	return stats

#






if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description='Generate new converage dataset',
									 formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('--1', type=str, metavar='<fastq1>',
						help='Specify fastq files for recalculation, comma seperated.',
						required=True)
	parser.add_argument('--2', type=str, metavar='<fastq2>',
						help='Specify fastq files for recalculation, comma seperated.',
						required=True)

	parser.add_argument('--coverage', type=float, metavar='<coverage>',
						help='Specify new coverage.',
						required=True)

	args = parser.parse_args()
	fastq_files = [args.fastq1, args.fastq2]
	print(fastq_files)
	total_frag_count = count_fragments(fastq_files)
	print(total_frag_count)
	new_frag_count = calculate_new_fragment_count(total_frag_count,args.coverage)
	print(new_frag_count)
	get_random_fragments(fastq_files,total_frag_count,new_frag_count)

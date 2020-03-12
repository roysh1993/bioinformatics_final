import os
import sys
import random
import subprocess
import subprocess
import assembly_stats
import numpy as np
import matplotlib.pyplot as plt
import shutil
import json
OUTPUT_FRAG_1_FILE =  "data/output_frag_1.fastq"
OUTPUT_FRAG_2_FILE =  "data/output_frag_2.fastq"
SPADES_EXE_LOCATION = "/data/roy/bio_informatics/SPAdes-3.12.0-Linux/bin/spades.py"
OUTPUT_DIR = "output_dir"


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
	return int(total_frag_count*new_converage)


def get_random_fragments(fastq_files,total_frag_count,new_frag_count):
	frag_random_indexes = sorted(random.sample(range(1, total_frag_count + 1), new_frag_count))
	# print ("random indices:", frag_random_indexes)
	frag_identifiers = get_fragments_by_index(fastq_files[0],frag_random_indexes)
	# print("##########")
	# print (frag_identifiers)
	# print("##########")

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
									# print(frag_indexes.pop(0))
									out.write(frag_data)
									out.write("\n")
									frag_indexes.pop(0)
									frag_identifiers.append(frag_identifier[:-1])


						else:
							i = 0
							frag_data = ""

	return frag_identifiers



def get_fragments_by_identifier(fastq_file,frag_identifiers,output_filename =OUTPUT_FRAG_2_FILE ):
	frag_count = 0

	frag_data = ""

	with open(fastq_file) as f,open(output_filename,"w") as out:
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
									out.write(frag_data)
									out.write("\n")
									frag_identifiers.pop(0)


						else:
							i = 0
							frag_data = ""


def single_ineration_per_corr(fastq_files, coverage_ratio,output_dir=OUTPUT_DIR, DELETE_FILES= True):
	if coverage_ratio != 1.0:
		semgented_file1 = OUTPUT_FRAG_1_FILE
		semgented_file2 = OUTPUT_FRAG_2_FILE
		total_reads_count = count_fragments(fastq_files)
		new_reads_count = calculate_new_fragment_count(total_reads_count,coverage_ratio)
		# get samples
		get_random_fragments(fastq_files,total_reads_count,new_reads_count)
		#run spades on sampled data
	else:
		semgented_file1 = "/data/roy/bio_informatics/Staphylococcus_aureus/frag_1.fastq"
		semgented_file2 = "/data/roy/bio_informatics/Staphylococcus_aureus/frag_2.fastq"

	res = os.system("python3 " + SPADES_EXE_LOCATION +" -1 " + semgented_file1 + " -2 " + semgented_file2 + " -o " + output_dir)
# 	get stats of per coverage
	stats = assembly_stats.calc_stats(os.path.join(output_dir, "scaffolds.fasta"))
	if DELETE_FILES:
		shutil.rmtree(output_dir)
		os.mkdir(output_dir)

# 	TODO add maybe more fields to stats
	return stats

def generate_plots(stats_list,N50_list, directory = OUTPUT_DIR):
	cov_list = []
	n_list = []
	for i in N50_list:
		cov,n50 = N50_list
		cov_list.append(cov)
		n_list.append(n50)

	res = plt.plot(cov_list, n_list)
	plt.xlabel("coverage")
	plt.ylabel("N50 avg")
	plt.title("N50 per cov")
	plt.savefig(os.path.join(output_dir,"N50_fig.png"))
	return



def save_stats(stats,directory_name):
	if not os.path.exists(directory_name):
		os.makedirs(directory_name)
	with open('stats.json', 'w') as fp:
		json.dump(stats, fp)
	return

def simualte_over_coverage(start,end,step,epochs,fastq_files):
	stats_per_cov = []
	N50_list =[]
	for cov in np.arange(start,end,step):
		print("testing for cov = {} ".format(cov))
		N50 = 0
		stats_list = []
		for i in range(epochs):
			stats = single_ineration_per_corr(fastq_files,cov)
			print("###################STATS##########################")
			print(stats)
			N50 += stats["Scaffold Stats"]["N50"]
			stats_list.append(stats)
		N50_list.append((cov,N50/epochs))
		stats_list.append((cov,stats_list))
	generate_plots(stats_list,N50_list)
	# save stats_list
	save_stats(stats,"EXP_1")
	print("FINISHED :)")
	return

#






if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description='Generate new converage dataset',
									 formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument("-1",'--file1', type=str, metavar='<fastq1>',
						help='Specify fastq files for recalculation, comma seperated.',
						required=True)

	parser.add_argument("-2",'--file2', type=str, metavar='<fastq2>',
						help='Specify fastq files for recalculation, comma seperated.',
						required=True)

	parser.add_argument('--coverage', type=float, metavar='<coverage>',
						help='Specify new coverage.',
						required=True)
	args = parser.parse_args()
	fastq_files = [args.file1, args.file2]



	# fastq_files = ["data/tiny_frag_1.fastq","data/tiny_frag_1.fastq"]


	simualte_over_coverage(0.2,1.0,0.2,1,fastq_files)


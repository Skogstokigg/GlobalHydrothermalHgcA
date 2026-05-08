#! /usr/bin/env python3

__version__ = "0.1"
__author__ = "Ester Maria Lopez Garcia"
__email__ = "estermlopezgarcia@gmail.com"


import subprocess # to import bash functions
import sys # from the system
import argparse # to write command-line interfaces
import glob # to retrieve files/pathnames matching a specified pattern
import re # regular expression
import os # get data from operating system
from collections import defaultdict # save a dicc into another dicc
import pandas as pd # data analysis
from multiprocessing import Pool
from datetime import datetime


### --------------------------------------------------------------------- ###
###                               ARGUMENTS	                              ###
### --------------------------------------------------------------------- ###

# Parser options
parser = argparse.ArgumentParser(
	prog="Monkseal",
	description="",
	formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="Written by Ester-Maria Lopez Garcia, 2026. ICM-CSIC"
)
parser.add_argument("-l", "--list_sra", help="include complete path to the list TXT file with SRA codes", type=str, metavar='')
parser.add_argument("-sra", "--get_sra", help="Download SRA files from NCBI. Provide a list [-l] with the SRA codes", action="store_true")
parser.add_argument("-mdata", "--get_metadata", help="Download metadata from NCBI and generate a CSV file", action="store_true")
parser.add_argument("-rd", "--reads_directory", help="Complete path to the directory where reads are located", type=str, metavar='') # By default, it is Marky-Coco's directory", type=str, default="/home/elopez/elopez/phd/bin/marky-coco/", metavar='')#, required=True)
parser.add_argument("-md", "--marky_coco_directory", help="Complete path to the directory where Marky-Coco program is located", type=str, metavar='') #default="/mnt/lustre/scratch/elopez/HV/metags_marky/marky-coco/" default="/home/elopez/elopez/phd/bin/marky-coco/",
parser.add_argument("-otherformat", "--no_ncbi_format", help="Indicate in case of having other than the SRA format as a metagenomic name", action="store_true")
parser.add_argument("-p", "--pair_end", help="Select pair-end reads (R1.fastq and R2.fastq)", action="store_true")
parser.add_argument("-s", "--single_end", help="Select single-end reads (R.fastq)", action="store_true")
parser.add_argument("-c", "--copy_to_marky", help="Copy the selected sequences into the Marky-Coco directory", action="store_true")
parser.add_argument("-cod", "--copy_in_output_directory", help="Copy the selected sequences back into the working directory (reads directory)", type=str, metavar='')
parser.add_argument("-e", "--erase_fastqs", help="Erase SRA and cleaned fastq", action="store_true")
parser.add_argument("-msa", "--marky_coco_single_assembly", help="Analyse the selected sequences with Marky-Coco (single assembly mode)", action="store_true")
parser.add_argument("-a", "--analyse", help="Analyse the output obtained from Marky-Coco", action="store_true")
parser.add_argument("-aout", "--analyses_output", help="Indicate the directory where you want to save the analyses output", type=str, metavar='')
parser.add_argument("-v", "--run_verifier", help="Indicate the path to Verifier program to check merA and merB genes. MerA and MerB reference protein sequences in FASTA format are needed in the folder", type=str, metavar='') 
parser.add_argument("-refmera", "--verifier_ref_merA", help="Indicate the path to the folder where MerA reference protein sequences in FASTA format is", type=str, metavar='') # default="/home/elopez/elopez/phd/data/kolumbo/04_Verifier/merA_bacillus_sp_RC607.fasta"
parser.add_argument("-refmerb", "--verifier_ref_merB", help="Indicate the path to the folder where MerB reference protein sequences in FASTA format is", type=str, metavar='') # default="/home/elopez/elopez/phd/data/kolumbo/04_Verifier/merB_E_coli_plasmid.fasta"
parser.add_argument("-vout", "--verifier_output", help="Indicate the directory where you want to save the verifier result output", type=str, metavar='') # default="/home/elopez/elopez/phd/data/kolumbo/04_Verifier/"
parser.add_argument("-emera", "--verifier_essential_aa_mera", help="List of positions of the essential aa, based on the reference protein (e.g. '[14, 20, 48]')", type=str, metavar='')
parser.add_argument("-emerb", "--verifier_essential_aa_merb", help="List of positions of the essential aa, based on the reference protein (e.g. '[14, 20, 48]')", type=str, metavar='')

parser.add_argument("-ntd", "--nucleotides", help="Indicate de file where saving the nucleotidic sequences of the contig from which hgcAB, merAB sequences were obtained", type=str, metavar='')

args = parser.parse_args()
if not args:
	parser.print_help()
	exit()


### --------------------------------------------------------------------- ###
###                               PROGRAM	                              ###
### --------------------------------------------------------------------- ###

assert sys.version_info >= (3, 7), "This program needs Python 3.7 or higher to run." # Check that the Python version is correct

def main():

	if args.get_sra:
		if not args.list_sra:
			print("Error: [-l] [--list_sra] must be indicated to download SRA files from NCBI")
			exit
		# Download all SRA files from the list
		download()

	if args.get_metadata:
		if not args.list_sra:
			print("Error: [-l] [--list_sra] must be indicated to download metadata associated to SRA files from NCBI")
			exit
		# Download metadata associated to SRA fro the list
		metadata_table = metadata()

	if args.reads_directory: 
		# see how many metagenomes are there downloaded. Grab all FASTQ files and collapse the number
		metag_list_all = []
		metag_string_all = ''
		metag_string_formated_all = ''
		s = 0
		metags = glob.glob(f'{args.reads_directory}*.fastq')

		for metag in metags:
			metag = metag.split("/")[-1] # Select only the name of the file from the path
			obj_complete_code = re.match(r"^(?!.*_cleaned)(.+?)(?:(?:_|_R)([1|2]))?\.fastq$", metag) # Select the code, not the end
			# print(obj_complete_code)
			if obj_complete_code:
				metag = obj_complete_code.group(1)
				# print("metag_no_cleaned: ", metag)
				if not metag in metag_list_all:
					metag_list_all.append(metag)
					metag_string_all += metag + "\n"
					metag_string_formated_all += metag + " "
					s += 1

		print(f'\nThere are {s} metagenomes in total:')
		print(metag_string_all)

	#if args.list:
	if args.single_end:
		if not args.reads_directory:
			print("\nError: [-rd] [--reads_directory] must be indicated to run the analyses")
			exit
		list_sra,flag = select_single_end()

	if args.pair_end:
		if not args.reads_directory:
			print("\nError: [-rd] [--reads_directory] must be indicated to run the analyses")
			exit
		list_sra,flag = select_pair_end()

	if args.copy_to_marky:
		if not args.list_sra:
			print("\nError: [-l] [--list_sra] must be indicated to run the analyses")
			exit
		if not (args.single_end or args.pair_end):
			print("\nError: [-p] [--pair_end] or [-s] [--single_end] must be indicated to run the analyses")
			exit
		copy_to_marky(list_sra,flag)

	if args.marky_coco_single_assembly:
		start = datetime.now()
		result = []
		parallel_function(list_sra, flag, result)
		print("End Time Apply Async:", (datetime.now() - start).total_seconds())

	if args.copy_in_output_directory:
		copy_from_marky(list_sra,flag)

	if args.erase_fastqs:
		erase_fastqs(list_sra,flag)

	if args.analyse:
		if not args.marky_coco_output:
			print("Error: [-mkout] [--marky_coco_output] must be indicated to run the analyses")
			exit
		if not args.marky_coco_directory:
			print("\nError: [-md] [--marky_coco_directory] must be indicated to run the analyses")
			exit
		if not args.analyses_output:
			print("Error: [-aout] [--analyses_output] must be indicated to run the analyses")
			exit
		# Create a dictionary with txid and complete taxon information

		# For each gene
		list_genes = ["hgcA", "hgcB", "merA", "merB"]
		for gene in list_genes:
			if gene == "hgcA":
				# Create master table with as much info as possible
				master_table = make_master_table(gene)


	if args.run_verifier:
		verifier()
	
	if args.nucleotides:
		# Add a column with the nucleotidic sequence given a list of protein_IDs
		retrieve_ntd()



### --------------------------------------------------------------------- ###
###                               FUNCTIONS	                              ###
### --------------------------------------------------------------------- ###


# Parallelise downloading metagenomes
result_down = []
if args.get_sra:
	print(result_down)

def collect_result_down(val):
	return result_down.append(val)

def parallel_function_down(list_metag, result_down):
	result_f = ''
	result_final = []

	with open(list_metag, "r") as IN:
		metags = []
		for metag in IN.readlines():
			metag = metag.strip("\n")
			metags.append(metag)

		pool = Pool(processes=100)
		for metag in metags:
			print(metag)
			result_f = pool.apply_async(download,args=(metag), callback=collect_result_down)
			result_final.append(result_f)

		pool.close()
		for f_res in result_final:
			# r = f_res.get(timeout=10)
			r = f_res.get()
			print(r)


# def download(sra_code):

# 	## PARALLEL
# 	sra_code = sra_code.strip("\n")
# 	print(f'Downloading {sra_code} metagenome from NCBI...')
# 	prefetch = f'~/elopez/phd/bin/sratoolkit.3.0.0-ubuntu64/bin/prefetch {sra_code} --max-size u --output-directory .'
# 	subprocess.call(prefetch, shell=True)
# 	fasterq_dump = f'~/elopez/phd/bin/sratoolkit.3.0.0-ubuntu64/bin/fasterq-dump --split-3 ./{sra_code}/{sra_code}.sra'
# 	subprocess.call(fasterq_dump, shell=True)
# 	print('----------\n')
# 	rm_folder = f'rm -r ./{sra_code}'
# 	subprocess.call(rm_folder, shell=True)


# Download metagenomes function 
def download():
	list_file = args.list_sra
	n_total = ''

	# Count the number of SRA codes that are in the input TXT file
	with open(list_file, "r") as IN:
		n_total = len(IN.readlines())

	with open(list_file, "r") as IN:
		n = 0

		for sra_code in IN.readlines():
			sra_code = sra_code.strip("\n")
			# single = f"{sra_code}.fasta"
			# paired = f"{sra_code}_2.fasta"
			# if (single or paired) not in os.path.isdir(f'{args.marky_coco_directory}'):
			n += 1
			print(f'Downloading {sra_code} metagenome from NCBI... ({n}/{n_total})')
			prefetch = f'/mnt/smart/users/elopez/sratoolkit.3.1.1-ubuntu64/bin/prefetch {sra_code} --max-size u --output-directory .'
			subprocess.call(prefetch, shell=True)
			fasterq_dump = f'/mnt/smart/users/elopez/sratoolkit.3.1.1-ubuntu64/bin/fasterq-dump --split-3 ./{sra_code}/{sra_code}.sra'
			subprocess.call(fasterq_dump, shell=True)
			print('----------\n')
			rm_folder = f'rm -r ./{sra_code}'
			subprocess.call(rm_folder, shell=True)
		print(f'DONE')


# Download metadata function
def metadata():
	list_file = args.list_sra

	with open(list_file, "r") as IN:
		string = ''
		for sra_code in IN.readlines():
			sra_code = sra_code.strip("\n") # remove change of line
			string += sra_code + " "
		print(string)
		# call Kingfisher
		kingfisher = f'~/elopez/phd/bin/kingfisher-download/bin/kingfisher annotate -r {string} -a -f csv > metadata_table.csv'
		subprocess.call(kingfisher, shell=True)


# Zip/unzip function
def zip_and_unzip():
	state_zipped = glob.glob(f'{args.reads_directory}/*.gz')
	if state_zipped:
		print(f'Unziping metagenomes before analysing them to Marky-Coco...')
		unzip = f'gzip -d *gz'
		subprocess(unzip, shell=True)
	state_unzipped = glob.glob(f'{args.reads_directory}/*.fastq')
	if state_unzipped:
		print(f'Ziping metagenomes after having analysed them with Marky-Coco...')
		unzip = f'gzip *fastq'
		subprocess(unzip, shell=True)


# Copy to Marky-Coco folder function
def copy_to_marky(list,flag):
	print(f'Copying {len(list)} metagenome(s) into the Marky-Coco directory...')

	if flag == "p":
		for sra_code in list:
			copy = f'cp {args.reads_directory}{sra_code}*1.fastq {args.reads_directory}{sra_code}*2.fastq {args.marky_coco_directory}'
			subprocess.call(copy, shell=True)
		print("DONE")

	if flag == "s":
		for sra_code in list:
			command = f'cp {args.reads_directory}{sra_code}.fastq {args.marky_coco_directory}'
			subprocess.call(command, shell=True)
		print("DONE")


# Parallelise Marky-Coco analysis
result = []
if (args.pair_end or args.single_end) and args.marky_coco_single_assembly:
	print(result)

def collect_result(val):
	return result.append(val)


def send_to_marky(metag_code,flag):
	print("Calling Marky-Coco...")
	print(f"Start process {metag_code}")
	# time.sleep(3)

	if flag == "p":
		# for sra_code in list:
			# print(sra_code)
		if os.path.exists(f'{args.marky_coco_directory}{metag_code}_1.fastq'):
			marky = f'cd {args.marky_coco_directory} ; rm -fr .snakemake/ ; bash {args.marky_coco_directory}marky_pe.sh {metag_code}' 
		if os.path.exists(f'{args.marky_coco_directory}{metag_code}_R1.fastq'):
			print("Esto no tendría que salir en pantalla")
			print(f'{args.marky_coco_directory}{metag_code}_R1.fastq')
			print(metag_code)
			marky = f'cd {args.marky_coco_directory} ; rm -fr .snakemake/ ; bash {args.marky_coco_directory}marky_pe.sh {metag_code}_R' 
		subprocess.call(marky, shell=True)
		print("DONE")

	if flag == "s":
		# for sra_code in list:
		marky = f'cd {args.marky_coco_directory} ; rm -fr .snakemake/ ; bash {args.marky_coco_directory}marky_se.sh {metag_code}'
		subprocess.call(marky, shell=True)
		print("DONE")

	print(f"End process {metag_code}")
	return marky


def parallel_function(list_metag,flag,result):
	result_f = ''
	result_final = []

	pool = Pool(processes=20)
	for metag in list_metag:
		result_f = pool.apply_async(send_to_marky,args=(metag,flag), callback=collect_result)
		result_final.append(result_f)

	pool.close()
	for f_res in result_final:
		# r = f_res.get(timeout=10)
		r = f_res.get()
		print(r)



# Erase metagenomes that have already been analysed
def erase_fastqs(list,flag):
	print(f'Erasing {len(list)} metagenome(s) FASTQ files from the Marky-Coco directory and from the tmp/ output folder')
	if flag == "p":
		if args.copy_in_output_directory:
			for sra_code in list:
				remove = f'rm {args.marky_coco_directory}{sra_code}_1.fastq {args.marky_coco_directory}{sra_code}_2.fastq \
						{args.copy_in_output_directory}Marky_Coco_OUTPUT/{sra_code}_tmp/{sra_code}_P1.fastq {args.copy_in_output_directory}Marky_Coco_OUTPUT/{sra_code}_tmp/{sra_code}_P2.fastq' # \
						# {args.copy_in_output_directory}Marky_Coco_OUTPUT/{sra_code}_tmp/{sra_code}.sam'
				subprocess.call(remove, shell=True)
		if args.marky_coco_output:
			for sra_code in list:
				remove = f'rm {args.marky_coco_directory}{sra_code}_1.fastq {args.marky_coco_directory}{sra_code}_2.fastq \
						{args.marky_coco_output}{sra_code}_tmp/{sra_code}_P1.fastq {args.marky_coco_output}{sra_code}_tmp/{sra_code}_P2.fastq' # \
						# {args.marky_coco_output}{sra_code}_tmp/{sra_code}.sam'
				subprocess.call(remove, shell=True)
	if flag == "s":
		if args.copy_in_output_directory:
			for sra_code in list:
				remove = f'rm {args.marky_coco_directory}{sra_code}_cleaned.fastq' # \
						# {args.copy_in_output_directory}{sra_code}_tmp/{sra_code}*.fastq \
						# {args.copy_in_output_directory}{sra_code}_tmp/{sra_code}*.sam' # Adapted for the new version of Marky-Coco, where single end metagenomes are first cleaned and renamed as "[SRAcode}_cleaned.fastq"
				subprocess.call(remove, shell=True)
		if args.marky_coco_output:
				remove = f'rm {args.marky_coco_directory}{sra_code}_cleaned.fastq' # \
						# {args.marky_coco_output}{sra_code}_tmp/{sra_code}*.fastq \
						# {args.marky_coco_output}{sra_code}_tmp/{sra_code}.sam' # Adapted for the new version of Marky-Coco, where single end metagenomes are first cleaned and renamed as "[SRAcode}_cleaned.fastq"
				subprocess.call(remove, shell=True)

# Copy output to another folder
def copy_from_marky(list,flag):
	if not os.path.isdir(f'{args.copy_in_output_directory}Marky_Coco_OUTPUT'):
		create_dir = f'cd {args.copy_in_output_directory} ; mkdir Marky_Coco_OUTPUT ; cd Marky_Coco_OUTPUT/'
		subprocess.call(create_dir, shell=True)

	# print(f'Copying {len(list)} metagenome(s) from the Marky-Coco directory into {args.copy_in_output_directory}Marky_Coco_OUTPUT...')ç
	print(f'Copying {len(list)} metagenome(s) from the Marky-Coco directory into {args.copy_in_output_directory}...')

	if flag == "p":
		# Just to have a control of which metagenomes were moved to the other folter
		list_file = f"{args.copy_in_output_directory}paired-end_list_generated_from_directory.txt"
		with open(list_file, "w") as IN:
			for sra_code in list:
				IN.write(sra_code, "\n")
				save_as_out = f'mv {args.marky_coco_directory}{sra_code}_outputs {args.marky_coco_directory}{sra_code}_tmp {args.copy_in_output_directory}Marky_Coco_OUTPUT/'
				subprocess.call(save_as_out, shell=True)
			print("DONE")

	if flag == "s":
		# Just to have a control of which metagenomes were moved to the other folter
		list_file = f"{args.copy_in_output_directory}single-end_list_generated_from_directory.txt"
		with open(list_file, "w") as IN2:
			for sra_code in list:
				print(sra_code)
				IN2.write(f"{sra_code}\n")
				save_as_out = f'mv {args.marky_coco_directory}{sra_code}_outputs {args.marky_coco_directory}{sra_code}_tmp {args.copy_in_output_directory}Marky_Coco_OUTPUT/'
				subprocess.call(save_as_out, shell=True)
			print("DONE")


# Choose PAIR-ENDS
def select_pair_end():
	codes_list = []
	# Read files from the folder
	forward_reads = glob.glob(f'{args.reads_directory}*[_|_R]1.fastq')
	reverse_reads = glob.glob(f'{args.reads_directory}*[_|_R]2.fastq')

	if args.reads_directory:
		for reverse_read in reverse_reads:
			reverse_read = reverse_read.split("/")[-1] # Select only the name of the file from the path
			# print(reverse_read)
			obj_complete_code = re.match(r"(.+?)(?:_|_R)2\.fastq$", reverse_read)
			# print(obj_complete_code)
			complete_code = obj_complete_code.group(1)
			codes_list.append(complete_code)
			# 	print(f'Metagenome code without "[1|2].fastq": {complete_code}')
		print(f'\nThere is/are {len(codes_list)} pair-end sequence(s) in total (R1/R2.fastq):')
		print("\n",codes_list)
		print(f'\nChoose [-msa] [--marky_coco_single assembly] or [-mco] [--marky_coco_coassembly] to analyse them with Marky-Coco and don\'t forget to indicate [-md] [--marky_coco_directory]\n')

	# send to marky
	return(codes_list,"p")


# Choose SINGLE-ENDS
def select_single_end():
	n = 0
	j = 0
	codes_list = []
	single_list = []
	single_paired_list = []

	# Read all single-end files from the folder
	single_reads = [f for f in os.listdir(f'{args.reads_directory}') if re.search(r'^(?!.*_cleaned)[^_]+(?:_(?!1\.fastq$|2\.fastq$)[^_]+)?\.fastq$', f)] # re.search("^[^_]*\.fastq$", f)]
	for single_read in single_reads:
		# Grab all numbers
		obj_complete_code = re.search(r'^([^_]+)\.fastq$', single_read) #re.search("^[a-zA-Z0-9]+(?=\.)", single_read)

		complete_code = obj_complete_code.group(1)
		codes_list.append(complete_code)
		j += 1
		single_list.append(single_read)
		print(single_read, complete_code)

		# Don't choose the single-end files which have paired-ends
		paired_reads = [f for f in os.listdir(f'{args.reads_directory}') if re.search(f"{complete_code}[_|_R]2.fastq", f)]
		if paired_reads:
			single_paired_list.append(paired_reads)
			n += 1

	print(f'There is/are {len(codes_list)} single-end sequence(s) in total, from which {n} have/has also paired-end sequences.')
	print(f'Total single-end: {single_list}\nAlso as paired-end: {single_paired_list}')
	print(f'\nChoose [-msa] [--marky_coco_single assembly] or [-mco] [--marky_coco_coassembly] to analyse them with Marky-Coco and don\'t forget to indicate [-md] [--marky_coco_directory]\n')

	# send to marky
	return(codes_list,"s")


# FUNCTIONS TO ANALYSE THE OUTPUT OF MARKY-COCO
# Generate a table with the relevant information
def plotting_tables(gene):
	# Create dictionaries to save info from IN
	norm_dict = defaultdict(dict) # key: metagenome, value: (key: taxa2, value: norm_abundance)
	stand_dict = defaultdict(dict) # key: metagenome, value: (key: taxa2, value: stand_abundance)

	# Output name
	table_normalised = f"{args.analyses_output}Table_TxidxMetag_Normalised_{gene}.tsv"
	table_standarised = f"{args.analyses_output}Table_TxidxMetag_Standarised_{gene}.tsv"

	# name specifications
	if "hgc" in gene:
		sufix = "final"
	if "mer" in gene:
		sufix = "homologs"

	with open(table_normalised, "w") as NORM_OUT:
		with open(table_standarised, "w") as STAND_OUT:

			master_table = f"{args.analyses_output}Master_Table_4.tsv"

			eof = False
			with open(master_table, "r") as IN:
				header = IN.readline()
				while not eof:
					line = IN.readline()
					if not line:
						eof = True
						break
					# Retrieve fields of interest
					metag = line.split("\t")[0]
					txid = line.split("\t")[5]
					taxa2 = line.split("\t")[14][:-1]
					# NORMALIZATION
					norm_rpoBb = line.split("\t")[12]
					# norm_rpoBa =
					# norm_rpoB =
					# STANDARIZATION
					# stand = line.split("\t")[]

					# Add info to dicts
					norm_dict[metag][taxa2] = norm_rpoBb
			print()
			# print(norm_dict)

			# Write in the file
			out_norm_dict = pd.DataFrame.from_dict(norm_dict)
			out_norm_dict.to_csv(NORM_OUT, sep="\t", na_rep="0")
			# out_stand_dict = pd.DataFrame.from_dict(stand_dict)
			# out_stand_dict.to_csv(STAND_OUT, sep="\t", na_rep="0")

			return(NORM_OUT, STAND_OUT)


def search_motif(gene,sequence,motive):
	if gene == "hgcA":
		if motive != "NA":
			# print(motive)
			if motive in sequence:
				# print("TRUE!!!", gene,motive, sequence)
				return(motive, "TRUE")
			else:
				return("","")
		if motive == "NA":
			# Search a similar motive (conserved Cys)
			# pseudo_motive = "NA"
			# Asp (N) and Glu (Q) may be substitutable
			obj_complete_code = re.search(r'[NQ][a-zA-Z]WC[a-zA-Z][a-zA-Z][a-zA-Z][a-zA-Z][a-zA-Z]', sequence) 
			if obj_complete_code:
				pseudo_motive = obj_complete_code.group(0)
				print(f"It looks like there's a pseudo-motive in here!!! -> {pseudo_motive}")
				return(pseudo_motive, "FALSE")
			else:
				#print("NO pseudo-motive... then putative. We will filter by e-value in Master_Table_6.tsv...")
				return("NA", "FALSE")

	if gene == "hgcB":
		motives = ["CMECGA", "CIECGA"]
		for motive in motives:
			if motive in sequence:
				return(motive,"TRUE")
			else:
				return("","FALSE")


def make_master_table(gene):
	# Output name
	outname = f"{args.analyses_output}Master_Table_1.tsv"

	# Name specifications
	if "hgc" in gene:
		sufix = "final"
	if "mer" in gene:
		sufix = "homologs"

	with open(outname, "w") as OUT:

		# Header line
		OUT.write(f"Metagenome\tGene_ID\tLength\tRead\tCoverage\ttxid\thgcA_Status_Motif\thgcA_Motif\thgcA_AA_Sequence\n")

		# INFILE 1: metag_hgcA_final.txt
		print("\n#### --------------------------------- MOTIVES INFO --------------------------------- ####")
		gene_final_files = glob.glob(f'{args.marky_coco_output}**/*{gene}_{sufix}.txt', recursive=True)
		for gene_final_file in gene_final_files:
			inname = gene_final_file.split('/')[-1]
			print(f'Analysing {inname}...')
			# Flags
			eof = False
			flag = True

			with open(gene_final_file, "r") as IN:
				# Numer of TRUE hgcA
				n = 0
				m = 0

				# Motives
				motives = ["NVWCAAGK","NVWCASGK","NVWCAGGK","NIWCAAGK","NIWCAGGK","NVWCSAGK"]

				if "hgc" in gene:
					header = IN.readline() # avoid 1st line (only hgcA and hgcB have header line, not merA nor merB)
					if "\t" in header:
						print(f"Some issues ocurred in the file format of {inname} that now should be solved")
						header = header.replace("\t", " ")
					n_fields_header = len(header.split(' '))
					header = header.split(' ')
				else:
					continue
				# Read all lines of the file
				while not eof:
					line = IN.readline().strip('\n')

					if not line:
						eof = True
						break

					# To solve format issue [Some files have "\t" as spacers, instead of " ". Even there are some with both types of spacers]
					if "\t" in line:  #and flag == True:
						line = line.replace("\t", " ")
						# print(line)
						# flag = False

					# Check and notify when there are some missing fields in the file (usually the sequence)
					n_fields = len(line.split(' '))
					# print(n_fields, line)
					if n_fields != n_fields_header:
						print(f"WARNING! Check file {inname} due to the lack of information (check: gene sequence, coverage). n_fields_shorter_line: {n_fields}; n_fields_header: {n_fields_header}")
						eof = True
						break
					fields = line.split(' ')
					sequence = fields[5]

					# Save final motif and flag
					mot_save = ""
					flag_mot_save = ""

					# Search hgcA STRICT motif
					true_n = 0
					for motive in motives:
						mot_a,flag_mot = search_motif(gene,sequence,motive) # flag_mot can be TRUE or FALSE, depending on if the motives for hgcA are present or not
						# If the sequence has the strict motif
						if flag_mot == "TRUE":
							n += 1
							mot_save = mot_a
							flag_mot_save = flag_mot
							print(mot_a)
							true_n += 1

					# print(true_n)

					# Search hgcA CONSERVED CYS motif after having searched for all the STRICT motives
					if true_n == 0:
						mot_a,flag_mot = search_motif(gene,sequence,"NA")
						# If the sequence has the motif with the conserved Cys
						mot_save = mot_a
						#print(f"flag_mot: {flag_mot}, mot_save: {mot_save}, mot_a: {mot_a}")
						flag_mot_save = "Putative"

					# Count all the putative hgcA genes found by Marky-Coco
					m += 1

					# Write the selected fields in the OUT file
					OUT.write(f"{inname}\t{fields[0]}\t{fields[1]}\t{fields[2]}\t{fields[3]}\t{fields[4]}\t{flag_mot_save}\t{mot_save}\t{sequence}\n")

				print(f"There are {n} TRUE hgcA genes and {m-n} putative hgcA genes in {inname}")
				print("-----")

	add_hmmer_info(gene) # Generate "plus_hmmer_info"
	add_rpob_info(gene) # Generate "plus_rpob_info"
	add_taxa_info(gene) # Generate "plus_taxa_info"
	add_taxonomic_levels(gene) # Add 7 columns
	master_file = add_hgcb_info(gene) # Add hgcB info
	# master_file = retrieve_ntd(gene) # Add nucleotidic sequence given a list of protein_IDs
	return(master_file)


def add_hmmer_info(gene):
	master_table_1 = f"{args.analyses_output}Master_Table_1.tsv"
	master_table_2 = f"{args.analyses_output}Master_Table_2.tsv"
	# print(f"File name: {master_table_1}")

	eof2 = False
	with open(master_table_1, "r") as IN:
		with open(master_table_2, "w") as OUT:

			# Write new header
			header = IN.readline().strip("\n") # and discard
			new_header = header.replace(header, f"{header}\tE-value_{gene}\tBIC_score_{gene}\n")
			#print(new_header)
			OUT.write(new_header)

			# INFILE 2: metag_hgcA_hmmer.out
			hmmerout_files = glob.glob(f'{args.marky_coco_output}**/*{gene}_hmmer.out', recursive=True)
			print("\n\nRetrieving probability information...")

			# With the multidimensional dictionary, we make sure that each gene has its corresponding value despite being present in >1 metagenome
			evalues_dict = defaultdict(dict) # key: metagenome, values: (key: gene; value: e-value)
			bic_dict = defaultdict(dict) # key: metagenome, values: (key: gene; value: BIT score)

			for hmmerout_file in hmmerout_files:
				inname = hmmerout_file.split('/')[-1]
				print(f'From file {inname}')
				general_adapted_obj = re.search("^.*?hgcA", inname)
				general_adapted = general_adapted_obj.group()
				adapted_inname = f"{general_adapted}_final.txt" 
				#print(adapted_inname)
				# Flags
				eof = False
				with open(hmmerout_file	, "r") as IN_hmmer:
					IN_hmmer.readline() # avoid 1st line
					IN_hmmer.readline() # avoid 2nd line
					IN_hmmer.readline() # avoid 3rd line
					# Read all lines of the file
					while not eof:
						line = IN_hmmer.readline().strip("\n")
						line = line.replace(" ", ",")
						if not line or line[0][0] == "#":
							eof = True
							break
						# Separate fields properly
						line_re = re.sub(r',{2,}', ',', line)
						fields = line_re.split(",")
						tab_separated_line = "\t".join(fields)
						#print(tab_separated_line)
						fields_tab = tab_separated_line.split("\t")

						# Get important fields
						gene_name = fields_tab[0]
						e_value = fields_tab[4]
						bit = fields_tab[5]
						# print(f"gene: {gene_name}, e_value_{gene}: {e_value}, BIC_score_{gene}: {bit}")

						# Fill dictionaries
						evalues_dict[adapted_inname][gene_name] = e_value
						bic_dict[adapted_inname][gene_name] = bit

			# print(evalues_dict)
			n=0
			while not eof2:
				n += 1
				line2 = IN.readline().strip("\n")
				# print(f'line {n}: {line2}')
				if not line2:
					eof2 = True
					break
				# Get first field (gene_name)
				metagenome_name2 = line2.split("\t")[0]
				gene_name2 = line2.split("\t")[1]
				# print(f"metagenome name: {metagenome_name2}, gene name: {gene_name2}")
				# print(len(line2.split("\t")))
				if not "hgcA_hom" in gene_name2 and len(line2.split("\t")) == 9 and gene_name2 in evalues_dict[metagenome_name2]:
					# Add the corresponding e-value and BIT score
					# print(f"{line2}\t{evalues_dict[metagenome_name2][gene_name2]}\t{bic_dict[metagenome_name2][gene_name2]}")
					OUT.write(f"{line2}\t{evalues_dict[metagenome_name2][gene_name2]}\t{bic_dict[metagenome_name2][gene_name2]}\n")



	erase = f"rm {args.analyses_output}Master_Table_1.tsv"
	subprocess.call(erase, shell=True)

	return(OUT)

# 2. Normalisation: hgcA cov/summed rpoB cov, hgcB cov/summed rpoB cov, merA cov/summed rpoB cov, merB cov/summed rpoB cov
def add_rpob_info(gene):

	master_table_2 = f"{args.analyses_output}Master_Table_2.tsv"
	master_table_3 = f"{args.analyses_output}Master_Table_3.tsv"

	with open(master_table_3, "w") as OUT:
		with open(master_table_2, "r") as IN:

			# Write new header
			header = IN.readline().strip("\n") # and discard
			# new_header = header.replace(header, f"{header}\trpoBb_Cov_Metag\tNormalised_Abundance_hgcA\n")
			new_header = header.replace(header, f"{header}\trpoBb_Cov_Metag\tNormalised_Abundance_Bacteria\trpoBa_Cov_Metag\tNormalised_Abundance_Archaea\trpoB_Summed_Cov_Metag\tNormalised_Abundance_hgcA\n")
			# print(new_header)
			OUT.write(new_header)

			# INFILE 3: rpoB_final.txt
			rpoBb_files = glob.glob(f'{args.marky_coco_output}**/*rpoBb_final.txt', recursive=True)
			rpoBa_files = glob.glob(f'{args.marky_coco_output}**/*rpoBa_final.txt', recursive=True)
			print("\n\nRetrieving constitutive genes information...")

			# Dictionary
			rpobb_cov_dict = {} # key: metagenome, values: summed rpoBb cov
			rpoba_cov_dict = {} # key: metagenome, values: summed rpoBa cov

			for rpoBb_file in rpoBb_files:
				summed_cov = 0

				inname = rpoBb_file.split('/')[-1]
				print(f'From file {inname}')
				general_adapted_obj = re.search("^.*?(?=rpoBb)", inname) # Search until the match, but don't include it
				general_adapted = general_adapted_obj.group()
				adapted_inname = f"{general_adapted}hgcA_final.txt"
				# Flags
				eof = False
				with open(rpoBb_file, "r") as IN_rpobb:
					IN_rpobb.readline() # avoid 1st line
					# Read all lines of the file
					while not eof:
						line = IN_rpobb.readline().strip("\n")
						if not line:
							eof = True
							break
						# Separate fields properly
						line = line.replace(" ", "\t")
						cov = line.split("\t")[-1]
						if cov != "length":
							summed_cov += float(cov)
				# print(f"summed_cov: {summed_cov}")
				# Write in the dictionary
				rpobb_cov_dict[adapted_inname] = summed_cov

			for rpoBa_file in rpoBa_files:
				summed_cov_2 = 0

				inname_2 = rpoBa_file.split('/')[-1]
				print(f'From file {inname_2}')
				general_adapted_obj_2 = re.search("^.*?(?=rpoBa)", inname_2) # Search until the match, but don't include it
				general_adapted_2 = general_adapted_obj_2.group()
				adapted_inname_2 = f"{general_adapted_2}hgcA_final.txt"
				# Flags
				eof2 = False
				with open(rpoBa_file, "r") as IN_rpoba:
					IN_rpoba.readline() # avoid 1st line
					# Read all lines of the file
					while not eof2:
						line2 = IN_rpoba.readline().strip("\n")
						if not line2:
							eof2 = True
							break
						# Separate fields properly
						line2 = line2.replace(" ", "\t")
						cov2 = line2.split("\t")[-1]
						if cov2 != "length":
							summed_cov_2 += float(cov2)
				# print(f"summed_cov: {summed_cov}")
				# Write in the dictionary
				rpoba_cov_dict[adapted_inname_2] = summed_cov_2

			# print(rpobb_cov_dict)
			# print(rpoba_cov_dict)

			eof3 = False
			while not eof3:
				line3 = IN.readline().strip("\n")
				# print(f'line {n}: {line2}')
				if not line3:
					eof3 = True
					break
				# Get first field (gene_name)
				metagenome_name = line3.split("\t")[0]
				gene_name = line3.split("\t")[1]
				hgc_cov = line3.split("\t")[4]
				# print(f"metagenome name: {metagenome_name2}, gene name: {gene_name2}, rpoB cov: {rpob_cov_dict[metagenome_name2]}, hgcA cov: {hgc_cov}")
				if hgc_cov == 0:
					print(f"hgcA = 0 is skiped: {metagenome_name}, {gene_name}, {hgc_cov}")
					continue
				# rpoBb
				if metagenome_name in rpobb_cov_dict.keys(): 
					if rpobb_cov_dict[metagenome_name] == 0:
						print(f"rpoBb = 0 is skiped: {metagenome_name}, {gene_name}, {rpobb_cov_dict[metagenome_name]}")
						normalised_rpobb = 0
					if rpobb_cov_dict[metagenome_name] != 0:
						normalised_rpobb = float(hgc_cov) / float(rpobb_cov_dict[metagenome_name])
				# rpoBa
				if metagenome_name in rpoba_cov_dict.keys():
					if rpoba_cov_dict[metagenome_name] == 0:
						print(f"rpoBa = 0 is skiped: {metagenome_name}, {gene_name}, {rpoba_cov_dict[metagenome_name]}")
						normalised_rpoba = 0
					if rpoba_cov_dict[metagenome_name] != 0:
						normalised_rpoba = float(hgc_cov) / float(rpoba_cov_dict[metagenome_name])
				# rpoB summed
				if metagenome_name in rpobb_cov_dict.keys() or metagenome_name in rpoba_cov_dict.keys():
					if rpobb_cov_dict[metagenome_name] != 0 or rpoba_cov_dict[metagenome_name] != 0:
						summed_rpob_coverages = float(rpobb_cov_dict[metagenome_name]) + float(rpoba_cov_dict[metagenome_name])
						normalised_all = float(hgc_cov) / summed_rpob_coverages
						# print(f"{line3}\t{rpobb_cov_dict[metagenome_name]}\t{rpoba_cov_dict[metagenome_name]}\t{summed_rpob_coverages}\t{normalised}")
						# Add the corresponding summed rpoB coverage and the normalised abundance
						OUT.write(f"{line3}\t{rpobb_cov_dict[metagenome_name]}\t{normalised_rpobb}\t{rpoba_cov_dict[metagenome_name]}\t{normalised_rpoba}\t{summed_rpob_coverages}\t{normalised_all}\n")


	erase = f"rm {args.analyses_output}Master_Table_2.tsv"
	subprocess.call(erase, shell=True)

	return(OUT)


def add_taxa_info(gene):

	master_table_3 = f"{args.analyses_output}Master_Table_3.tsv"
	master_table_4 = f"{args.analyses_output}Master_Table_4.tsv"

	with open(master_table_4, "w") as OUT:
		with open(master_table_3, "r") as IN:

			# Write new header
			header = IN.readline().strip("\n") # and discard
			new_header = header.replace(header, f"{header}\tTaxo1\tTaxo2\n")
			#print(new_header)
			OUT.write(new_header)

			# INFILE 4: Open the txID-to-taxon association file
			db_txid = f'{args.marky_coco_directory}db/db_txid_220220.txt' 
			print("\n\nRetrieving taxonomical information...")

			# Dictionary
			taxo1_dict = {} # key: txid, values: taxo1
			taxo2_dict = {} # key: txid, values: taxo2

			# Flags
			eof = False
			with open(db_txid, "r") as IN_txid:
				# Read all lines of the file
				while not eof:
					line = IN_txid.readline()
					if not line:
						eof = True
						break
					# Split rows by space
					txid = line.split('\t')[0]
					taxo1 = line.split('\t')[1]
					taxo2 = line.split('\t')[2][:-1]
					# Fill the dictionary
					taxo1_dict[txid] = taxo1
					taxo2_dict[txid] = taxo2

			eof2 = False
			n = 0
			# info_txids = {}
			while not eof2:
				line2 = IN.readline().strip("\n")
				# print(f'line {n}: {line2}')
				if not line2:
					eof2 = True
					break
				# Get fields
				metag = line2.split("\t")[0]
				gene2 = line2.split("\t")[1]
				txid2 = line2.split("\t")[5]

				if txid2 == "2629452" or txid2 == "2638683" or txid2 == "2220" or txid2 == "67817" or txid2 == "2206" or txid2 not in taxo1_dict.keys():
					n += 1
					print(f"Txid {txid2} from {metag} ({gene2}) wasn't in the taxonomy file from the Marky-Coco DB")
					OUT.write(f"{line2}\tNA\tNA\n")
					continue
				else:
					#print(f"{metag} ({gene2}): {taxo1_dict[txid2]}\t{taxo2_dict[txid2]}")
					# Add the taxonomy
					OUT.write(f"{line2}\t{taxo1_dict[txid2]}\t{taxo2_dict[txid2]}\n")

			# print(f"The following {n} txids did not be found in the DB:\n {}")

	erase = f"rm {args.analyses_output}Master_Table_3.tsv"
	subprocess.call(erase, shell=True)

	return(OUT)

def file_exists_partial_name(partial_name):
	# print("foca monje")
	print(os.listdir(args.analyses_output))
    # for filename in os.listdir(f'{args.analyses_output}'):
    #     if partial_name in filename:
    #         print("True")
	# 		return True

	# print("False")
	# return False

def add_taxonomic_levels(gene):

	# Retrieve the most updated taxonomic information from NCBI with ncbitax2lin
	folder_files = glob.glob(f'{args.analyses_output}*')
	lineages_complete_name = ''
	flag = False
	for folder_file in folder_files:
		if "ncbi_lineages" and ".csv" in folder_file:
			lineages_complete_name = folder_file
			print(f"The lineage folder exists: {lineages_complete_name}")
			flag = True
	if not flag:
		print("Installing ncbitax2lin to retrieve the complete lineage associated to the txID code...")
		# Install ncbtax2lin
		install_ncbitax2lin = f"pip install --upgrade pip; pip install -U ncbitax2lin"
		subprocess.call(install_ncbitax2lin, shell=True)
		# Generate linage file
		get_linage = f"wget -N ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz; mkdir -p taxdump && tar zxf taxdump.tar.gz -C ./taxdump; ncbitax2lin --nodes-file taxdump/nodes.dmp --names-file taxdump/names.dmp; gzip -d ncbi_lineages_*.csv.gz; mv ncbi_lineages_*.csv {args.analyses_output}"
		subprocess.call(get_linage, shell=True)
		for folder_file in folder_files:
			if "ncbi_lineages" and ".csv" in folder_file:
				lineages_complete_name = folder_file
				print(f"The lineage folder NOW exists: {lineages_complete_name}")
				flag = True

	# Add info to the Master Table
	master_table_4 = f"{args.analyses_output}Master_Table_4.tsv"
	master_table_5 = f"{args.analyses_output}Master_Table_5.tsv"

	with open(master_table_5, "w") as OUT:
		with open(master_table_4, "r") as IN:

			# Write new header
			header = IN.readline().strip("\n") # and discard
			new_header = header.replace(header, f"{header}\tSuperkingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n")
			#print(new_header)
			OUT.write(new_header)

			# Add 7 levels of taxonomy
			# INFILE 5: Open the ncbitax2lin file
			in_lineages_complete_name = f"{lineages_complete_name}"

			# txid_lineages = f"{args.analyses_output}ncbi_lineages_*.csv"
			print("\n\nAssigning lineages...")

			# Dictionary
			superkingdom_dict = {} # key: txid, values: superkingdom
			phylum_dict = {} # key: txid, values: phylum
			class_dict = {} # key: txid, values: class
			order_dict = {} # key: txid, values: order
			family_dict = {} # key: txid, values: family
			genus_dict = {} # key: txid, values: genus
			species_dict = {} # key: txid, values: species

			# Flags
			eof = False
			with open(in_lineages_complete_name, "r") as IN_txid:
				# Read all lines of the file
				while not eof:
					line = IN_txid.readline()
					if not line:
						eof = True
						break
					# Split rows by comma
					txid = line.split(',')[0]
					superkingdom = line.split(',')[1]
					phylum = line.split(',')[2]
					class_ = line.split(',')[3]
					order = line.split(',')[4]
					family = line.split(',')[5]
					genus = line.split(',')[6]
					species = line.split(',')[7]
					# Fill the dictionary
					superkingdom_dict[txid] = superkingdom
					phylum_dict[txid] = phylum
					class_dict[txid] = class_
					order_dict[txid] = order
					family_dict[txid] = family
					genus_dict[txid] = genus
					species_dict[txid] = species

			# print(superkingdom_dict)
			# print(phylum_dict)

			eof2 = False
			n = 0
			# info_txids = {}
			while not eof2:
				line2 = IN.readline().strip("\n")
				# print(f'line {n}: {line2}')
				if not line2:
					eof2 = True
					break
				# Get fields
				metag = line2.split("\t")[0]
				gene2 = line2.split("\t")[1]
				txid2 = line2.split("\t")[5]

				if txid2 not in superkingdom_dict.keys():
					n += 1
					print(f"Txid {txid2} from {metag} ({gene2}) wasn't in the DB")
					OUT.write(f"{line2}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
					continue
				else:
					# Add the lineage
					OUT.write(f"{line2}\t{superkingdom_dict[txid2]}\t{phylum_dict[txid2]}\t{class_dict[txid2]}\t{order_dict[txid2]}\t{family_dict[txid2]}\t{genus_dict[txid2]}\t{species_dict[txid2]}\n")

	erase = f"rm {args.analyses_output}Master_Table_4.tsv"
	subprocess.call(erase, shell=True)

	return(OUT)


def add_hgcb_info(gene):

	master_table_5 = f"{args.analyses_output}Master_Table_5.tsv"
	master_table_6 = f"{args.analyses_output}Master_Table_6.tsv"

	with open(master_table_6, "w") as OUT:
		with open(master_table_5, "r") as IN:

			# Write new header
			header = IN.readline().strip("\n") # and discard
			# new_header = header.replace(header, f"{header}\trpoBb_Cov_Metag\tNormalised_Abundance_hgcA\n")
			new_header = header.replace(header, f"{header}\tHgcB_Presence\tHgcB_Cov\tHgcB_AA_Sequence\n")
			# print(new_header)
			OUT.write(new_header)

			# INFILE 6: sample_hgcB_final.txt
			hgcB_files = glob.glob(f'{args.marky_coco_output}**/*hgcB_final.txt', recursive=True)
			print("\n\nRetrieving hgcB gene information...")

			# Dictionary
			hgcb_dict = defaultdict(dict) # key: metagenome, values: (key: gene; value: hgcB protein)
			cov_dict = defaultdict(dict) # key: metagenome, values: (key: gene; value: hgcB coverage)
			seq_dict = defaultdict(dict) # key: metagenome, values: (key: gene; value: hgcB sequence)

			for hgcB_file in hgcB_files:

				inname = hgcB_file.split('/')[-1]
				print(f'From file {inname}')
				general_adapted_obj = re.search("^.*?(?=hgcB)", inname) # Search until the match, but don't include it
				general_adapted = general_adapted_obj.group()
				adapted_inname = f"{general_adapted}hgcA_final.txt"
				# Flags
				eof = False
				with open(hgcB_file, "r") as IN_hgcb:
					IN_hgcb.readline() # avoid 1st line
					# Read all lines of the file
					while not eof:
						line = IN_hgcb.readline().strip("\n")
						if not line:
							eof = True
							break
						# Separate fields properly
						line = line.replace(" ", "\t")
						line = line.split("\t")
						
						# Get info
						gene_name_protein = line[0]
						gene_name_contig_1 = gene_name_protein.split("_")[0]
						gene_name_contig_2 = gene_name_protein.split("_")[1]
						gene_name_contig = f"{gene_name_contig_1}_{gene_name_contig_2}" 
						
						# print(f"gene_name_contig: {gene_name_contig}, gene_name_protein: {gene_name_protein}")
						cov = line[3]
						seq = line[4]
						# print(f"gene: {gene_name}, cov: {cov}, seq: {seq}")

						# Fill dictionaries
						hgcb_dict[adapted_inname][gene_name_contig] = gene_name_protein
						cov_dict[adapted_inname][gene_name_contig] = cov
						seq_dict[adapted_inname][gene_name_contig] = seq

			eof2 = False
			while not eof2:
				line2 = IN.readline().strip("\n")
				if not line2:
					eof2 = True
					break
				# Get first field (gene_name)
				metagenome_name = line2.split("\t")[0]
				gene_name = line2.split("\t")[1]
				gene_name_contig2_1 = gene_name.split("_")[0]
				gene_name_contig2_2 = gene_name.split("_")[1]
				gene_name_contig2 = f"{gene_name_contig2_1}_{gene_name_contig2_2}" 

				if gene_name_contig2 in hgcb_dict[metagenome_name]:
					# print(f"{line2}\t{hgcb_dict[metagenome_name][gene_name_contig2]}\t{cov_dict[metagenome_name][gene_name_contig2]}\t{seq_dict[metagenome_name][gene_name_contig2]}\n")
					# Add the corresponding hgcB presence and gene_id
					OUT.write(f"{line2}\t{hgcb_dict[metagenome_name][gene_name_contig2]}\t{cov_dict[metagenome_name][gene_name_contig2]}\t{seq_dict[metagenome_name][gene_name_contig2]}\n")
				else:
					OUT.write(f"{line2}\tNA\tNA\tNA\n")


	erase = f"rm {args.analyses_output}Master_Table_5.tsv"
	#subprocess.call(erase, shell=True)

	return(OUT)


def retrieve_ntd():

	master_table_6 = f"{args.nucleotides}" # -ntd Master_Table_6_Cys_Protein_ID.tsv (hgcA); (merB); (merA)
	master_table_ntd = f"{args.analyses_output}Master_Table_Ntd_Paired_End_Smart.tsv"

	with open(master_table_ntd, "w") as OUT:
		with open(master_table_6, "r") as IN:
			
			# Write new header
			header = IN.readline().strip("\n") # and discard
			# new_header = header.replace(header, f"{header}\trpoBb_Cov_Metag\tNormalised_Abundance_hgcA\n")
			new_header = header.replace(header, f"{header}\tContig_Sequence\n")
			# print(new_header)
			OUT.write(new_header)

			# INFILE (7): final.contigs.fa
			contigs_files = glob.glob(f'{args.marky_coco_output}**/final.contigs.fa', recursive=True)
			print("\n\nRetrieving contigs...")

			# Dictionary
			contig_dict = defaultdict(dict) # key: metagenome, values: (key: gene; value: contig nucleotidic sequence)

			for contigs_file in contigs_files:
				
				# print(contigs_file)
				sufix = contigs_file.split('/')[-1]
				metag_name1 = contigs_file.split('/')[-2]
				metag_name = metag_name1.split('_')[-2]
				# print(f"metag_name: {metag_name}")
				new_inname1 = f"{metag_name}_{sufix}"
				part1 = new_inname1.split(".")[-3]
				part2 = new_inname1.split(".")[-2] 
				part3 = new_inname1.split(".")[-1]
				new_inname = f"{part1}_{part2}.{part3}"
				print(f'From file {new_inname}')

				# Flags
				eof = False
				with open(contigs_file, "r") as IN_contigs:

					# Flag
					n_header = 0
					n_seq = 0

					# Read all lines of the file
					while not eof:
						line = IN_contigs.readline().strip("\n")
						if not line:
							eof = True
							break

						# Separate fields properly
						if ">" in line[0]:
							n_header += 1
							contig_id = line.split(" ")[0][1:]
							# print(contig_id)
						else:
							n_seq += 1
							contig_seq = line

						if n_header == n_seq:
							# print(n_header, n_seq)
							# print(f"{contig_id}: {contig_seq}")
							contig_dict[metag_name][contig_id] = contig_seq

			# print(contig_dict)			
			
			eof2 = False
			while not eof2:
				line2 = IN.readline().strip("\n")
				if not line2:
					eof2 = True
					break
				# Get first field (gene_name)
				metagenome_name = line2.split("\t")[0]
				gene_name = line2.split("\t")[1]
				gene_name_contig2_1 = gene_name.split("_")[0]
				gene_name_contig2_2 = gene_name.split("_")[1]
				gene_name_contig2 = f"{gene_name_contig2_1}_{gene_name_contig2_2}" 
				# print(metagenome_name, gene_name, gene_name_contig2)

				# print(gene_name_contig2, contig_dict[metagenome_name])
				
				if gene_name_contig2 in contig_dict[metagenome_name]:
					# print(f"{line2}\t{contig_dict[metagenome_name][gene_name_contig2]}\n")
			 		# Add the corresponding contig
					print("funsiona")
					OUT.write(f"{line2}\t{contig_dict[metagenome_name][gene_name_contig2]}\n")
				else:
					OUT.write(f"{line2}\tNA\n")


	#erase = f"rm {args.analyses_output}Master_Table_5.tsv"
	#subprocess.call(erase, shell=True)

	#return(OUT)	


# def make_master_table_merA(gene):



# 3. Filter sequences to match only those with the pattern
def hgcab_sequence_collocated():

	file_genes_id = sys.argv[1] # gene_id hgcA (or the smallest from the "hgcA_final" and "hgcB_final" files)
	fileA = sys.argv[2]
	fileB = sys.argv[3]

	gene_id_list = [] # gene_id
	with open(file_genes_id, "r") as in1:
		allfile1 = in1.readlines()[1:]
		for gene_id in allfile1:
			#print(gene_id, end="")
			gene_id = gene_id.strip('\n')
			gene_id_list.append(gene_id)

	genes_ids_2 = []
	genes_ids_3 = []
	with open(fileA, "r") as in2:
		with open(fileB, "r") as in3:
			allfile2 = in2.readlines()[1:]
			allfile3 = in3.readlines()[1:]

			for line2 in allfile2:
				line2 = line2.strip('\n')
				gene_id_2 = line2.split(' ')[0]
				hgc_status_2 = line2.split(' ')[-1]
				genes_ids_2.append(gene_id_2)

			for line3 in allfile3:
				line3 = line3.strip('\n')
				gene_id_3 = line3.split(' ')[0]
				hgc_status_3 = line3.split(' ')[-1]
				genes_ids_3.append(gene_id_3)

	print(len(genes_ids_2))
	print(len(genes_ids_3))
	for elem in gene_id_list:
		for id2 in genes_ids_2:
			if elem in id2:
				print(elem, id2, "hgcA_hom")
		for id3 in genes_ids_3:
			if elem in id3:
				print(elem, id3, "hgcB_hom\tTrue?")



# OPTIONAL, this is not strictly needed for the analyses
def standarise_by_reads(dictionary):

	# Create a dictionary to save normalized data
	normalization_dict = defaultdict(dict)

	outname = "Table_2_Standarisation_Reads.tsv"
	with open(outname, "w") as OUT:

		bowtie2_files = glob.glob(f'{args.marky_coco_directory}/**/*bowtie2.log')
		for bowtie2_file in bowtie2_files:

			inname = bowtie2_file.split('/')[-1]
			print(f'Retrieved the total number of reads mapped from {inname}...')
			eof = False

			with open(bowtie2_file, "r") as IN:
				# Read all lines of the file
				while not eof:
					line = IN.readline()
					if not line:
						eof = True
						break

					# Edit dictionary (to normalize values)
					#for key in dictionary:


# 4. Check with verifier merA and merB
def verifier():

	# Convert Marky-Coco merA and merB files to FASTA files
	# merA_final_files = glob.glob(f'{args.marky_coco_output}**/*merA_.txt', recursive=True)
	# merB_final_files = glob.glob(f'{args.marky_coco_output}**/*merB_.txt', recursive=True)
	# awk '{print ">ERR2834528_"$1"\n"$6}' ERR2834528_outputs/ERR2834528_merB_homologs.txt > ../04_Verifier/ERR2834528_merB_homologs.fasta
	# awk '{print ">ERR2834528_"$1"\n"$6}' ERR2834528_outputs/ERR2834528_merA_homologs.txt > ../04_Verifier/ERR2834528_merA_homologs.fasta


	# INFILE 1: metag_mer[A|B]_homologs.fasta
	merAs_to_check = glob.glob(f'{args.verifier_output}/*merA_homologs.fasta', recursive=True)
	merBs_to_check = glob.glob(f'{args.verifier_output}/*merB_homologs.fasta', recursive=True)
	print(merAs_to_check)
	print(merBs_to_check)

	activate_mafft = f"module load mafft"
	subprocess.call(activate_mafft, shell=True)

	for merA_to_check in merAs_to_check:
		print(merA_to_check)
		name = merA_to_check.split("/")[-1]
		name_no_fasta = name.split(".")[0]
		print("hola", name, name_no_fasta)
		check_merA = f"python3 {args.run_verifier} -ref {args.verifier_ref_merA} -t {merA_to_check} -e {args.verifier_essential_aa_mera} -o {args.verifier_output} -n {name_no_fasta}_merAs_verified -TPF"
		subprocess.call(check_merA, shell=True)

	for merB_to_check in merBs_to_check:
		print(merB_to_check)
		name = merB_to_check.split("/")[-1]
		name_no_fasta = name.split(".")[0]
		print(name)
		check_merB = f"python3 {args.run_verifier} -ref {args.verifier_ref_merB} -t {merB_to_check} -e {args.verifier_essential_aa_merb} -o {args.verifier_output} -n {name_no_fasta}_merBs_verified -TPF"
		subprocess.call(check_merB, shell=True)


main()

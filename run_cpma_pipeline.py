import sys
import subprocess
import os
import argparse

def get_edit_freq(input, sample, chr, mutation, output):
	chr = "chr" + chr
	beginning = int(mutation)-500
	end = int(mutation)+500
	file1 = os.path.join(input, sample + "_R1_001.fastq.gz")
	file2 = os.path.join(input, sample + "_R2_001.fastq.gz")
	bam = sample + ".bam"
	bam_mutation = f'{sample}_{chr}_{beginning}_{end}.bam'
	# bam_mutation = sample + "_" + chr + "_" + str(beginning) + "_" + str(end) + ".bam"

	output_path = os.path.join(output, sample)
	if not os.path.isdir(output_path):
		os.mkdir(output_path)
	# subprocess.call(["mkdir", sample])

	# minimap_cmd = f'minimap2 -t 10 -ax sr hg38.mni {file1} {file2}'.split(' ')
	minimap = subprocess.Popen(["minimap2", "-t", "10", "-ax", "sr", "hg38_sr.mni", file1, file2], stdout=subprocess.PIPE)
	sam_view = subprocess.Popen(["samtools", "view", "-b"], stdin=minimap.stdout, stdout=subprocess.PIPE)
	# minimap.stdout.close()

	bam_path = os.path.join(output_path, bam)
	print(bam_path)
	with open(bam_path, "w") as f:
		sam_sort = subprocess.Popen(["samtools", "sort", "-@10"], stdin=sam_view.stdout, stdout=f).wait()
	print("Finished sam_sort")
	# sam_view.stdout.close()

	bam_mutation_path = os.path.join(output_path, bam_mutation)
	print(bam_mutation_path)
	subprocess.call(["samtools", "index", bam_path])
	view_window = chr + ":" + str(beginning) + "-" + str(end)
	with open(bam_mutation_path, "w") as f:
		subprocess.Popen(["samtools", "view", "-b", bam_path, chr, view_window], stdout=f).wait()

	#sam_index_cmd = f'minimap2 -t 10 -ax sr hg38.mni {file1} {file2}'.split(' ')

	readcounts_raw = os.path.join(output_path, sample + "_readcounts_raw.txt")
	readcounts_path = os.path.join(output_path, readcounts_raw)
	with open(readcounts_path, "w") as f:
		subprocess.call(["bam-readcount", "-f", "/storage/cynthiawu/iSeq/hg38.fa", "-w", "1", bam_mutation_path], stdout=f)
		
	readcounts = os.path.join(output_path, sample + "_readcounts.txt")
	subprocess.call(["python", "prepare_readcounts.py", "-i", readcounts_path, "-o", readcounts])

def main():
	parser = argparse.ArgumentParser()

	parser.add_argument('-i', '--input', type=str, dest='input', help='Input folder with sample')
	parser.add_argument('-n', '--name', type=str, dest='sample', help='Sample Name')
	parser.add_argument('-c', '--chr', type=str, dest='chr', help='Chromosome where mutation is')
	parser.add_argument('-m', '--mutation', type=int, dest='mutation', help='Position where mutation is')
	parser.add_argument('-o', '--output', type=str, dest='output', help='Destination folder for output')
	input_values = parser.parse_args()

	get_edit_freq(input_values.input, input_values.sample, input_values.chr, input_values.mutation, input_values.output)


if __name__ == '__main__':
	main()
'''	
bam_read = subprocess.Popen(["bam-readcount", "-f", ""], stdin=minimap.stdout, stdout=subprocess.PIPE)
minimap.stdout.close()
with open(sample/bam, "w") as f:
	sam_sort = subprocess.Popen(["samtools", "sort", "-@10"], stdin=sam_view.stdout, stdout=f)
	
#output = sam_sort.communicate()
\
'''

'''
minimap2 -t 10 -ax sr GRCh38_latest_genomic.mni file1 file2 | samtools view -b | samtools sort -@10 > sample/bam
samtools view -b sample/bam chr,_GRCh38.p12_Primary_Assembly:beginning-end >  sample/bam_mutation
samtools index  sample/bam_mutation


/storage/cynthiawu/soft/bam-readcount/bam-readcount/bin/bam-readcount -f GRCh38_latest_genomic.fna -w 1 sample/bam_mutation >  HEP-GFP-R76W_S2_L001_readcounts_raw.txt
python prepare_readcounts.py -i HEP-GFP-R76W_S2_L001_readcounts_raw.txt -o HEP-GFP-R76W_S2_L001_readcounts.txt
'''

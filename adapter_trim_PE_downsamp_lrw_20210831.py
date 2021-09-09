# adapter trim with reads downsampling

import sys,os
from random import sample

def reverse_complement(seq):
	complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c'}
	bases = list(seq)
	bases = [complement[base] for base in bases]
	return ''.join(bases[::-1])


def clean_fastq_wPrimer(input_name1, input_name2, output_name1):
	from itertools import izip
	import gzip
	with open(input_name1, "rb") as r1, open(input_name2, "rb") as r2, open(output_name1, "wb") as w1:
	#gzip.open(input_name1, "rb") as r1, gzip.open(input_name2, "rb") as r2, gzip.open(output_name1, "wb") as w1:
		counter = 0
		last = False
		bunchmax = 1000000
		bunchsize = 0
		bunchnum = 0
		bunch1 = []
		#bunch2 = []
		minlen = 60

		downsamp_counter = 0
		samplist_counter = 0

		totalreads = 1441570# total reads of the library
		samplereads = 500000 # how many reads to use, if same as total reads, no downsampling

		allreadslist1 = range(1,totalreads+1)
		samplelist = sample(allreadslist1,samplereads)
		samplelist.sort()


		for l1, l2 in izip(r1, r2):
			if counter%4 == 0: #name
				#create two read elements and 
				read1 = ['','','','']
				#read2 = ['','','','']
				read1[0] = l1
				#read2[0] = l2
				downsamp_counter += 1
			elif counter%4 == 1:#seq
				seq1 = l1
				seq2 = l2
				to_write = 1
				if 'N'in l1 or 'N' in l2:
					to_write = 0
				else:
					seq2find = reverse_complement(seq1[:10])
					cut_pos = seq2.find(seq2find)
					if cut_pos == -1:
						to_write = 0
						#read1[1] = seq1
						#read2[1] = seq2
					else:
						new_seq1 = seq1[:cut_pos + 10]+'\n'
						new_seq2 = seq2[:cut_pos + 10]+'\n'

						if not (reverse_complement(new_seq2[:-1]) == new_seq1[:-1]):
							to_write = 0
						else:
							cut_seq1 = new_seq1[19:] #includes \n
							umi_seq1 = new_seq1[4:19] #not includes \n
							read1[1] = cut_seq1
							spacepos = read1[0].find(' ')
							read1[0] = read1[0][:spacepos]+'UMI'+umi_seq1+'\n'
							#read2[1] = new_seq2
			elif counter%4 == 2:#+
				read1[2] = l1
				#read2[2] = l2
			elif counter%4 == 3:#Qscore
				if samplist_counter < samplereads: 
					if downsamp_counter == samplelist[samplist_counter]: #decide whether this read belong to the sample list
						samplist_counter += 1

						if to_write == 1:
							if cut_pos == -1:
								read1[3] = l1
								#read2[3] = l2
							else:
								read1[3] = l1[19:cut_pos + 10]+'\n'
								#read2[3] = l2[:cut_pos + 10]+'\n'

							if cut_pos > (minlen-10) or cut_pos == -1:
								if to_write == 1:
									bunch1.extend(read1)
									#bunch2.extend(read2)
									bunchsize += 1

			if bunchsize > bunchmax:
				w1.writelines(bunch1)
				#w2.writelines(bunch2)
				bunchsize = 0 
				bunchnum += 1
				bunch1 = []
				#bunch2 = []
			counter += 1 
		w1.writelines(bunch1)
		#w2.writelines(bunch2)
		#print '%s\t%s\n' % (counter/4, bunchnum*bunchmax+bunchsize)
		with open('trim_result.txt', 'a') as f:
			f.write('%s\t%s\n' % (counter/4, bunchnum*bunchmax+bunchsize))


clean_fastq_wPrimer('lib1_R1.fastq','lib1_R2.fastq','trim1.fastq')

# lib_start = 3
# lib_stop = 4
# for i in range(lib_start, lib_stop+1):
# 	clean_fastq_wPrimer('../raw_reads/lib%d_R1_001.fastq.gz' % i, 
# 		'../raw_reads/lib%d_R2_001.fastq.gz' % i, 
# 		'../trimmed_reads/lib%d_trim.fastq.gz' % i)
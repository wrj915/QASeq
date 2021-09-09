import pysam
from collections import Counter

def umi2num (umi_str):
	#umi_str = upper(umi_str)
	bases = list(umi_str)
	complement = {'A': 0, 'C': 1, 'G': 2, 'T': 3} 
	bases = [complement[base] for base in bases]
	len_bases = 15
	out_num = 0
	for i in range(15):
		out_num = out_num + 4**(len_bases-i-1)*bases[i]
	return out_num

def num2umi(umi_num):
	complement = { 0:'A', 1:'C', 2:'G', 3:'T'} 
	bases = ['']*15
	for i in range(15):
		base = umi_num%4
		bases[15-i-1] = complement[base]
		umi_num = umi_num/4
	return ''.join(bases)

def write_SNP_file(file_name, umi_list):
	with open(file_name, 'w') as f:
		for umi in umi_list:
			if sum(umi_list[umi].values()) > 2:
				for seq, seq_count in umi_list[umi].most_common():
					f.write("{}\t{}\t{}\n".format(umi, seq, seq_count))
				f.write('\n')
primer_info ={}
primer_info[1-1] = ['CCTTAGACAACTACCTTTCTACGGAC','GTGCAGGGGGCAGACGA']
primer_info[2-1] = ['ATTCCAGTGGCCATCAAAGTGT','AGGGTGGAGGGGCTTACG']
primer_info[3-1] = ['CAGGAAGCATACGTGATGGCT','CCCAGAAGGCGGGAGACATA']
primer_info[4-1] = ['GATGAGCTACCTGGAGGATGTG','GTTCCGAGCGGCCAAGTC']
primer_info[5-1] = ['GATGGCGCTGGAGTCCATT','ACACATCACTCTGGTGGGTGAA']
primer_info[101-1] = ['GGGACCCACTCCATCGAGA','CTCACAGTAAAAATAGGTGATTTTGGTCT']
primer_info[102-1] = ['GGAGTATTTCATGAAACAAATGAATGATGC','AAGATCCAATCCATTTTTGTTGTCCAG']
primer_info[103-1] = ['CAATTTCTACACGAGATCCTCTCTCT','CTGTGACTCCATAGAAAATCTTTCTCC']
primer_info[104-1] = ['GCCGCCAGGTCTTGATGTACT','GTCTGACGGGTAGAGTGTGC']
primer_info[105-1] = ['CTCACCATCGCTATCTGAGCAG','CACATGACGGAGGTTGTGAGG']
primer_info[106-1] = ['TCGTCAAGGCACTCTTGCCTA','CTGAAAATGACTGAATATAAACTTGTGGTAGT']

def analysis(infile_name):
	samfile = pysam.AlignmentFile(infile_name, "rb" )
	#readcount = 0 
	pre_ref_id = 0 
	umi_counter = Counter()
	boo = True
	while boo:
		try:
			#readcount = readcount + 1
			read = samfile.next()
		except:
			boo = False
			break

		ref_id = read.reference_id
		umi = read.query_name[-15:]
		
		if ref_id < 0:
			outfile_name = 'UMIquantlib'+infile_name[3:-4]+'tar'+str(pre_ref_id + 1)+'.txt'
			with open(outfile_name,'w') as f:
				for k,v in umi_counter.most_common():
					f.write("{}\t{}\n".format(k,v))
			boo = False
			break


		if ref_id != pre_ref_id: 
			# write previous count information
			outfile_name = 'UMIquantlib'+infile_name[3:-4]+'tar'+str(pre_ref_id + 1)+'.txt'
			with open(outfile_name,'w') as f:
				for k,v in umi_counter.most_common():
					f.write("{}\t{}\n".format(k,v))
			# reset pre_ref_id
			pre_ref_id = ref_id
			# create another space for counting
			umi_counter = Counter()
			print ref_id

		umi_counter[umi] += 1


lib_start = 1
lib_stop = 1
for i in range(lib_start, lib_stop+1):
	samname = "lib"+str(i)+"_sorted.bam"
	analysis(samname)

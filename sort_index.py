import pysam

lib_start = 1
lib_stop = 1
for i in range(lib_start, lib_stop+1):
	pysam.sort('-o', 'lib%d_sorted.bam' % i, 'alignlib%d.sam' % i)
	pysam.index('lib%d_sorted.bam' % i)
	
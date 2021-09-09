echo trimming...
python adapter_trim_PE_downsamp_lrw_20210831.py

echo building index...
bowtie2-build all2seq_new.fasta all2seq_new

# different panels:
#bowtie2-build all179seq.fasta all179seq
#bowtie2-build all226seq.fasta all226seq

for i in 1
do
	echo lib$i
	bowtie2 -x all2seq_new -U trim$i.fastq -S alignlib$i.sam # may need to change index name
done

echo sorting...
python sort_index.py

echo UMI grouping...
python UMIcounter.py

echo UMI counting...
export PATH=$PATH:/Applications/MATLAB_R2019a.app/bin/
matlab -nodisplay -nosplash -nodesktop -r "run('analysis_umisorted_v2_Dynamic.m');exit;"
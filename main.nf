#!/usr/bin/env nextflow

seqlist='s3://pipe.scratch.3/resources/accessions.txt'
seqdata= Channel.fromPath(seqlist)

chlamyref='s3://pipe.scratch.3/resources/Chlamy23s.fasta'
refseq=Channel.fromPath(chlamyref)


chunks=36

process runfasta {
	
	input:
	each x from 1..chunks
  	path accession from seqdata
	
	output:
	file "*_{1,2}.fastq" into dumpouts
	
	
	"""
  	fil=`head -n $x $accession | tail -n 1`
	fastq-dump --split-3 \$fil
	"""

}


refseq.into{bwachlamy; consensuschlamy}

process bwamap {
	
	input:
  	path sradump from dumpouts
	path chlamy from bwachlamy
	
	output:
	file "*.sorted.bam" into bams
	
	
	"""
	bwa index $chlamy
	bwa mem $chlamy *1.fasta *2.fasta > bwa_mapped.sam
	samtools view -bS bwa_mapped.sam > bwa_mapped.bam
	samtools sort bwa_mapped.bam -o bwa_mapped.bam.sort
	move bwa_mapped.bam.sort > "basename(*1.fastq _1.fastq).sorted.bam"
	"""

}

process consensus {
	
	input:
  	path map from bams
	path chlamy from consensuschlamy
	
	output:
	file "consensus.fa" into consensus
	
	
	"""
	bcftools mpileup -Ou -f $chlamy $map | bcftools call -mv -Oz -o calls.vcf.gz
	bcftools index calls.vcf.gz
	cat $chlamy | bcftools consensus calls.vcf.gz > consensus.fa
	"""

}







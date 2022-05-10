#!/usr/bin/env nextflow

myFile = 's3://pipe.scratch.3/resources/accessions.txt'
allLines = myFile.readLines()

seqdata= Channel.fromPath(allLines)

chlamyref='s3://pipe.scratch.3/resources/Chlamy23s.fasta'
refseq=Channel.fromPath(chlamyref)



process runfasta {
	
	input:
  	val accession from seqdata
	
	output:
	file "*_{1,2}.fastq" into dumpouts
	
	
	"""
	fastq-dump --split-3 $accession
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
	bwa mem $chlamy *1.fastq *2.fastq > bwa_mapped.sam
	samtools view -bS bwa_mapped.sam > bwa_mapped.bam
	samtools sort bwa_mapped.bam -o bwa_mapped.bam.sort
	mv bwa_mapped.bam.sort > "basename(*1.fastq _1.fastq).sorted.bam"
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







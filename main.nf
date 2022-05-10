#!/usr/bin/env nextflow

sraLines1=file('s3://pipe.scratch.3/resources/accessions.txt')
    .readLines()
    .each { println it }


chlamyref='s3://pipe.scratch.3/resources/Chlamy23s.fasta'


process runfasta {
	
	input:
  	val accession from sraLines1
	
	output:
	tuple val(accession), file("${accession}_1.fastq"), file("${accession}_2.fastq") into dumpout
	
	
	"""
	fastq-dump --split-3 $accession
	"""

}



process bwamap {
	
	input:
  	tuple val(accession), file(R1), file(R2) from dumpout
	path chlamy from chlamyref
	
	
	output:
	file "${accession}.sorted.bam" into bams
	
	
	"""
	bwa index $chlamy
	bwa mem $chlamy $R1 $R2 > bwa_mapped.sam
	samtools view -bS bwa_mapped.sam > bwa_mapped.bam
	samtools sort bwa_mapped.bam -o bwa_mapped.bam.sort
	mv bwa_mapped.bam.sort "${accession}.sorted.bam"
	"""

}

process consensus {
	
	input:
  	path map from bams
	path chlamy from chlamyref
	
	output:
	file "basename($map .sorted.bam)_cns.fastq" into consensus
	
	
	"""
	bcftools mpileup -Ou -f $chlamy $map | bcftools call -mv -Oz -o calls.vcf.gz
	bcftools mpileup -Ou -f $chlamy $map | bcftools call -c | vcfutils.pl vcf2fq > "basename($map .sorted.bam)_cns.fastq"
	"""

}







#!/usr/bin/env nextflow

sraLines=file('s3://pipe.scratch.3/resources/accessions.txt')
    .readLines()
    .each { println it }


chlamyref='s3://pipe.scratch.3/resources/Chlamy23s.fasta'
refseq=Channel.fromPath(chlamyref)


refseq.into{bwachlamy; consensuschlamy}

process runfasta {
	
	input:
  	val accession from sraLines
	
	output:
	file "*_1.fastq" into dumpout1
	file "*_2.fastq" into dumpout2
	
	
	"""
	fastq-dump --split-3 $accession
	"""

}



process bwamap {
	
	input:
  	path sradump1 from dumpout1
	path sradump2 from dumpout2
	path chlamy from bwachlamy
	
	output:
	file "basename($sradump1 _1.fastq).sorted.bam" into bams
	
	
	"""
	bwa index $chlamy
	bwa mem $chlamy $sradump1 $sradump2 > bwa_mapped.sam
	samtools view -bS bwa_mapped.sam > bwa_mapped.bam
	samtools sort bwa_mapped.bam -o bwa_mapped.bam.sort
	mv bwa_mapped.bam.sort > "basename($sradump1 _1.fastq).sorted.bam"
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







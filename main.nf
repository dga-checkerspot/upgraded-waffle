#!/usr/bin/env nextflow

sraLines1=file('s3://pipe.scratch.3/resources/accessions.txt')
    .readLines()
    .each { println it }

sraLines2=file('s3://pipe.scratch.3/resources/accessions.txt')
    .readLines()
    .each { println it }


chlamyref='s3://pipe.scratch.3/resources/Chlamy23s.fasta'
refseq=Channel.fromPath(chlamyref)


refseq.into{bwachlamy; consensuschlamy}

process runfasta {
	
	input:
  	val accession from sraLines1
	
	output:
	file "*.fastq" into dumpout
	
	
	"""
	fastq-dump --split-3 $accession
	"""

}



process bwamap {
	
	input:
	val accession from sraLines2
  	path sradump1 from dumpout.collect()
	path chlamy from bwachlamy
	
	output:
	file "$accession.sorted.bam" into bams
	
	
	"""
	bwa index $chlamy
	bwa mem $chlamy "$accession_1.fastq" "$accession_2.fastq" > bwa_mapped.sam
	samtools view -bS bwa_mapped.sam > bwa_mapped.bam
	samtools sort bwa_mapped.bam -o bwa_mapped.bam.sort
	mv bwa_mapped.bam.sort > "$accesion.sorted.bam"
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







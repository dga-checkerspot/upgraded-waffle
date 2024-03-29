#!/usr/bin/env nextflow

sraLines1=file('s3://pipe.scratch.3/resources/BigNCBISearchSeqsPedinophyceae.txt')
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


dumpout.into{dumpout1; dumpoutAssemble1; dumpoutAssemble2}




process bwamap {


	errorStrategy 'retry'
	memory '16G'
	
	input:
  	tuple val(accession), file(R1), file(R2) from dumpout1
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

	errorStrategy 'retry'
	memory '16G'
	
	input:
  	path map from bams
	path chlamy from chlamyref
	
	output:
	file "${map.baseName}.cns.fastq" into consensus
	
	
	"""
	bcftools mpileup -Ou -f $chlamy $map | bcftools call -c | vcfutils.pl vcf2fq > "${map.baseName}.cns.fastq"
	"""

}


process fasta {

	errorStrategy 'retry'

	input:
  	path fastqfile from consensus
	
	output:
	file "${fastqfile.baseName}_cns.fasta" into fasta


	"""
	seqtk seq -aQ64 -q15 -n N $fastqfile > "${fastqfile.baseName}_cns.fasta"
	"""
}


params.results = "s3://pipe.scratch.3/resources/ConsensusOutput/"

myDir = file(params.results)


fasta.subscribe { it.copyTo(myDir) }

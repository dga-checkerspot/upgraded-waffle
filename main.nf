#!/usr/bin/env nextflow

sraLines1=file('s3://pipe.scratch.3/resources/BigNCBISearchSeqs.txt')
    .readLines()
    .each { println it }


chlamyref='s3://pipe.scratch.3/resources/Chlamy23s.fasta'


process runfastadump {
	
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



process pileup {

	memory '64G'
	errorStrategy 'retry'
	
	input:
  	path map from bams
	path chlamy from chlamyref
	
	output:
	file "${map.baseName}.pileup.gz" into pileup
	
	"""
	bcftools mpileup -Ob -f $chlamy $map > ${map.baseName}.pileup.gz
	"""

}


process call {
	
	input:
  	path bcf from pileup
	
	output:
	file "${bcf.baseName}.vcf" into vcf
	
	
	"""
	bcftools call -mv -Ob  -o ${bcf.baseName}.vcf $bcf 
	"""

}



process fastq { 

	input:
  	path vcffile from vcf
	
	output:
	file "${vcffile.baseName}_cns.fastq" into consensus


	"""
	vcfutils.pl vcf2fq $vcffile > "${vcffile.baseName}_cns.fastq"
	"""
	
}


process fasta {

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


process read_data {
  
  input:
  path fastafile from fasta
  
  output:
  file "${fastafile.baseName}_final.fasta" into read

  """
  mv $fastafile "${fastafile.baseName}_final.fasta"
  """

}

read.subscribe { it.copyTo(myDir) }





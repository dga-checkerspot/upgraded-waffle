#!/usr/bin/env nextflow

seqlist='s3://pipe.scratch.3/resources/accessions.txt'
seqdata= Channel.fromPath(seqlist)

chlamyref='s3://pipe.scratch.3/resources/chlamy23s.fasta'
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




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


dumpout.into{dumpout1; dumpoutAssemble1; dumpoutAssemble2}




process bwamap {

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
	
	input:
  	path map from bams
	path chlamy from chlamyref
	
	output:
	file "${map.baseName}_cns.fastq" into consensus
	
	
	"""
	bcftools mpileup -Ou -f $chlamy $map | bcftools call -c | vcfutils.pl vcf2fq > "${map.baseName}_cns.fastq"
	"""

}

pairInt='s3://transcriptomepipeline/PairInterleaves.sh'


process cutadapt11 {
	memory '16G'
	
	input:
	tuple val(accession), file(R1), file(R2) from dumpoutAssemble1
	
	output:
	file("${accession}_1.fastq") into reads11
	file("${accession}_2.fastq") into reads12
	
	"""
	cutadapt --rename='{id}/1' $R1 -j 0 -o "${accession}_1.fastq"
	cutadapt --rename='{id}/2' $R2 -j 0 -o "${accession}_1.fastq"
	"""
}


process bbnorm {

	memory '196G'
	
        input:
        path seq1 from reads11
        path seq2 from reads12
        
        output:
        file("${seq1.baseName}.mid.fastq") into ReadTrimNorm1

	"""
	bbnorm.sh in=$seq1 in2=$seq2 outlow=low.fq outmid=mid.fq outhigh=high.fq passes=1 lowbindepth=6 highbindepth=150 -Xmx192g
	mv mid.fq "${seq1.baseName}.mid.fastq"
	"""
}



process pairInt {

	memory '4G'

	input:
	path 'pairInt' from pairInt
	path 'Intpair' from ReadTrimNorm1

	output:
	file '${Intpair.baseName}.R1reads.fastq' into R1Tofastq
	file '${Intpair.baseName}.R2reads.fastq' into R2Tofastq

	"""
	chmod 744 $pairInt
	./$pairInt < $Intpair "${Intpair.baseName}.R1reads.fastq" "${Intpair.baseName}.R2reads.fastq"
	"""

}


process fastqpair2 {

	memory '32G'

	input:
	path R1p from R1Tofastq
	path R2p from R2Tofastq

	output:
	file '$R1p.paired.fq' into pairR1T
	file '$R2p.paired.fq' into pairR2T
	//For now not even bothering with unpaired

	"""
	fastq_pair -t 100000000 $R1p $R2p
	"""
}

pairR1T.into{P1NormSpades; P1NormTrinity}
pairR2T.into{P2NormSpades; P2NormTrinity}


process SpadeAssemble {
	
  memory '24G'

  input:
    path R1Norm from P1NormSpades
    path R2Norm from P2NormSpades

  output:
    file '${R1Norm.baseName}.spades.tar.gz' into Spades
    
    """
    rnaspades.py  --pe1-1 $R1Norm --pe1-2 $R2Norm  -o spades_output
    tar -zcvf ${R1Norm.baseName}.spades.tar.gz spades_output 
    
    """
    
    
}


process TrinityAssemble {
	
  memory '196G'
	
  input:
	path R1pair from P1NormTrinity
	path R2pair from P2NormTrinity
	
  output:
	file '${R1pair.baseName}.trinity.tar.gz' into Trinity
	
  """
	Trinity --seqType fq --left $R1pair --right $R2pair --max_memory 190G --output trinity_output
	tar -zcvf ${R1pair.baseName}.trinity.tar.gz trinity_output 
	"""

}












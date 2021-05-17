
configfile: "config.yaml"

import pandas as pd

sample_table = pd.read_table(config["runtable"], sep="\t")


samples = list(sample_table["sampleName"].unique())

RENAME = list(sample_table["sraName"].unique())

Rplot = sample_table['bedgraph'] == "x"
Rplot = sample_table[Rplot]
Rplot = list(Rplot["sampleName"].unique())


# everything that need network connection for
localrules:

rule all:
	input:
		expand("data/RSEM/{SRS}.genes.results", SRS=samples),
		expand("data/bedgraph/{SRS}.BedGraph", SRS=samples),
		expand("data/bedgraph/{SRS}.txt", SRS=Rplot)


rule index:
	input:
		fasta="genomeFiles/ASM294v2.27.fasta",
		gtf="genomeFiles/sortedStrandedAllLTRs.gtf"
	output:
		"data/index/Genome",
		"data/index/genome.transcripts.fa",
		"data/index/chrNameLength.txt"
	shell:
		"""
		rsem-prepare-reference --gtf {input.gtf} --star {input.fasta} data/index/genome
		"""


rule rename:
	input:
		expand("data/sraFiles/{rename}.fastq.gz", rename=RENAME)
	output:
		"data/sraFiles/{SRS}.fastq.gz"
	shell:
		"""
		python rename.py
		"""


rule STAR:
	input:
		file1="data/sraFiles/{SRS}.fastq.gz",
		index="data/index/Genome"
	output:
		temp("data/SAMBAM/{SRS}.Aligned.toTranscriptome.out.bam"),
		temp("data/SAMBAM/{SRS}.Aligned.sortedByCoord.out.bam")

	threads: 16

	shell:
		"""
		STAR --genomeDir data/index \
		--runThreadN {threads} \
		--readFilesIn <(gunzip -c {input.file1}) \
		--sjdbGTFfile genomeFiles/sortedStrandedAllLTRs.gtf \
		--outFileNamePrefix data/SAMBAM/{wildcards.SRS}. \
		--quantMode GeneCounts TranscriptomeSAM \
		--outSAMtype BAM SortedByCoordinate
		"""


rule RSEM:
	input:
		alignedBam="data/SAMBAM/{SRS}.Aligned.toTranscriptome.out.bam",
		index="data/index/Genome"
	output:
		"data/RSEM/{SRS}.genes.results",
		temp("data/RSEM/{SRS}.transcript.bam")

	threads: 16

	shell:
		"""
		rsem-calculate-expression --strandedness reverse --bam {input.alignedBam} -p {threads} data/index/genome data/RSEM/{wildcards.SRS}
		"""


rule bedtools:
	input:
		"data/SAMBAM/{SRS}.Aligned.sortedByCoord.out.bam"
	output:
		"data/bedgraph/{SRS}.BedGraph"
	shell:
		"""
		TmpScale=$(bc <<< "scale=6;1000000/$(samtools view -f 0 -c {input})")
		bedtools genomecov -ibam {input} -bg -split -scale $TmpScale  > data/bedgraph/{wildcards.SRS}.BedGraph
		"""


rule bedtoolsForR:
	input:
		"data/SAMBAM/{SRS}.Aligned.sortedByCoord.out.bam"
	output:
		"data/bedgraph/{SRS}.txt"
	shell:
		"""
		TmpScale=$(bc <<< "scale=6;1000000/$(samtools view -f 0 -c {input})")
		bedtools genomecov -ibam {input} -d -split -scale $TmpScale  > data/bedgraph/{wildcards.SRS}.txt
		"""


onsuccess:
        shell("rm *.out")

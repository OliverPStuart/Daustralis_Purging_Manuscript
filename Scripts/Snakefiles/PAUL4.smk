### Snakefile for PAUL4 analysis. The pipeline takes reads, maps and
### deduplicates them, calls variants in parallel with FreeBayes, and filters
### them. Then, it measures genome wide H in 1mbp bins with 500kbp sliding windows,
### estimates ROH coordinates with bcftools.

### May add a rule to estimate recombination map using ismc.

HOME_DIR = "/mnt/data/dayhoff/home/scratch/groups/mikheyev/LHISI"
READ_DIR = HOME_DIR + "/Reads/WGS"
ALN_DIR = HOME_DIR + "/Alignments/WGS"
REF_DIR = HOME_DIR + "/References"
OUT_DIR = HOME_DIR + "/Analysis/PAUL4"
SOFT_DIR = "/mnt/data/dayhoff/home/u6905905"
BCFTOOLS = SOFT_DIR + "/bcftools/bcftools"
REF = REF_DIR + "/LHISI_Scaffold_Assembly.fasta"


### Modify this line to exclude or include scaffolds as necessary
### Only include autosomes for this
CHROMS = ["CM056993.1","CM056994.1","CM056995.1","CM056996.1","CM056997.1","CM056998.1","CM056999.1","CM057000.1","CM057001.1","CM057002.1","CM057003.1","CM057004.1","CM057005.1","CM057006.1","CM057007.1","CM057008.1"]
#CHROMS = ["CM057008.1"] # tester line

### Define chunks for calling later
NCHUNKS = 30
CHUNKS = list(range(1,NCHUNKS+1))

### Define intervals and reps for PSMC testing
#INTERVALS = [15, 20, 25]
NREPS = 100
REPS = list(range(1,NREPS+1))

### Final rule
rule all:
	input:
		OUT_DIR + "/PAUL4_WGS.regions.bed.gz",
		OUT_DIR + "/GENOME_WIDE_H.txt",
		OUT_DIR + "/GENOME_WIDE_ROH_04_BCFTOOLS.txt",
		OUT_DIR + "/GENOME_WIDE_ROH_05_BCFTOOLS.txt",
		OUT_DIR + "/GENOME_WIDE_ROH_06_BCFTOOLS.txt"

### Map the files to the reference and follow samtools pipeline to remove duplicates
rule map:
	input:
		ref = REF,
		r1 = READ_DIR + "/PAUL4_WGS_R1.fastq.gz",
		r2 = READ_DIR + "/PAUL4_WGS_R2.fastq.gz"
	output:
		ALN_DIR + "/PAUL4_WGS.bam"
	params:
		dir = ALN_DIR
	threads: 46
	shell:
		"""
		bwa mem -t {threads} -R "@RG\\tID:PAUL4\\tSM:PAUL4\\tLB:NEXTFLEX\\tPL:ILLUMINA" \
		{input.ref} {input.r1} {input.r2} |
		samtools fixmate -m -r -u - - | \
		samtools sort -u -@ {threads} -T {params.dir}/PAUL4_WGS_1 - | \
		samtools markdup -@ {threads} -T {params.dir}/PAUL4_WGS_2 --reference {input.ref} -r -f {params.dir}/PAUL4_WGS_dupstats.txt - - | \
		samtools view -bh > {output}
		"""

### Index the bam files
rule index:
	input:
		ALN_DIR + "/PAUL4_WGS.bam"
	output:
		ALN_DIR + "/PAUL4_WGS.bam.bai"
	threads: 46
	shell:
		"""
		samtools index -@ {threads} {input}
		"""

### Rule to get depth statistics
rule get_depth:
	input:
		bam = ALN_DIR + "/PAUL4_WGS.bam",
		index = ALN_DIR + "/PAUL4_WGS.bam.bai"
	output:
		OUT_DIR + "/PAUL4_WGS.regions.bed.gz"
	threads: 46
	params:
		out=OUT_DIR
	shell:
		"""
		cd {params.out}
		mosdepth -t {threads} -n -Q 10 -b 100000 PAUL4_WGS {input.bam}
		rm *mosdepth*
		"""

### Rule to make bed files to split up calling
rule generate_freebayes_regions:
	input:
		ref_idx = REF,
		index = REF_DIR + "/LHISI_Scaffold_Assembly.fasta.fai",
	output:
		regions = expand(REF_DIR + "/bed_files/genome.{chrom}.region.{i}.bed", chrom=CHROMS, i=CHUNKS)
	params:
		chroms = CHROMS,
		chunks = NCHUNKS,
		soft = SOFT_DIR,
		ref_dir = REF_DIR
	shell:
		"""
		{params.soft}/freebayes/scripts/fasta_generate_regions.py \
		--chunks \
		--chromosomes {params.chroms} \
		--bed {params.ref_dir}/bed_files/genome \
		{input.index} \
		{params.chunks}
		"""

### Rule to do calls
rule variant_calling_freebayes:
	input:
		bam = ALN_DIR + "/PAUL4_WGS.bam",
		index = ALN_DIR + "/PAUL4_WGS.bam.bai",
		ref = REF,
		regions = REF_DIR + "/bed_files/genome.{chrom}.region.{i}.bed"
	output:
		OUT_DIR + "/{chrom}/variants.{i}.vcf"
	wildcard_constraints:
		i="[0-9]+"
	params:
		soft = SOFT_DIR
	threads: 1
	wildcard_constraints:
		i="[0-9]+"
	shell:
		"""
		{params.soft}/freebayes_1.3.6 \
		-f {input.ref} \
		-t {input.regions} \
		--use-best-n-alleles 2 \
		 --report-monomorphic \
		-0 --min-coverage 10 -g 30 \
		{input.bam} > {output}
		"""

### Rule to filter vcfs, treating variant and invariant sites separately
rule filter_vcfs:
	input: OUT_DIR + "/{chrom}/variants.{i}.vcf"
	output:
		out_het = OUT_DIR + "/{chrom}/variants.{i}.het.vcf",
		out_hom = OUT_DIR + "/{chrom}/variants.{i}.hom.vcf"
	threads:2
	params:
		dir = OUT_DIR,
		ref_dir = REF_DIR
	wildcard_constraints:
		i="[0-9]+"
	shell:
		"""
		cd {params.dir}

		# Intersect specific bedfile with repeats bed file to reduce file-reading for vcftools

		bedtools intersect -a {params.ref_dir}/bed_files/genome.{wildcards.chrom}.region.{wildcards.i}.bed -b {params.ref_dir}/Annotation/repeats_major.bed > filter_{wildcards.chrom}_{wildcards.i}.bed

		# Filter by depth first
		# This will significantly reduce the burden of later sections
		# Permissive filter

		vcftools --vcf {input} \
		--min-meanDP 15 --max-meanDP 25 \
		--remove-indels \
		--exclude-bed filter_{wildcards.chrom}_{wildcards.i}.bed \
		--recode --recode-INFO-all --out depthfilt_{wildcards.chrom}_{wildcards.i}

		rm filter_{wildcards.chrom}_{wildcards.i}.bed

		# Split into header and tailer

		grep "#" depthfilt_{wildcards.chrom}_{wildcards.i}.recode.vcf > {wildcards.chrom}_{wildcards.i}_header
		grep -v "#" depthfilt_{wildcards.chrom}_{wildcards.i}.recode.vcf > {wildcards.chrom}_{wildcards.i}_tailer
		rm depthfilt_{wildcards.chrom}_{wildcards.i}.recode.vcf

		# Split tailer into hom and het

		grep "0/0\|1/1\|0|0\|1|1" {wildcards.chrom}_{wildcards.i}_tailer > {wildcards.chrom}_{wildcards.i}_hom
		grep "0/1\|1/0\|0|1\|1|0" {wildcards.chrom}_{wildcards.i}_tailer > {wildcards.chrom}_{wildcards.i}_het

		# Make full hom vcf

		cat {wildcards.chrom}_{wildcards.i}_header {wildcards.chrom}_{wildcards.i}_hom > {output.out_hom}

		# Get full het vcf for filtering

		cat {wildcards.chrom}_{wildcards.i}_header {wildcards.chrom}_{wildcards.i}_het > {wildcards.chrom}_{wildcards.i}_het.vcf

		# Remove tmp files

		rm {wildcards.chrom}_{wildcards.i}_tailer {wildcards.chrom}_{wildcards.i}_het {wildcards.chrom}_{wildcards.i}_hom {wildcards.chrom}_{wildcards.i}_header

		# Filter variant VCF, easy part

		vcffilter -f 'AB > 0.2 & AB < 0.8 & QUAL / DP > 0.25' {wildcards.chrom}_{wildcards.i}_het.vcf | grep -v "#" > {output.out_het}

		# Remove other tmp files

		rm {wildcards.chrom}_{wildcards.i}_het.vcf
		"""

### Rule to concat calls into chromosome specific directories
rule concat_vcfs:
	input:
		calls_het = expand(OUT_DIR + "/{{chrom}}/variants.{i}.het.vcf", i=CHUNKS),
		calls_hom = expand(OUT_DIR + "/{{chrom}}/variants.{i}.hom.vcf", i=CHUNKS),
		index = REF_DIR + "/LHISI_Scaffold_Assembly.fasta.fai"
	output:
		out_het = OUT_DIR + "/PAUL4_Het_{chrom}.vcf",
		out_hom = OUT_DIR + "/PAUL4_Hom_{chrom}.vcf"
	threads:4
	params:
		dir = OUT_DIR
	shell:
		"""
		cd {params.dir}

		# We have to separate processing the header into two lines to avoid 141 error
		# This occurs when a pipe is still streaming after the following command exits

		cat {input.calls_hom} | awk '/#/' > tmp2_{wildcards.chrom}
		awk '1;/CHROM/{{exit}}' tmp2_{wildcards.chrom} >> head_{wildcards.chrom}

		cat {input.calls_het} | awk '!/#/' | sort -k2,2n > het_{wildcards.chrom}
		cat head_{wildcards.chrom} het_{wildcards.chrom} | vcfuniq > {output.out_het}

		cat {input.calls_hom} | awk '!/#/' | sort -k2,2n > hom_{wildcards.chrom}
		cat head_{wildcards.chrom} hom_{wildcards.chrom} | vcfuniq > {output.out_hom}

		rm head_{wildcards.chrom} het_{wildcards.chrom} hom_{wildcards.chrom} tmp2_{wildcards.chrom}
		"""

### Rule to get table of per-scaffold heterozygosity in sliding windows
rule H:
	input:
		het = OUT_DIR + "/PAUL4_Het_{chrom}.vcf",
		hom = OUT_DIR + "/PAUL4_Hom_{chrom}.vcf"
	output: OUT_DIR + "/H_count_{chrom}.txt"
	params:
		ref_dir = REF_DIR,
		dir = OUT_DIR
	shell:
		"""
		cd {params.dir}

		# Get length of scaffold
		cut -f1,2 {params.ref_dir}/LHISI_Scaffold_Assembly.fasta.fai | grep {wildcards.chrom} > {wildcards.chrom}_length

		# Make bedfile
		bedtools makewindows -g {wildcards.chrom}_length -w 1000000 -s 500000 > {wildcards.chrom}_windows.bed

		# Get file of site positions for variant, invariant
		awk '!/#/' {input.het} | cut -f1,2 | awk '{{print $1"\t"$2-1"\t"$2}}'	> {wildcards.chrom}_het.bed
		awk '!/#/' {input.hom} | cut -f1,2 | awk '{{print $1"\t"$2-1"\t"$2}}'	> {wildcards.chrom}_hom.bed

		# Now we intersect each to get the count of variant and invariant sites per window
		bedtools intersect \
		-a {wildcards.chrom}_windows.bed \
		-b {wildcards.chrom}_het.bed \
		-c | \
		cut -f4 > {wildcards.chrom}_het_count

		bedtools intersect \
		-a {wildcards.chrom}_windows.bed \
		-b {wildcards.chrom}_hom.bed \
		-c > {wildcards.chrom}_hom_count

		# And paste them together into one table
		paste {wildcards.chrom}_hom_count {wildcards.chrom}_het_count > {output}

		# And clean up
		rm {wildcards.chrom}_*bed {wildcards.chrom}_*count {wildcards.chrom}_length
		"""

### Rule to combine H estimates per chromosome into single table for analysis
rule combine_H:
	input: expand(OUT_DIR + "/H_count_{chrom}.txt",chrom=CHROMS)
	output: OUT_DIR + "/GENOME_WIDE_H.txt"
	params:
		dir = OUT_DIR
	shell:
		"""
		cd {params.dir}

		cat {input} > {output}

		rm {input}
		"""

### Rule to run bcftools ROH analysis
rule bcftools_roh:
	input:
		het = OUT_DIR + "/PAUL4_Het_{chrom}.vcf",
		hom = OUT_DIR + "/PAUL4_Hom_{chrom}.vcf"
	output:
		o04 = OUT_DIR + "/{chrom}_ROH_04.txt",
		o05 = OUT_DIR + "/{chrom}_ROH_05.txt",
		o06 = OUT_DIR + "/{chrom}_ROH_06.txt"
	params:
		dir = OUT_DIR,
		bcftools = BCFTOOLS
	shell:
		"""
		cd {params.dir}

		# Get head and het/hom tails
		grep "#" {input.het} > {wildcards.chrom}_head
		grep -v "#" {input.het} > {wildcards.chrom}_tmp1
		grep -v "#" {input.hom} > {wildcards.chrom}_tmp2

		# Combine het and hom tails, then cat into single vcf
		cat {wildcards.chrom}_tmp1 {wildcards.chrom}_tmp2 | sort -k2,2n > {wildcards.chrom}_tail
		cat {wildcards.chrom}_head {wildcards.chrom}_tail > PAUL4_BOTH_{wildcards.chrom}.vcf

		# Clean up
		rm {wildcards.chrom}_tmp1 {wildcards.chrom}_tmp2 {wildcards.chrom}_head {wildcards.chrom}_tail

		# Perform bcftools ROH estimation with three AF settings
		{params.bcftools} roh --AF-dflt 0.4 -G 30 -o {wildcards.chrom}_tmp_o04  PAUL4_BOTH_{wildcards.chrom}.vcf
		{params.bcftools} roh --AF-dflt 0.5 -G 30 -o {wildcards.chrom}_tmp_o05 PAUL4_BOTH_{wildcards.chrom}.vcf
		{params.bcftools} roh --AF-dflt 0.6 -G 30 -o {wildcards.chrom}_tmp_o06 PAUL4_BOTH_{wildcards.chrom}.vcf

		# Clean up
		rm PAUL4_BOTH_{wildcards.chrom}.vcf

		# Clean up table
		grep RG {wildcards.chrom}_tmp_o04 | sed 's/# //g' | cut -f3- | sed 's/\[[0-9]\]//g' | sed 's/ //g' > {output.o04}
		grep RG {wildcards.chrom}_tmp_o05 | sed 's/# //g' | cut -f3- | sed 's/\[[0-9]\]//g' | sed 's/ //g' >{output.o05}
		grep RG {wildcards.chrom}_tmp_o06 | sed 's/# //g' | cut -f3- | sed 's/\[[0-9]\]//g' | sed 's/ //g' > {output.o06}

		# Clean up
		rm {wildcards.chrom}_tmp_o04 {wildcards.chrom}_tmp_o05 {wildcards.chrom}_tmp_o06
		"""

### Rule to combine bcftools ROH estimates into single table
rule bctools_roh_combine:
	input:
		i04 = expand(OUT_DIR + "/{chrom}_ROH_04.txt",chrom=CHROMS),
		i05 = expand(OUT_DIR + "/{chrom}_ROH_05.txt",chrom=CHROMS),
		i06 = expand(OUT_DIR + "/{chrom}_ROH_06.txt",chrom=CHROMS)
	output:
		o04 = OUT_DIR + "/GENOME_WIDE_ROH_04_BCFTOOLS.txt",
		o05 = OUT_DIR + "/GENOME_WIDE_ROH_05_BCFTOOLS.txt",
		o06 = OUT_DIR + "/GENOME_WIDE_ROH_06_BCFTOOLS.txt"
	params:
		dir = OUT_DIR
	shell:
		"""
		cd {params.dir}

		cat {input.i04} | sort | uniq > {output.o04}
		cat {input.i05} | sort | uniq > {output.o05}
		cat {input.i06} | sort | uniq > {output.o06}

		rm {input.i04} {input.i05} {input.i06}
		"""

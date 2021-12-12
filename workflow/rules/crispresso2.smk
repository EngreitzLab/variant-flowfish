## Run CRISPResso for each FASTQ file

# To do: consider adding  --max_paired_end_reads_overlap for FLASH overlap (calculate from read lengths, and length of amplicon)


if single_end:
	rule run_crispresso:
		input:
			read1=lambda wildcards: samplesheet.at[wildcards.SampleID,'fastqR1']
		output:
			directory('results/crispresso/CRISPResso_on_{SampleID}/')
		params:
			amplicon_id=lambda wildcards: samplesheet.at[wildcards.SampleID,'AmpliconID'],
			amplicon_seq=lambda wildcards: samplesheet.at[wildcards.SampleID,'AmpliconSeq'],
			q=config['crispresso_min_average_read_quality'],
			s=config['crispresso_min_single_bp_quality'],
			qws=lambda wildcards: samplesheet.at[wildcards.SampleID,'QuantificationWindowStart'],
			qwe=lambda wildcards: samplesheet.at[wildcards.SampleID,'QuantificationWindowEnd'],
		shell:
			"""
			bash -c '
				. $HOME/.bashrc 
				conda activate crispresso2_v2.2.6
				CRISPResso \
					-r1 {input.read1} \
					-o results/crispresso/ \
					--amplicon_seq {params.amplicon_seq} \
					--amplicon_name {params.amplicon_id} \
					--name {wildcards.SampleID} \
					--quantification_window_coordinates {params.qws}-{params.qwe} \
					--exclude_bp_from_left 0 --exclude_bp_from_right 0 \
					--plot_window_size 0 --min_frequency_alleles_around_cut_to_plot 0.05 \
					-q {params.q} -s {params.s} || true'
			"""

else:
	rule run_crispresso:
		input:
			read1=lambda wildcards: samplesheet.at[wildcards.SampleID,'fastqR1'],
			read2=lambda wildcards: samplesheet.at[wildcards.SampleID,'fastqR2']
		output:
			directory('results/crispresso/CRISPResso_on_{SampleID}/')
			#'crispresso/CRISPResso_on_{SampleID}/{AmpliconID}.Alleles_frequency_table.txt'
		params:
			amplicon_id=lambda wildcards: samplesheet.at[wildcards.SampleID,'AmpliconID'],
			amplicon_seq=lambda wildcards: samplesheet.at[wildcards.SampleID,'AmpliconSeq'],
			q=config['crispresso_min_average_read_quality'],
			s=config['crispresso_min_single_bp_quality'],
		#conda:
		#    "envs/CRISPResso.yml"  
		## 4/14/21 JE - Specifying the conda environment here is not working, and I am not sure why. Snakemake builds the conda environment, but then the conda environment doesn't work properly (CRISPResso not on the path)
		#   (This was on Sherlock, running snakemake from EngreitzLab conda envrionment).
		#  So, instead used the syntax below to activate the already installed conda env
		shell:
			"""
			bash -c '
				. $HOME/.bashrc 
				conda activate crispresso2_v2.2.6
				CRISPResso \
					-r1 {input.read1} \
					-r2 {input.read2} \
					-o results/crispresso/ \
					--amplicon_seq {params.amplicon_seq} \
					--amplicon_name {params.amplicon_id} \
					--name {wildcards.SampleID} \
					-q {params.q} -s {params.s} || true'
			"""

	## TO do — add in quantification window for this command, as in single-end section
	#				--quantification_window_coordinates {params.qws}-{params.qwe} \
	#				--exclude_bp_from_left 0 --exclude_bp_from_right 0 \
	#				--plot_window_size 0 --min_frequency_alleles_around_cut_to_plot 0.05 \

## Run CRISPRessoAggregate to collect mapping statistics across all runs
rule run_crispresso_aggregate:
	input:
		['results/crispresso/CRISPResso_on_{s}/'.format(s=s) for s in samplesheet['SampleID']]
	output:
		directory('results/crispresso/CRISPRessoAggregate_on_Aggregate/')
	shell:
		"""
		bash -c '
			. $HOME/.bashrc 
			conda activate crispresso2_v2.2.6
			cd results/crispresso/
			CRISPRessoAggregate --name Aggregate --prefix ./CRISPResso_on_ || true'
		"""



## Align all FASTQs to all amplicons to check that they're the correct samples
import io

rule create_fasta:
  input:
  output:
    fasta='amplicons/Amplicons.fa'
  run:
    f = open(output.fasta, 'w')
    for i, row in samplesheet[['AmpliconID','AmpliconSeq']].drop_duplicates().iterrows():
      f.write('>' + row['AmpliconID'] + "\n" + row['AmpliconSeq'] + "\n")
    f.close()


rule create_bowtie2_index:
  input:
    fasta='amplicons/Amplicons.fa'
  output:
    index='amplicons/Amplicons.fa.1.bt2'
  shell:
    """
    bash -c '
      . $HOME/.bashrc 
      conda activate crispresso2_210104
      bowtie2-build {input.fasta} {input.fasta}
    '
    """

rule run_bowtie2:
  input:
    read1=lambda wildcards: samplesheet.at[wildcards.SampleID,'fastqR1'],
    read2=lambda wildcards: samplesheet.at[wildcards.SampleID,'fastqR2'],
    fasta='amplicons/Amplicons.fa',
    index='amplicons/Amplicons.fa.1.bt2'
  output:
    bam='results/aligned/{SampleID}/{SampleID}.bam',
    bai='results/aligned/{SampleID}/{SampleID}.bam.bai',
    unaligned_R1='results/aligned/{SampleID}/{SampleID}_unaligned.fastq.1.gz',
    unaligned_R2='results/aligned/{SampleID}/{SampleID}_unaligned.fastq.2.gz'
  params:
    amplicon_seq=lambda wildcards: samplesheet.at[wildcards.SampleID,'AmpliconSeq'],
    guide=lambda wildcards: samplesheet.at[wildcards.SampleID,'GuideSpacer'],
    q=config['crispresso_min_average_read_quality'],
    s=config['crispresso_min_single_bp_quality'],
    unaligned='results/aligned/{SampleID}/{SampleID}_unaligned.fastq.gz',
    codedir=codedir
  #conda:
  #    "envs/CRISPResso.yml"  
  ## 4/14/21 JE - Specifying the conda environment here is not working, and I am not sure why. Snakemake builds the conda environment, but then the conda environment doesn't work properly (CRISPResso not on the path)
  #   (This was on Sherlock, running snakemake from EngreitzLab conda envrionment).
  #  So, instead used the syntax below to activate the already installed conda env
  shell:
    """
    bash -c '
      . $HOME/.bashrc 
      conda activate EngreitzLab
      mkdir -p tmp
      bowtie2 -x {input.fasta} \
          -1 {input.read1} \
          -2 {input.read2} \
          --un-conc-gz {params.unaligned} \
          --very-sensitive-local \
          | samtools sort -T tmp/sort.{wildcards.SampleID} -O bam -o {output.bam} - && samtools index {output.bam}'
    """


## Count the number of reads aligned to each amplicon, using samtools idxstats
rule count_bowtie2_alignments:
  input:
    expand("results/aligned/{SampleID}/{SampleID}.bam", SampleID=samplesheet['SampleID'])
  output:
    table="results/summary/alignment.counts.tsv"
  run:
    resultString = "SampleID\t" + next(shell("samtools idxstats {} | cut -f 1 | tr '\n' '\t'".format(input[0]), iterable=True)) + "\n"
    for file in input:
      currOut = next(shell("samtools idxstats {} | cut -f 3 | tr '\n' '\t'".format(file), iterable=True))
      resultString = resultString + os.path.splitext(os.path.basename(file))[0] + "\t" + currOut + "\n"
    result = pd.read_csv(io.StringIO(resultString), sep='\t')
    result = result.merge(samplesheet.reset_index(drop=True))
    result.to_csv(output.table, sep='\t', header=True, index=False)

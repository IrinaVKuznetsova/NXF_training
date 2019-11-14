/*
 * 13 November 2019 |
 */


/*
 * Step 1. Define the default parameters
 */
params.genome     = "/home/user1/rnaseq_mito/reference/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.reads      = "/home/user1/rnaseq_mito/Inputdata/*_{1,2}.fq" 
params.results    = "results" 


/*
 *  Step 2. Parse the input parameters
 */
genome_file  = file(params.genome)
Channel.fromFilePairs(params.reads).into{reads_ch; reads_ch2}
/*
 *  RUN to check: nextflow run main_mt.nf
 */


/*
 * PROCESS 1 | Step 3. Create a FASTA genome index with samtools
 * use existing container from biocontainers :  biocontainers/samtools:v1.9-4-deb_cv1
 */

process 'prepare_indexfile_genome_samtools' { 

  container 'biocontainers/samtools:v1.9-4-deb_cv1'
  input:
      file genome from genome_file 

  output:
      tuple file(genome), file("${genome}.fai") into genome_index_ch 

  script:
  """
  samtools faidx ${genome} 
  """
}

/*
 * PROCESS 2 | Step 4. Check quality with FASTQC
 */

process 'check_Fa_quality_FASTQC' {
    publishDir 'results_gual'
    container 'biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1'

  input:
      tuple val(sampleId), path(reads) from reads_ch

  output:
      file "fastqc_results/*_fastqc.{zip,html}" into results_gual
      
      
  script:
  """
  mkdir fastqc_results
  fastqc $reads -o fastqc_results
  """
}




/*
 * PROCESS 3 | Step 5. Trimming
 */
process 'Trimming_cutadapt' {
    publishDir 'results_trim'
    container 'quay.io/biocontainers/cutadapt:2.6--py36h516909a_0'

  input:
      tuple val(sampleId), path(reads) from reads_ch2

  output:
      file("*_trimmed.fq") into results_trim
 

  script:
    """
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${reads[0]}_trimmed.fq -p ${reads[1]}_trimmed.fq $reads
    """
}



/*
 * PROCESS 4 | Step 6. Alignment with HISAT2
 
process 'Reads_Alignment' {
    publishDir 'results_align'
    container 'bioconda/label/cf201901'
    1.3.1--h0592bc0_6

  input:
      

  output:
      file() into 
 
  script:
    """
    
    """
}
*/
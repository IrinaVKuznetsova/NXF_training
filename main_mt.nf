/*
 * 13 November 2019 |
 */


// USEFUL HINTS: 
// https://hub.docker.com/r/makaho/hisat2-zstd/tags
// reads_ch.view() - to run script that is  above this command
// return
// cd work/e7/  - to see content of the channel



/*
 * Step 1. Define the default parameters
 */
params.genome     = "/home/user1/NXF_training/test-data/reference/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.reads      = "/home/user1/NXF_training/test-data/reads/*_{1,2}.fq" 
params.results    = "results" 


/*
 *  Step 2. Parse the input parameters
 */
genome_file  = file(params.genome)   // value channel because if 10 processes ask for this value it will find it 10 times | in is the same out
Channel.fromFilePairs(params.reads).into{reads_ch; reads_ch2} // channelFrom is queue channel  | beacuse inpt is queued output is queued
/*
 *  RUN to check: nextflow run main_mt.nf
 */





/*
 * PROCESS 2 | Step 3. Check quality with FASTQC
 */
process 'check_Fa_quality_FASTQC' {
    container 'biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1'
  input:
      tuple val(sampleId), path(reads) from reads_ch
  output:
      file "*_fastqc.{zip,html}" into results_qual_ch
  script:
  """
   fastqc $reads -o .
  """
}






/*
 * PROCESS 3 | Step 3. Trimming
 */
process 'Trimming_cutadapt' {
    container 'quay.io/biocontainers/cutadapt:2.6--py36h516909a_0'
  input:
      tuple val(sampleId), path(reads) from reads_ch2
  output:
      file("*_trimmed.fq") into results_trim_ch
   script:
    """
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${reads[0]}_trimmed.fq -p ${reads[1]}_trimmed.fq $reads
    """
}


/*
 * PROCESS 4 | Step 4. Get INDEX file with HISAT2-build
 * docker hub
 * hisat2-build [options]* <reference_in> <ht2_base>
 * tuple is used as we need Fa and Fai together in future
 */
process 'Get indexed file HISAT2' { 
  container 'makaho/hisat2-zstd:latest'
  input:
      file(genome) from genome_file 
  output:
      file("indexed.*.ht2") into(genome_index_ch)  //, genome_aligned_ch2)
  script:
  """
  hisat2-build $genome indexed 
  """
}


/*
 * PROCESS 5 | Step 5. Aligne with HISAT2
 * hisat2-build [options]* <reference_in> <ht2_base>
 */
process 'Align HISAT2' { 
  container 'makaho/hisat2-zstd:latest'
  input:
      file(reads) from results_trim_ch 
      file(indexed) from genome_index_ch
  output:
      tuple file("OUT_aligned_ggal.sam"), file("stats.log") into sam_channel
  script:
  """
  hisat2 -p 20 --no-spliced-alignment -x indexed -1 ${reads[0]} -2 ${reads[1]} -S OUT_aligned_ggal.sam &> stats.log
  """
}


/*
 * PROCESS 6 | Step 6. Sam to Bam
 * samtools view -S -b sample.sam > sample.bam
 */
process 'Convert SAM to BAM' { 
  container 'jweinstk/samtools:latest'
  input:
      file(aligned_sam) from sam_channel 
  output:
      file("OUT_aligned_ggal.bam") into bam_channel
  script:
  """
  samtools view -S -b aligned_sam > OUT_aligned_ggal.bam
  
  """
}




///samtools view -S -b OUT_aligned_ggal.sam > OUT_aligned_ggal.bam

// publishDir - specify directory where you want to copy files


/*
 * PROCESS 4 | Step 4. Get INDEX file with HISAT2-build
 * docker hub
 * hisat2-build [options]* <reference_in> <ht2_base>
 * tuple is used as we need Fa and Fai together in future
 
process 'Get indexed file HISAT2' { 
  container 'makaho/hisat2-zstd:latest'
  input:
      file(genome) from genome_file 
      file(reads) from results_trim_ch
  output:
      tuple file(genome), file("indexed.*.ht2") into(genome_index_ch)  //, genome_aligned_ch2)
      tuple file("OUT_aligned_ggal.sam") into sam_channel  
  script:
  """
  hisat2-build $genome indexed 
  hisat2 -p 20 --no-spliced-alignment -x indexed -1 ${reads[0]} -2 ${reads[1]} -S OUT_aligned_ggal.sam 
  samtools SAM to BAM
  """
}
*/

//reads_ch.view() 
//return

// hisat2 -p 20 --no-spliced-alignment -x $index -1 $ -2 $ -S OUT_aligned_ggal.sam &> OUT_STATS_ggal.log

/*
 * PROCESS 4 | Step 6. Alignment with HISAT2
 process 'Reads_Alignment' {
    container 'makaho/hisat2-zstd:latest'
 input:
    file trimmed_reads from results_trim_ch
      

  output:
      file() into 
 
  script:
    """
    hisat2 -p 20 --no-spliced-alignment -x $index -1 $ -2 $ -S aligned_ggal.sam &> aligned_ggal_stats.log
    """
}
*/
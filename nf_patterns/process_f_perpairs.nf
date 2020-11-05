/*
* This script shows how to fetch all FASTQ file pairs from a directory
* passed in the command line
*
* USAGE: nexftlow process_f_perpairs.nf --dir my_data/
*/

Channel
    .fromFilePairs(params.dir+"*_{1,2}.fastq.gz")
    .set { samples_ch }

process foo {
  input:
  set sampleId, file(reads) from samples_ch

  script:
  """
  echo your_command --sample $sampleId --reads $reads
  """
}
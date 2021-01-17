/*
* This script shows how to fetch all FASTQ file pairs from a directory
* passed in the command line
*
* USAGE: nexftlow process_f_perpairs.nf --dir my_data/
*/

Channel
    .fromFilePairs(params.dir+"*_{1,2}.fastq.gz", checkIfExists:true)
    .set { samples_ch }

Channel
    .fromFilePairs(params.dir+"*_{1,2}.fastq.gz", flat:true, checkIfExists:true)
    .set { samples_ch1 }

process foo {
  input:
  set sampleId, file(reads) from samples_ch

  script:
  """
  echo your_command --sample $sampleId --reads $reads
  """
}

process bar {
  /*
  This process will use the flat: true option
  so you can access each of the files in the pair
  independently
  */

  input:
  set sampleId, file(R1), file(R2) from samples_ch1

  script:
  """
  echo your_command --sample $sampleId -1 $R1 -2 $R2
  """
}

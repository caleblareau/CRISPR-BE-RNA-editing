library(dplyr)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(reshape2)
library(stringr)
library(data.table)
library(diffloop)
library(seqinr)

# Import annotation data only if it is missing from the environment
if(!exists("exome_hepg2")){
  source("00_helperdata.R")
}

# Import functions to facilitate data preparation
source("01_essential_functions.R")

do_all_CBE <- function(simple_name, full_bamreadcount_library_file, cell_type = "HEK293"){
  
  if(cell_type == "HEK293"){
    dna_variants=exome_hek293t
  }else if(cell_type == "HepG2"){
    dna_variants=exome_hepg2
  }
  
  # Create new variables based on input
  step1_br_filter_file <- paste0("bam-readcount/" , simple_name, "-HQcounts.tsv")
  
  # Step 1 - Determine universe of possible edits / non-edits from bam-readcount data
  call_bam_readcount_CBE(full_bamreadcount_library_file , step1_br_filter_file)
  
  # Step 2 - Determine universe of possible edits / non-edits from bam-readcount data
  new_fastas <- processFastaSample_CBE(simple_name, step1_br_filter_file, dna_variants, gtf_forward, gtf_reverse, pad = 50)
  
  # Step 3 - Annotate secondary structure for all files generated in step 2
  new_structure_files <- secondaryStructureLaunch(unlist(new_fastas), simple_name)
  
  # Step 4 - remove/compress for economy
  system(paste0("rm ", step1_br_filter_file))
  simple_name
}


if(FALSE){
  
  do_all_CBE(simple_name = "161B",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp7-hiseq/bam-readcount/161B.all_positions.txt.gz",
             cell_type = "HEK293")
  
  do_all_CBE(simple_name = "161F",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp7-hiseq/bam-readcount/161F.all_positions.txt.gz",
             cell_type = "HEK293")
  
  do_all_CBE(simple_name = "161G",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp7-hiseq/bam-readcount/161G.all_positions.txt.gz",
             cell_type = "HEK293")
  
  do_all_CBE(simple_name = "161H",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp7-hiseq/bam-readcount/161H.all_positions.txt.gz",
             cell_type = "HEK293")
  
}





if(TRUE){
  do_all_CBE(simple_name = "nova160A",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp6-novaseq/bam-readcount/160A.all_positions.txt.gz",
             cell_type = "HEK293")
  do_all_CBE(simple_name = "nova160B",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp6-novaseq/bam-readcount/160B.all_positions.txt.gz",
             cell_type = "HEK293")
  do_all_CBE(simple_name = "nova160C",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp6-novaseq/bam-readcount/160C.all_positions.txt.gz",
             cell_type = "HEK293")
  do_all_CBE(simple_name = "nova160D",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp6-novaseq/bam-readcount/160D.all_positions.txt.gz",
             cell_type = "HEK293")
  do_all_CBE(simple_name = "nova160E",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp6-novaseq/bam-readcount/160E.all_positions.txt.gz",
             cell_type = "HEK293")
  do_all_CBE(simple_name = "nova160F",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp6-novaseq/bam-readcount/160F.all_positions.txt.gz",
             cell_type = "HEK293")
  do_all_CBE(simple_name = "nova160G",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp6-novaseq/bam-readcount/160G.all_positions.txt.gz",
             cell_type = "HEK293")
  do_all_CBE(simple_name = "nova160H",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp6-novaseq/bam-readcount/160H.all_positions.txt.gz",
             cell_type = "HEK293")
}


if(FALSE){
  do_all_CBE(simple_name = "146A",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp4-hiseq/bam-readcount/146A.all_positions.txt.gz",
             cell_type = "HEK293")
  
  do_all_CBE(simple_name = "146B",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp4-hiseq/bam-readcount/146B.all_positions.txt.gz",
             cell_type = "HEK293")
  
  do_all_CBE(simple_name = "146C",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp4-hiseq/bam-readcount/146C.all_positions.txt.gz",
             cell_type = "HEK293")
  
  do_all_CBE(simple_name = "146D",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp4-hiseq/bam-readcount/146D.all_positions.txt.gz",
             cell_type = "HEK293")
}

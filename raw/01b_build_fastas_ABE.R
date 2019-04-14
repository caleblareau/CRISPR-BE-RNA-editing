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

do_all_ABE <- function(simple_name, full_bamreadcount_library_file, cell_type = "HEK293"){
  
  if(cell_type == "HEK293"){
    dna_variants=exome_hek293t
  }else if(cell_type == "HepG2"){
    dna_variants=exome_hepg2
  }
  
  # Create new variables based on input
  step1_br_filter_file <- paste0("bam-readcount/" , simple_name, "-HQcounts.tsv")
  
  # Step 1 - Determine universe of possible edits / non-edits from bam-readcount data
  call_bam_readcount_ABE(full_bamreadcount_library_file , step1_br_filter_file)
  
  # Step 2 - Determine universe of possible edits / non-edits from bam-readcount data
  new_fastas <- processFastaSample_ABE(simple_name, step1_br_filter_file, dna_variants, gtf_forward, gtf_reverse, pad = 50)
  
  # Step 3 - Annotate secondary structure for all files generated in step 2
  new_structure_files <- secondaryStructureLaunch(unlist(new_fastas), simple_name)
  
  # Step 4 - remove/compress for economy
  system(paste0("rm ", step1_br_filter_file))
  simple_name
}




if(FALSE){
  
  do_all_ABE(simple_name = "156A",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp7-hiseq/bam-readcount/156A.all_positions.txt.gz",
             cell_type = "HEK293"
  )
  
  do_all_ABE(simple_name = "156B",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp7-hiseq/bam-readcount/156B.all_positions.txt.gz",
             cell_type = "HEK293"
  )
  
  do_all_ABE(simple_name = "157A",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp7-hiseq/bam-readcount/157A.all_positions.txt.gz",
             cell_type = "HEK293"
  )
  
  do_all_ABE(simple_name = "157B",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp7-hiseq/bam-readcount/157B.all_positions.txt.gz",
             cell_type = "HEK293"
  )
  
  do_all_ABE(simple_name = "158A",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp7-hiseq/bam-readcount/158A.all_positions.txt.gz",
             cell_type = "HEK293"
  )
  
  do_all_ABE(simple_name = "158B",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp7-hiseq/bam-readcount/158B.all_positions.txt.gz",
             cell_type = "HEK293"
  )
}


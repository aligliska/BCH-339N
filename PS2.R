## BCH 339N Problem Set 2

## Alice Gee, ag67642
## For this assignment please only use base R and do not rely on any packages. 


## Q1 - 3pt
## GC-content of a piece of DNA is defined as the percentage of Gs + Cs. 
## For example "ACCTGCA" has a 57.1% GC-content
## Write a function that will take a text file in FASTA format as input. 
## FASTA format: A sequence record in a FASTA format consists of a single-line description (sequence name), followed by line(s) of sequence data. 
## The first character of the description line is a greater-than (">") symbol.

# Example input/output:
# The file "example_dna.fa" contains:

# >T7_promoter_1
# TAATACGACTCACTATAGGG
# >T7_promoter_2
# TAATACGACTCACTATAGGGG

# gc_calculator("example_dna.fa")
# [1] T7_promoter_1
# [1] 0.4

## The output should be the sequence id with the lowest GC content followed by the calculated value 
gc_calculator = function (input_fasta_name) { 
  path = setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # had to look this code up in the past, since I always have issues with my working directory
  fasta_file = readLines(paste(path, input_fasta_name, sep = "/"))
  sequence_names = c()
  sequences = c()
  for (i in seq(1, length(fasta_file), by = 2)){
    sequence_names = c(sequence_names, substr(fasta_file[i], 2, nchar(fasta_file[i])))
    sequences = c(sequences, toupper(fasta_file[i+1]))
  }
  gc_count = c()
  for (j in (1:length(sequences))){
    temp_seq = unlist(strsplit(sequences[j], split = ""))
    temp_count = 0
    for (k in 1:length(temp_seq)){
      if (temp_seq[k] == "G" | temp_seq[k] == "C") {
        temp_count = temp_count + 1
      }
    }
    gc_count[j] = round((temp_count / length(temp_seq)), digits = 3)
  }
  
  lowest = which.min(gc_count)
  for (r in (1: length(lowest))){
    out = sprintf("The sequence %s had the lowest GC_content at %s (i.e. %s).", 
                  sequence_names[r], gc_count[r], paste(gc_count[r]*100, "%", sep = ""))
  }
  return(out)
}

# function will work if the fasta file is in the same working directory as this file! 
## edit this variable match the title of input/test fasta file 
name_of_fasta = "example_dna.fa"
gc_calculator(name_of_fasta)

## Q2 - 6pt
## Every amino acid in a protein is encoded by three nucleotides. 
## Execute the following two lines to get a list of all codons and corresponding amino acids
codons = c('UUU','UUC','UUA','UUG','UCU','UCC','UCA','UCG','UAU','UAC','UAA','UAG','UGU','UGC','UGA','UGG','CUU','CUC','CUA','CUG','CCU','CCC','CCA','CCG','CAU','CAC','CAA','CAG','CGU','CGC','CGA','CGG','AUU','AUC','AUA','AUG','ACU','ACC','ACA','ACG','AAU','AAC','AAA','AAG','AGU','AGC','AGA','AGG','GUU','GUC','GUA','GUG','GCU','GCC','GCA','GCG','GAU','GAC','GAA','GAG','GGU','GGC','GGA','GGG')
amino_acids = c('F','F','L','L','S','S','S','S','Y','Y','*','*','C','C','*','W','L','L','L','L','P','P','P','P','H','H','Q','Q','R','R','R','R','I','I','I','M','T','T','T','T','N','N','K','K','S','S','R','R','V','V','V','V','A','A','A','A','D','D','E','E','G','G','G','G' )

## Write a function that will take a coding region sequence as input. You can assume the sequence is starting with AUG and is a multiple of three nucleotides. 
## The output should be the corresponding protein sequence. Report up to the first stop codon. 

# Example input/output:
# translate_rna_to_protein("AUGCUGGUGUAGUCGUGGUUAUUCUUU")
# [1] "MLV"

# "MLV*SWLFF" will be acceptable but
#  the protein sequence should ideally end at stop codons.
translate_rna_to_protein = function(rna) {
  rna = toupper(rna)
  sep_codons = c()
  for (i in seq(1, nchar(rna), by=3)){
    sep_codons = c(sep_codons, substr(rna, i, i+2))
  }
  
  aa_seq = c()
  for (j in 1:length(sep_codons)){
    aa_seq = c(aa_seq, amino_acids[which(codons == sep_codons[j])])
  }
  protein_seq = ""
  for (k in 1:length(aa_seq)){
    if (aa_seq[k] == "*"){
      break
    }else{
      protein_seq = paste0(protein_seq, aa_seq[k])
    }
  }
  return(protein_seq)
}

translate_rna_to_protein("AUGCUGGUGUAGUCGUGGUUAUUCUUU")

## Q3 - 3pt
## Let's define
## F(n) = F(n-1) + 3* F(n-2)
## F(1) = 1
## F(0) = 0
## Write a function that takes any positive integer k as input and 
## prints the value of smallest "n" such that F(n) >= k

# Example input/output:
# smallest_n_finder(20)
# [1] 6

smallest_n_finder = function ( k) { 
  if (k == 0){return(0)}
  else if (k == 1){return(1)}
  else{
    store = c(0, 1)
    n = length(store)
    while (store[length(store)] < k){
      n = n + 1
      temp = store[n-1] + 3*store[n-2]
      store = c(store, temp)
    }
  }
  return(n-1)
}
smallest_n_finder(0)
smallest_n_finder(1)
smallest_n_finder(20)


## Q4 - 9 pt
## In class, we have discussed dynamic programming for global alignment of two sequences. 
## Please implement this algorithm. 
## Specifically, the function should take two strings in addition to scores for match, mismatch, gap penalty
## Note that all mismatches will have the same score for this question. 
## We will use a linear gap penalty (no separate penalty for gap open). 
## The output should be the optimal alignment and the associated score

## Expected Output Example: 
## ATT-AGC
## ATTCAGG
## Score: 18 

alignment_function = function (str1, str2, match_score, mismatch_score, gap_penalty) { 
  str1 = unlist(strsplit(str1, split = ""))
  str2 = unlist(strsplit(str2, split = ""))
  
  str1 = c(0,str1)
  str2 = c(0, str2)
  
  scoring_mat <- matrix(NA, nrow = length(str1), ncol = length(str2))
  scoring_mat[,1] <- sapply(1:length(str1)-1, function(x) x * gap_penalty)
  scoring_mat[1,] <- sapply(1:length(str2)-1, function(x) x * gap_penalty)
  direction_mat <- matrix(NA, nrow = length(str1), ncol = length(str2))
  direction_mat[,1] <- "a"
  direction_mat[1,] <- "l"
  direction_mat[1,1] <- "*"
  rownames(scoring_mat) <- str1
  colnames(scoring_mat) <- str2
  rownames(direction_mat) <- str1
  colnames(direction_mat) <- str2
  
  for (i in 2:(length(str1))){
    for (j in 2:(length(str2))){
      if (str1[i] == str2[j]){
        scoring_mat[i,j] = scoring_mat[i-1, j-1] + match_score
        direction_mat[i,j] = "d"
      } else{
        mismatched = scoring_mat[i-1, j-1] + mismatch_score
        x_gap = scoring_mat[i, j-1] + gap_penalty
        y_gap = scoring_mat[i-1, j] + gap_penalty
        possible_value = c(mismatched, y_gap, x_gap)
        scoring_mat[i,j] = max(possible_value)
        if (mismatched == scoring_mat[i,j]){
          direction_mat[i,j] = "d"
        } else if(x_gap == scoring_mat[i,j]){
          direction_mat[i,j] = "l"
        } else if(y_gap == scoring_mat[i,j]){
          direction_mat[i,j] = "a"
        }
      }
    }
  }
  
  i = length(str1)
  j = length(str2)
  seq_align1 = c()
  seq_align2 = c()
  while (direction_mat[i,j] != "*"){
    if (direction_mat[i,j] == "d"){
      seq_align1 = c(str1[i], seq_align1)
      seq_align2 = c(str2[j], seq_align2)
      i = i - 1
      j = j - 1
    } else if (direction_mat[i,j] == "l"){
      seq_align1 = c("-", seq_align1)
      seq_align2 = c(str2[j], seq_align2)
      j = j - 1
    } else if (direction_mat[i,j] == "a"){
      seq_align1 = c(str1[i], seq_align1)
      seq_align2 = c("-", seq_align2)
      i = i - 1
    }
  }
  
  seq_align1 = paste(seq_align1, collapse = "")
  seq_align2 = paste(seq_align2, collapse = "")
  seq_align1 = substr(seq_align1, 1, nchar(seq_align1))
  seq_align2 = substr(seq_align2, 1, nchar(seq_align2))
  
  seq_align1 = paste(seq_align1, collapse = "")
  seq_align2 = paste(seq_align2, collapse = "")
  cat(seq_align1, sep = "\n")
  cat(seq_align2, sep = "\n")
  cat("Score: ")
  cat(scoring_mat[length(str1),length(str2)], sep = "\n")
}

str1 = "ATTAGC"
str2 = "ATTCAGG"
match_score = 6
mismatch_score = -4
gap_penalty = -8
alignment_function(str1, str2, match_score, mismatch_score, gap_penalty)







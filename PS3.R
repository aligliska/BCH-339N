## Alice Gee
## ag67642


## Q1 - 6pt
## Sequences for the 2019-nCoV are available through NCBI: https://www.ncbi.nlm.nih.gov/genbank/2019-ncov-seqs/
## Wuhan_nCov19.fasta (see Wuhan_nCov19-1.fasta  Download Wuhan_nCov19-1.fastaattached  Download attached) contains the complete sequence of one isolate of the 2019 Coronavirus
## An open reading frame (ORF) is a stretch of sequence that starts with a start codon (ATG) 
## and ends with a stop codon (TAG, TAA, TGA) without any other stop codons in between. 
## By enumerating ORFs, we can create a list of candidate proteins encoded by this coronavirus. 
## Write a function that will find all ORFs with at least 59 aa. 
## Your function should output start and stop position of all such ORFs. 

fasta_file = 'Wuhan_nCov19.fasta'
path = setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
ncov = readLines(paste(path, fasta_file, sep = "/"))

## We can collapse the multiline FASTA as follows: 
sequence = paste(ncov[2:length(ncov)], collapse= "" )

ncov_orf_finder = function (query_sequence, orf_length_threshold = 60) { 
  start_codon = "ATG"
  stop_codons = c("TAG", "TAA", "TGA")
  store = c()
  for (frame in 0:2){
    start_positions = c()
    stop_positions = c()
    for (i in seq(1 + frame,nchar(query_sequence), 3)){
      if (substr(query_sequence, i, i + 2) == start_codon){
        start_positions = c(start_positions, i)
      } else if (substr(query_sequence, i, i + 2) %in% stop_codons){
        stop_positions = c(stop_positions, i)
      }
    } 
    for (j in 1:length(stop_positions)){
      if (length(start_positions) == 0){
        break
      }
      if (stop_positions[j] < start_positions[1]){
        next
      }else {
        temp_start = start_positions[start_positions < stop_positions[j]]
        for (k in 1:length(temp_start)){
          if (stop_positions[j] - temp_start[k] >= orf_length_threshold*3){
            temp_string = paste(temp_start[k], stop_positions[j] + 2, sep = "-")
            store = c(store, temp_string)
            break
          }
        }
        start_positions = start_positions[!(start_positions %in% temp_start)]
      }
    }
  }
return(store)
}

output = ncov_orf_finder (sequence)
for (pair in output){
  print(pair)
}

## Example Output: 

## Note that the order doesn't need to be identical to the example output.
# [1] "28734-28955"
# [1] "28284-28577"
# [1] "27894-28259"
# [1] "26523-27191"
# [1] "21936-22199"
# [1] "10215-10400"
# [1] "6156-6350"
# [1] "2958-3206"
# [1] "28274-29533"
# [1] "21536-25384"
# [1] "15461-15667"
# [1] "266-13483"
# [1] "27394-27759"
# [1] "27202-27387"
# [1] "26245-26472"
# [1] "25393-26220"
# [1] "13768-21555"



## Q2 - 6pt
## A potential ORF we identified in Q1 spans the nucleotides "26523-27191"
substr(sequence, 26523, 27191)
## This is one of the proteins encoded by the new 2019 coronavirus genome. 
## Use blastn to search for similar sequences 
## You will notice that the first few hits match the Wuhan seafood market pneumonia virus isolate 2019.
## What other organisms did you find? 
## Click on at least one such alignment and explain in your own words the output. Include a description of 
## Score, Expect and Identities and include a screenshot. 

# I clicked on the "Bat coronavirus isolate BANAL-20-236/Laos/2020, complete genome." The output 
# shows the query sequence alignment against the target alignment (i.e. the specific sample).  
# This particular entry has a score of 1110 bits (577), which indicates that normalized raw score (in parentheses)
# of the sequence alignment score (i.e. that was generated using the respective matrix and gap penalties). 
# This bit score allows for comparison across different alignment scores calculated by other substitution
# matrices. The expect for this alignment was 0.0, which indicates that there is expected to be 0 other hits
# in the database (of size S) that would have a similar score by chance. The identity of this sequence alignment 
# was 631/658(96%). This value indicates how similar the query sequence was to the target sequence
# The value of 96% means that the target sequence and the query sequence are 96% similar. The overall sequence also
# has 0 gaps. 


## Q3 - 3pt
## We will use the same sequence for this question
substr(sequence, 26523, 27191)
## Use blastx to search the database. 
## Describe how blastx differs from blastn.
## Given the results page, what do you think is the most likely function of this protein? 

# Blastn searches through a nucleotide database with an input nucleotide sequence to query. 
# Blastx searches through a protein database by translating the input nucleotde query. 
# While both blastn and blastx take in an input query sequence, the output is different 
# in that blastn returns matches for other similar nucleotide sequences while blastx
# returns matches for similar proteins coded by the input nucleotide sequence. Given 
# the results page, the function of this protein is most likely a membrane glycoprotein, 
# which allows for cell-to-cell recognition and binding to other molecules. 


## Q4 -6pt 
## Enumerate k-mers
## BLAST uses heuristics to speed up the task of searching a query sequence against a large database. 
## For example, given a word size (default 3), BLAST will create a table of all possible short words (k-mers) contained in the query. 
## Write a function that will create this table of all possible words of given size. 
## For example, given  a a word size 3 and query sequence (LRITSLRI): 
## we will have LRI, RIT, ITS, TSL, SLR, LRI => Hence 5 distinct words are possible. 

test_sequence = "LRITSLRI"
test_sequence2 = "LRITSLRIK"

enumerate_words = function(query_seq, word_size = 3) { 
  n = nchar(query_seq)
  store = c()
  for (i in 1:(n-word_size + 1)){
    temp = substr(query_seq, i, i + 2)
    if (!temp %in% store){
      store = c(store, temp)
    }
  }
  return(length(store))
}

enumerate_words (test_sequence) 
## Expected Output
## [1] 5
enumerate_words (test_sequence2) 
## Expected Output
## [1] 6
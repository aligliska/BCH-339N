# Alice Gee
# ag67642

## BCH 339N
## PS4 - Hidden Markov Models

## We have covered the Viterbi and Forward Algorithms in class. 
## These methods are implemented in the R package - HMM

library(HMM)

## In this homework, we are going to develop an HMM that can predict secondary structure of proteins.
## Kaggle (https://www.kaggle.com) organizes data-science and machine learning challenges for the community
## They have many challenges related to bioinformatics. 
## We will explore the https://www.kaggle.com/alfrandom/protein-secondary-structure/data challenge. 
## Read the associated introduction to familiarize yourself with this dataset. 

## Q1 -3pt
## What is protein secondary structure? 
## We are going to work with the three state secondary structure encoding (sst3). 
## Explain the biochemical features of the three states that we will model (C, E, and H ).

### Protein secondary structure is the local three dimensional segments that make up the overall macromolecule. 
### These structures are largely known to be made of alpha helices and beta sheets, with less common 
### structures like beta turns and omega loops also contributing to secondary structure. 
### C (i.e. loops and irregular elements) help connect two secondary structures together in a protein. 
### They can also hold active site residues that help harbor enzyme reactions and ligand binding. 
### Moreover, protein loops allow polypeptide chains to change directions. E (β-strand) are made
### of nearly linear, but zigzagging polypeptide chains that can form hydrogen bonds with itself
### to then create beta sheets (which are typically key components to a protein's structure). H (α-helix) 
### is a structural motif in most proteins, characterized by its tight helical structure that is stabilized 
### by internal hydrogen bonds. α-helix are key components in determining a protein's overall structure, 
### which thus influences the protein's function. 

## Q2 - 3pt 
## We will use the data with the strict quality controls for this homework. Here  Download Hereis the file for convenience. 
## First, we will remove all examples with non-standard amino acids 
prot_sec_data = read.csv('./2018-06-06-pdb-intersect-pisces.csv', stringsAsFactors = F)
str(prot_sec_data)
prot_sec_data_std = prot_sec_data[as.logical(prot_sec_data$has_nonstd_aa) == F, ]

## Print the number of rows and columns of the remaining data. 
nrow(prot_sec_data_std)
ncol(prot_sec_data_std)

## Next, we are going to split our dataset into two parts. 
## The first part should contain 10% of the sequences and will be used for training our HMM
## The remaining 90% will be used as a test set to assess how well our HMM does in predicting secondary structure. 
## To ensure reproducibility of your code, we are going to use the set.seed function. 
set.seed(3)

train  = sample(nrow(prot_sec_data_std), size = floor(nrow(prot_sec_data_std)*.1 ) )
test =  setdiff(1:nrow(prot_sec_data_std), train  )

train_data = prot_sec_data_std[train, ]
test_data = prot_sec_data_std[test, ]


## Q3 - 6pt                                   
## We will use the training data to infer the parameters of the HMM model. 
## Specifically, use sst3 and seq columns of the data and determine the 
## transition, emission and initial probabilities.
## Write a function that will take the training data as input. 
## The output should be a list with names
## c(“initial_probs”, “emission_probs”, “transition_probs”)
## Note that emission_probs and transition_probs should be matrices of the form defined in initHMM. 
symbols = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y") 
states =  c("C", "E", "H")

get_hmm_probs<- function(train_data, symbols, states){
  # initial probabilities 
  init_counts = matrix(0, 1, length(states))
  colnames(init_counts) = states
  num = nrow(train_data)
  for (i in 1:num){
    start = substr(train_data$sst3[i], 1, 1)
    init_counts[,start] = init_counts[,start] + 1
  }
  init_counts = init_counts / sum(init_counts)
  
  trans_matrix = matrix(0, length(states), length(states), byrow = TRUE)
  rownames(trans_matrix) <- states
  colnames(trans_matrix) <- states
  
  emis_matrix = matrix(0, length(states), length(symbols), byrow = TRUE)
  rownames(emis_matrix) <- states
  colnames(emis_matrix) <- symbols
  
  # transition probabilities 
  states_seq = ""
  aa_seq = ""
  for (j in 1:num){
    states_seq = paste(states_seq, train_data$sst3[j], sep = "")
    aa_seq = paste(aa_seq, train_data$seq[j], sep = "")
  }
  for (k in 2:nchar(states_seq)){
    previous = substr(states_seq, k-1, k-1)
    current = substr(states_seq, k, k)
    trans_matrix[previous, current] = trans_matrix[previous, current] + 1
  }
  trans_matrix[1,] = trans_matrix[1,]/sum(trans_matrix[1,])
  trans_matrix[2,] = trans_matrix[2,]/sum(trans_matrix[2,])
  trans_matrix[3,] = trans_matrix[3,]/sum(trans_matrix[3,])
  
  for (k in 1:nchar(aa_seq)){
    curr_aa = substr(aa_seq, k, k)
    curr_state = substr(states_seq, k, k)
    emis_matrix[curr_state, curr_aa] = emis_matrix[curr_state, curr_aa] + 1
  }
  
  emis_matrix[1,] = emis_matrix[1,] / sum(emis_matrix[1,])
  emis_matrix[2,] = emis_matrix[2,] / sum(emis_matrix[2,])
  emis_matrix[3,] = emis_matrix[3,] / sum(emis_matrix[3,])
  
  return(list(initial_probs=init_counts, emission_probs=emis_matrix, transition_probs=trans_matrix)) 
}

hmm_values =  get_hmm_probs(train_data, symbols, states)
hmm_values

## The expected output should be as follows. 
## Do not worry if yours is slightly different. 
# $initial_probs
# [1] 1 0 0
# 
# $emission_probs
# [,1]       [,2]       [,3]       [,4]       [,5]       [,6]       [,7]
# [1,] 0.06585424 0.01165331 0.07881398 0.05734681 0.03039266 0.11675254 0.03626634
# [2,] 0.06147861 0.01684346 0.03523090 0.04828457 0.05855106 0.04896633 0.02452327
# [3,] 0.11237240 0.01019533 0.05272974 0.09006702 0.04173868 0.03504955 0.02049012
# [,8]       [,9]      [,10]      [,11]      [,12]      [,13]      [,14]
# [1,] 0.03188721 0.05509976 0.06162143 0.02146717 0.05971928 0.07833321 0.03431193
# [2,] 0.09332077 0.04435443 0.10398829 0.02221732 0.02648834 0.02263841 0.02845341
# [3,] 0.05793930 0.06470303 0.11973293 0.02671922 0.03325915 0.02469259 0.04735854
# [,15]      [,16]      [,17]      [,18]      [,19]      [,20]
# [1,] 0.04505597 0.07720446 0.05752448 0.04383315 0.01026327 0.02659880
# [2,] 0.04784344 0.05181368 0.06474705 0.12821078 0.01862806 0.05341782
# [3,] 0.05923237 0.04983277 0.04145271 0.05986647 0.01582762 0.03674048
# 
# $transition_probs
# [,1]        [,2]       [,3]
# [1,] 0.80925703 0.110084193 0.08065877
# [2,] 0.20296365 0.783140502 0.01389585
# [3,] 0.09893198 0.004737097 0.89633093


## Q4 - 3pt
## Look at how hmms are defined in the HMM package
## We will Use the inferred parameters from Q3 to define an HMM. 
sec_struct_hmm = initHMM(States = states, Symbols = symbols,
                         startProbs = hmm_values$initial_probs, 
                         transProbs = hmm_values$transition_probs, 
                         emissionProbs = hmm_values$emission_probs)


## Next, We are going to assess the performance on the test data
## For each example in the test data, use the given sequence (seq column) to predict the most likely path of hidden states. 
## You can use the appropriate function in the hmm package for this step. 
## For each example, compare the predicted most likely hidden state path with the experimentally identified values (sst3 column). 
## Output a vector named percent_correct containing the percentage of amino acids whose secondary structure was correctly predicted for each example in test data.

percent_correct = c()
for (i in 1:nrow(test_data)){
  test_seq = unlist(strsplit(test_data[i, "seq"], ""))
  test_sst3 = unlist(strsplit(test_data[i, "sst3"], ""))
  viterbi_out = viterbi(sec_struct_hmm, test_seq)
  correct_ratio = mean(test_sst3 == viterbi_out)
  percent_correct = c(percent_correct, correct_ratio)
}
print(percent_correct)

## Q5 - 3pt
## Plot the distribution of percent_correct 
hist(percent_correct, 40, main = "Distribution of Percent Correct", xlab = "Percent Correct (as proportion)")

## What is the mean, median and 90th percentile of this distribution? 
print(paste("Mean:",mean(percent_correct)))
print(paste("Median:",median(percent_correct)))
print(paste("90th percentile:", quantile(percent_correct, probs = 0.9)))

# Write a sentence about your interpretation of these results.

## The mean percent correct is 0.456, indicating the average proportion of correct predictions based from 
## the generated HMM model compared to the actual sequence (across all examples). The median percent correct is
## 0.446, indicating the middle value (i.e. center value of the data set) of proportions of correct predictions based from the generated HMM
## model compared to the actual sequence (from all the examples). The 90th percentile is 0.709, indicating that 
## the percent correct of 0.709 is greater than 90% of all of the percent correct scores for the entire set of 
## examples. 

## Q6 - 3pt 
## Identify the training examples where we succeed less than 1%. Output the pdb_ids of these examples
## Examine the predicted hidden state for these and the actual secondary structure. 
## Notice that these test examples have a high percentage of "E" in their secondary structure. 
## Calculate for each test example, the percentage of residues with "E" secondary structure state.
## Plot the relationship between this value and percent_correct.
## What is the correlation coefficient? 

ids = which(percent_correct < 0.01)
pdb_ids_less = test_data$pdb_id[ids]
print(pdb_ids_less)

percent_e_less = c()
for (i in ids){
  test_sst3 = unlist(strsplit(test_data[i, "sst3"], ""))
  e_correct = mean(test_sst3 == "E")
  percent_e_less = c(percent_e_less, e_correct)
}
print(percent_e_less)

percent_e = c()
for (j in 1:nrow(test_data)){
  test_sst3 = unlist(strsplit(test_data[j, "sst3"], ""))
  e_correct = mean(test_sst3 == "E")
  percent_e = c(percent_e, e_correct)
}
print(percent_e)

plot(percent_e, percent_correct)
print(paste("correlation coefficient:", cor(percent_e, percent_correct)))

## Explain why you think our HMM doesn't work for this certain class of proteins?
## Hint: Think about the biology of the "E" state and the definition of a Markov chain

### The "E" state refers to beta strands. HMM does not work accurately for this class of proteins because of the 
### structure of beta strands. Beta strands are components that make up beta sheets, which are made up of long-range
### associations. This non-local component of beta strands and beta sheets make it difficult to predict the associations
### that would indicate the secondary structure is a beta sheet from the primary structure. 

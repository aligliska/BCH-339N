## BCH 339N Problem Set 1

## Q1 -3pt
## Write a function that will simulate the number of tails in a series of coin tosses.
## The input to the function is the number of coin tosses
## The output should be a single number equivalent to the number of tails in the series

number_of_tails = function (coin_tosses) {
  num_tails <- sample(c(0,1), replace = TRUE, size = coin_tosses, prob = c(0.5, 0.5))
  return (sum(num_tails))
}

## Q2 -3pt
## Using the function you wrote in Q2, generate 5000 experiments each with 40 coin tosses
## Plot the distribution of the outputs as a histogram
store <- c()
for (i in 1:5000)
  store[i] <- number_of_tails(40)
hist(store, main = "Distribution of # of Tails over 5000 Experiments", 
     xlab = "Number of Tails")


## Q3 -3pt
## Values and objects have different data types in R
## You are given two DNA sequences of equal length
## find out how many bases are different between the two

diff_bases = function(seq1, seq2) {
  s1 = unlist(strsplit(seq1, split = ""))
  s2 = unlist(strsplit(seq2, split = ""))
  n = length(s1)
  diff_count = 0
  for (i in 1:n){
    if (s1[i] != s2[i]){
      diff_count = diff_count + 1
    } 
  }
  return(diff_count)
}

## Q4 -3pt
## A purine to purine or a pyrimidine to pyrimidine change is called a transition
## Transversion is a change from a pyrimidine to a purine or vice versa.
## As in Q3, you are given two DNA sequences of equal length
## Write a function that will output the ratio of transversions to transitions

transversion_transition_ratio = function (dna1, dna2) {
  dna1 = toupper(dna1)
  dna2 = toupper(dna2)
  purine = c("A", "G")
  pyrimidine = c("C", "T")
  s1 = unlist(strsplit(dna1, split = ""))
  s2 = unlist(strsplit(dna2, split = ""))
  n = length(s1)
  transitions = 0
  transversions = 0
  for (i in 1:n){
    if ((s1[i] != s2[i]) & (s1[i] %in% purine) & (s2[i] %in% purine)) {
      transitions = transitions + 1
    } else if ((s1[i] != s2[i]) & (s1[i] %in% pyrimidine) & (s2[i] %in% pyrimidine)){
      transitions = transitions + 1
    } else if ((s1[i] == s2[i])){
      next
    } else {
      transversions = transversions + 1
    }
  }
  output = paste(transversions, transitions, sep = ":" )
  return(output)
}


## Q5 -3pt
## You are given a matrix of unknown dimensions.
## Write a function that will replace all zeros in the matrix with 1s.

## Example solution
replace_zeros = function(matrix) {
  num_rows = dim(matrix)[1]
  num_cols = dim(matrix)[2]
  for (row in 1:num_rows){
    for (col in 1:num_cols){
      if (matrix[row,col] == 0){
        matrix[row,col] = 1
      }
    }
  }
  return(matrix)
}

## Q6 -3pt
## Download the dslabs package (https://cran.r-project.org/web/packages/dslabs/index.html)
## execute "data(brca)" to load the BRCA data into your session
## execute "str(brca)
## Write a few sentence to describe what the output means in your own words as a comment

## Example solution
#install.packages("dslabs")
library(dslabs)
data(brca)
str(brca)
## answer: The str() function details the structure of the variable and outputs the data [from brca] in 
##        a more compact and readable format. The first line indicates the number of rows and columns,
##        and previews the first 5 entries of the first column. The second line tells us that the data
##        contains 2 different lists. The third and fourth lines are the row-names and column-names,
##        respectively. The fifth line gives details on the second list, that contains two types of 
##        values: "B" and "M". 

## Q7 -3pt
## This is an open ended question so many answers are possible.
## Look at the help page related to brca using "?brca"
## Explore the relationship between biopsy features and whether the tissue is malignant or benign
## Your answer will receive full points as long as it has at least one plot.

# Example solution
?brca
## answer: The list "x" contains 30 different columns/variables that are used to predict whether a person 
##        has breast tissue contains malignant or benign cancer. Based on the mean texture, malignant
##        tissue has an average texture rating higher than benign tissue. Likewise, malignant tissue
##        has a larger average mean than benign tissue. Malignant tissue also tend to be more smooth, 
##        on average, compared to benign tissue, though there is significantly more overlap than the
##        previous two variables. There is no distinct difference between the malignant and benign
##        tissue in terms of fractal dim, where there exists significant overlap in values that it 
##        would not be an accurate predictor variable on its own. 

store = data.frame(condition = brca$y, texture_mean = brca$x[,"texture_mean"], 
                   area_mean = brca$x[,"area_mean"], smoothness_mean = brca$x[,"smoothness_mean"],
                   fractal_dim_mean = brca$x[,"fractal_dim_mean"])
plot(store[,"condition"], store[,"texture_mean"], main = "Condition vs. Texture", 
     xlab = "Condition", ylab = "Texture Mean")
plot(store[,"condition"], store[,"area_mean"], main = "Condition vs. Area", 
     xlab = "Condition", ylab = "Area Mean")
plot(store[,"condition"], store[,"smoothness_mean"], main = "Condition vs. Smoothness", 
     xlab = "Condition", ylab = "Smoothness Mean")
plot(store[,"condition"], store[,"fractal_dim_mean"], main = "Condition vs. Fractal Dim", 
     xlab = "Condition", ylab = "Fractal Dim Mean")



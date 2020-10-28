#!/usr/bin/env Rscript

######################################################################
# This script will read an xml file from BEAST2 using MASCOT (unclear
# if it will work with other models) and shuffle the trait levels 
# assigned to each tip keeping their numbers. Usage:
#
# Rscript traitPermutator.R mascot.xml
#
# Written by G. Yebra, 28/10/20
######################################################################

# Import the name of the xml file
args = commandArgs(trailingOnly=TRUE)

inputname <- args[1]

# Read the lines in the input file
lines <- readLines(inputname)

# Extract the line containing the traits
traits_index <- grep("<typeTrait",lines)
traits_line <- lines[traits_index]

# Access the traits
traits_line_split1 <- unlist(strsplit(traits_line, " "))[17]
traits_line_split2 <- unlist(strsplit(traits_line_split1, '"'))[2]
traits_line_split3 <- unlist(strsplit(traits_line_split2, ","))

# loop over traits to extract each ID and its corresponding trait
ids <- vector(mode = "list", length = length(traits_line_split3))
original_traits <- vector(mode = "list", length = length(traits_line_split3))
for (i in 1:length(traits_line_split3)){
  ids[i] <- unlist(strsplit(traits_line_split3[i], "="))[1]
  original_traits[i] <- unlist(strsplit(traits_line_split3[i], "="))[2]
}

# shuffle traits
rand <- sample(original_traits)

# assign back to ids
df <- as.data.frame(cbind(ids, rand))
new_traits <- paste(df$ids, df$rand, sep="=")
new_traits <- paste(new_traits, collapse = ",")

# reconstruct again the line containing the traits
traits_line_merge1 <- paste(unlist(strsplit(traits_line_split1, '"'))[1], 
                            new_traits, 
                            unlist(strsplit(traits_line_split1, '"'))[3], 
                            sep = '"')

traits_line_merge2 <- paste(paste(unlist(strsplit(traits_line, " "))[1:16],collapse =" "),
                            traits_line_merge1,
                            sep = " ")
traits_line_merge3 <- noquote(traits_line_merge2)

# substitute the modified traits and create the new xml file
lines[traits_index] <- traits_line_merge3
fileConn <- file(paste(gsub("\\.xml", "", inputname), "_randomTrait.xml", sep = ""))
writeLines(lines, fileConn)
close(fileConn)

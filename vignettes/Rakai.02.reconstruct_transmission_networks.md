# Reconstruct transmission networks

## Introduction
[In the previous tutorial](Rakai.01.run_phyloscanner.md), we performed *phyloscanner* analyses to reconstruct phylogenetic relationships between individuals that together form a potential transmission network. The analysis proceeded in batches, in which each transmission network was given a batch number, and small networks were combined into a single batch. Each batch was identified by the prefix `ptyrX_` where `X` is the batch number. 

**In this tutorial, we will take a look at the output.** We will:
1. Collect pairs of individuals between whom phylogenetic linkage cannot be excluded based on the distance and topological relationship of viral reads from both individuals. 
2. Reconstruct transmission networks from these pairs. Unlike typical phylogenetic clusters, the transmission networks that we reconstruct from deep-sequence data contain information into the direction of transmission. 
3. Calculate the proportion of male-female, male-male, and female-female links in the most likely transmission chains. 

## Setting up the analysis
copy pairwise_relationships.rda to new folder
 
## Analysis
we now combine results from analysis of all potential transmission networks
keep all pairs of individuals between whom linkage cannot be excluded

get linked pairs
plot phyloscans for some pairs
get networks
plot networks
plot ML chains



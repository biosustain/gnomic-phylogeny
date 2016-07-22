Gnomic-phylogeny is an NPM package for generating distance matrix or phylogenetic trees based on [Gnomic](https://www.npmjs.com/package/gnomic-grammar) - grammar for describing genotypes and phenotypes of microbial strains.
# Installation

To install Neighbor-joining package with NPM use: `npm install neighbor-joining`

# Description

Gnomic-phylogeny operates on gnomic `Genotype` objects. It has two main functions:
* `buildDistanceMatrix(gnomicGenotypes)` - generates two dimensional array containing distances between genotypes
* `buildPhylogeneticTree(taxa, gnomicGenotypes, newick=false)` - generates a phylogenetic tree as an object or a string in Newick format. First, the function calls `buildDistanceMatrix` function and then it uses [Neighbour-joining](https://www.npmjs.com/package/neighbor-joining) NPM library to create a phylogenetic tree. To get more

Currently, features and plasmids are distinguished only by looking at their `name` property. The following table presents the approach used to calculate the distance between two genotypes.
|           | Genotype 1 | Genotype 2 | Distance |
|-----------|:----------:|:----------:|----------|
| feature A |  inserted  |   deleted  |     2    |
|           |  inserted  |     N/A    |     1    |
|           |     N/A    |  inserted  |     1    |
|           |     N/A    |     N/A    |     0    |
|           |  inserted  |  inserted  |     0    |
|           |   deleted  |   deleted  |     0    |
| plasmid A |  inserted  |     N/A    |     1    |
|           |     N/A    |  inserted  |     1    |
|           |     N/A    |     N/A    |     0    |
|           |  inserted  |  inserted  |     0    |
# Gnomic-phylogeny

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
| plasmid B |  inserted  |     N/A    |     1    |
|           |     N/A    |  inserted  |     1    |
|           |     N/A    |     N/A    |     0    |
|           |  inserted  |  inserted  |     0    |

# Example
```javascript
var Genotype = require("gnomic-grammar").Genotype;
var _gnomicPhylogeny = require("gnomic-phylogeny");
var buildDistanceMatrix = _gnomicPhylogeny.buildDistanceMatrix;
var buildPhylogeneticTree = _gnomicPhylogeny.buildDistanceMatrix;
var taxa = [
            { name: "A",
              genotype: "p1{geneA geneB}::cc+ geneC>geneD" },
            { name: "B",
              genotype: "geneC>geneD" },
            { name: "C",
              genotype: "p1{geneA geneF} geneC>geneE" },
            { name: "D",
              genotype: "p2{geneG geneC} +geneG" }
        ];
var gnomicGenotypes = taxa.map(d => d.genotype).map(Genotype.parse);
var D = buildDistanceMatrix(gnomicGenotypes);
var treeObject = buildPhylogeneticTree(taxa, gnomicGenotypes);
var treeNewick = buildPhylogeneticTree(taxa, gnomicGenotypes, true);
```
As a result, `D`, `treeObject` and `treeNewick` will keep the following information:

`D`:
```javascript
[[0, 1, 2, 5],
 [1, 0, 3 ,4],
 [2, 3, 0, 5],
 [5, 4, 5, 0]]
```
`treeObject`:
```javascript
{
    "taxon": null,
    "length": null,
    "children": [{
        "taxon": null,
        "length": 1.75,
        "children": [{
            "taxon": null,
            "length": 0.5,
            "children": [{
                "taxon": {
                    "name": "A",
                    "genotype": "p1{geneA geneB}::cc+ geneC>geneD"
                },
                "length": 0.5,
                "children": []
            }, {
                "taxon": {
                    "name": "B",
                    "genotype": "geneC>geneD"
                },
                "length": 0.5,
                "children": []
            }]
        }, {
            "taxon": {
                "name": "C",
                "genotype": "p1{geneA geneF} geneC>geneE"
            },
            "length": 1.5,
            "children": []
        }]
    }, {
        "taxon": {
            "name": "D",
            "genotype": "p2{geneG geneC} +geneG"
        },
        "length": 1.75,
        "children": []
    }]
}
```
`treeNewick`:
```javascript
"(((A:0.5,B:0.5):0.5,C:1.5):1.75,D:1.75);"
```
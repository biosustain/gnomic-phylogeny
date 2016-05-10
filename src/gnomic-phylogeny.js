import { RapidNeighborJoining, allocateSquareMatrix } from 'neighbor-joining';
import {Mutation, Plasmid} from 'gnomic-grammar';

export function buildPhylogeneticTree(taxa, gnomicGenotypes, newick=false) {
    let D = buildDistanceMatrix(gnomicGenotypes),
        RNJ = new RapidNeighborJoining(D, taxa);

    RNJ.run();
    return newick ? RNJ.getAsNewick() : RNJ.getAsObject();
}

export function buildDistanceMatrix(gnomicGenotypes) {
    let genotypes = [],
        allFeatures = {},
        localFeatures,
        genotype1,
        genotype2,
        dist,
        k = 0,
        n = gnomicGenotypes.length,
        n1 = n - 1;

    for(let genotype of gnomicGenotypes) {
        let changes = genotype.changes();
        localFeatures = [];
        for(let change of changes) {
            if (change instanceof Plasmid) {
                if (allFeatures[change.name] === undefined) {
                    allFeatures[change.name] = k++;
                }
                localFeatures.push(allFeatures[change.name]);
            }
            else if (change instanceof Mutation) {
                let del = change.before ? "del#" : "",
                    featureTree = change.before || change.after;
                for (let item of featureTree.features()) {
                    let name = del + item.name + (item.variant || "");
                    if (allFeatures[name] === undefined) {
                        allFeatures[name] = k++;
                    }
                    localFeatures.push(allFeatures[name]);
                }
            }
        }
        localFeatures.sort((a, b) => a - b);
        genotypes.push(localFeatures);
    }
    let D = allocateSquareMatrix(n, 0);
    for (let i = 0; i < n1; i++) {
        genotype1 = genotypes[i];
        for (let j = i + 1; j < n; j++) {
            genotype2 = genotypes[j];
            dist = calculateDistance(genotype1, genotype2);
            D[i][j] = D[j][i] = dist;
        }
    }
    return D;
}

function calculateDistance(seq1, seq2) {
    var i1 = 0, i2 = 0,
        val1, val2,
        len1 = seq1.length, len2 = seq2.length,
        totalLength = len1 + len2,
        dist = 0;

    if (!len1) return len2;
    else if(!len2) return len1;

    while(i1 + i2 < totalLength) {
        val1 = seq1[i1];
        val2 = seq2[i2];

        if (val1 !== val2) dist++;

        if (val1 < val2) i1++;
        else if (val1 > val2) i2++;
        else { i1++; i2++; }

        if (i1 === len1) {
            dist += len2 - i2;
            break;
        }
        if (i2 === len2) {
            dist += len1 - i1;
            break;
        }
    }
    return dist;
}

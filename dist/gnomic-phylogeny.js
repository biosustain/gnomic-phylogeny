"use strict";

Object.defineProperty(exports, "__esModule", {
    value: true
});
exports.buildPhylogeneticTree = buildPhylogeneticTree;
exports.buildDistanceMatrix = buildDistanceMatrix;

var _neighborJoining = require("neighbor-joining");

function buildPhylogeneticTree(taxa, gnomicGenotypes) {
    var newick = arguments.length <= 2 || arguments[2] === undefined ? false : arguments[2];

    var D = buildDistanceMatrix(gnomicGenotypes),
        RNJ = new RapidNeighbourJoining(D, taxa);

    RNJ.run();
    return newick ? RNJ.getAsNewick() : RNJ.getAsObject();
}

function buildDistanceMatrix(gnomicGenotypes) {
    var genotypes = [],
        allFeatures = {},
        localFeatures = void 0,
        genotype1 = void 0,
        genotype2 = void 0,
        dist = void 0,
        k = 0,
        n = gnomicGenotypes.length,
        n1 = n - 1;

    var _iteratorNormalCompletion = true;
    var _didIteratorError = false;
    var _iteratorError = undefined;

    try {
        for (var _iterator = gnomicGenotypes[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
            var genotype = _step.value;

            var changes = genotype.changes();
            localFeatures = [];
            var _iteratorNormalCompletion2 = true;
            var _didIteratorError2 = false;
            var _iteratorError2 = undefined;

            try {
                for (var _iterator2 = changes[Symbol.iterator](), _step2; !(_iteratorNormalCompletion2 = (_step2 = _iterator2.next()).done); _iteratorNormalCompletion2 = true) {
                    var change = _step2.value;

                    if (change instanceof Plasmid) {
                        if (allFeatures[change.name] === undefined) {
                            allFeatures[change.name] = k++;
                        }
                        localFeatures.push(allFeatures[change.name]);
                    } else if (change instanceof Mutation) {
                        var del = change.before ? "del#" : "",
                            featureTree = change.before || change.after;
                        var _iteratorNormalCompletion3 = true;
                        var _didIteratorError3 = false;
                        var _iteratorError3 = undefined;

                        try {
                            for (var _iterator3 = featureTree.features()[Symbol.iterator](), _step3; !(_iteratorNormalCompletion3 = (_step3 = _iterator3.next()).done); _iteratorNormalCompletion3 = true) {
                                var item = _step3.value;

                                var name = del + item.name + (item.variant || "");
                                if (allFeatures[name] === undefined) {
                                    allFeatures[name] = k++;
                                }
                                localFeatures.push(allFeatures[name]);
                            }
                        } catch (err) {
                            _didIteratorError3 = true;
                            _iteratorError3 = err;
                        } finally {
                            try {
                                if (!_iteratorNormalCompletion3 && _iterator3.return) {
                                    _iterator3.return();
                                }
                            } finally {
                                if (_didIteratorError3) {
                                    throw _iteratorError3;
                                }
                            }
                        }
                    }
                }
            } catch (err) {
                _didIteratorError2 = true;
                _iteratorError2 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion2 && _iterator2.return) {
                        _iterator2.return();
                    }
                } finally {
                    if (_didIteratorError2) {
                        throw _iteratorError2;
                    }
                }
            }

            localFeatures.sort(function (a, b) {
                return a - b;
            });
            genotypes.push(localFeatures);
        }
    } catch (err) {
        _didIteratorError = true;
        _iteratorError = err;
    } finally {
        try {
            if (!_iteratorNormalCompletion && _iterator.return) {
                _iterator.return();
            }
        } finally {
            if (_didIteratorError) {
                throw _iteratorError;
            }
        }
    }

    var D = (0, _neighborJoining.allocateSquareMatrix)(n, 0);
    for (var i = 0; i < n1; i++) {
        genotype1 = genotypes[i];
        for (var j = i + 1; j < n; j++) {
            genotype2 = genotypes[j];
            dist = calculateDistance(genotype1, genotype2);
            D[i][j] = D[j][i] = dist;
        }
    }
    return D;
}

function calculateDistance(seq1, seq2) {
    var i1 = 0,
        i2 = 0,
        val1,
        val2,
        len1 = seq1.length,
        len2 = seq2.length,
        totalLength = len1 + len2,
        dist = 0;

    if (!len1) return len2;else if (!len2) return len1;

    while (i1 + i2 < totalLength) {
        val1 = seq1[i1];
        val2 = seq2[i2];

        if (val1 !== val2) dist++;

        if (val1 < val2) i1++;else if (val1 > val2) i2++;else {
            i1++;i2++;
        }

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
//- Create controls for the power iterations

//- Residual on the eigenvalue
scalar lambdaRes(1.);
scalar psiRes(1.);
scalar psiAdjRes(1.);
scalar XRes(1.);
label currPowerIter(0);
label currCellIter(0);
bool powerConverged(false);
bool cellConverged(false);
//scalar AitkenRelax(1.0);
dictionary cellDict(mesh.solutionDict().subDict("cellProblem"));

//- Tolerance for the eigenvalue
scalar eigenValTol
(
    readScalar(cellDict.lookup("eigenvalueTol"))
);

//- Tolerance for the eigenfunction
scalar eigenFunTol
(
    readScalar(cellDict.lookup("eigenfunctionTol"))
);

//- Tolerance for the first order corrector
scalar corrTol
(
    readScalar(cellDict.lookup("cellProblemTol"))
);

//- Non orthogonal corrections
label nNonOrthCorr
(
    readLabel(cellDict.lookup("nNonOrthogonalCorrectors"))
);

//- Max iterations
label maxPowerIter
(
    readLabel(cellDict.lookup("maxPowerIterations"))
);

label maxCellProbIter
(
    readLabel(cellDict.lookup("maxCellProblemIterations"))
);

vector Xref
(
    cellDict.lookup("Xref")
);

Switch setXref
(
    cellDict.lookup("setXref")
);

//label XrefCell
//(
//    readLabel(cellDict.lookup("XrefCell"))
//);

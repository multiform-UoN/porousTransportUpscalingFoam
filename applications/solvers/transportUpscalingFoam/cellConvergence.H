{
    currCellIter++;

    cellConverged =
        (
            XRes < corrTol
        )
        ||
        (
            currCellIter > maxCellProbIter
        );

    Info<<"\nResidual : " << XRes << " Converged : " << cellConverged
        <<endl;

}

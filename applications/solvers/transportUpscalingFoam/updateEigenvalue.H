{
    dimensionedScalar one("one",dimless/dimTime,1.);
    dimensionedScalar smallV
    (
        "smallV",
        psi.dimensions()*psi.dimensions()*dimVol,
        small
    );

    scalar lambdaOld = lambda.value();

    //- The new Eigenvalues are calculated using the Rayleight quotient
    //  such that lambda corresponds to the spectral radius of the ADE operator
    lambda = lambda *
        (
            fvc::domainIntegrate(psi.prevIter()*psi)
            /
            ( fvc::domainIntegrate(psi*psi) + smallV)
        );

    lambdaAdj = lambdaAdj *
        (
            fvc::domainIntegrate(psiAdj.prevIter()*psiAdj)
            /
            (fvc::domainIntegrate(psiAdj*psiAdj) + smallV)
        );

    //- Calculate new residual
    lambdaRes = mag(lambdaOld - lambda.value())/mag(lambdaOld);

}
    //- Evaluate error on eigenvalues
    dimensionedScalar eigenErr = mag(lambdaAdj -lambda)/(mag(lambda));

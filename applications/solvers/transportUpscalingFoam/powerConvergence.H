// - Check convergence
{
    currPowerIter++;

    powerConverged =
        (
            (
                lambdaRes < eigenValTol
                //||
                //Required for convergence with Neumann BCs
            //    lambda.value() < scalar(1e-6)
            )
          &&
            (
             psiRes < eigenFunTol
            &&
              psiAdjRes < eigenFunTol
        )   )
        ||
        (
            currPowerIter > maxPowerIter
        );

    Info<<"\nPower iteration " << currPowerIter << endl
        <<"Residual on eigenvalue : " << lambdaRes << " "
        <<"Residual on psi : " << psiRes << " "
        <<"Residual on psiAdj : " << psiAdjRes << endl
        << "Converged : " << powerConverged << endl;

    Info<< "Eigenvalue from psi : " << lambda.value() <<  endl;
    Info<< "Eigenvalue from psiAdj : " << lambdaAdj.value() <<  endl;
    Info<< "Relative discrepancy in the eigenvalues : " << eigenErr.value()
        << endl;
    Info<< "Average eigenvalue : " << lambda.value() <<  endl;


}


    //- Equation for the adjoint eigenfunction psi
    fvScalarMatrix psiEqnAdj
    (
        fvm::div(phi,psiAdj)
      + fvm::laplacian(D, psiAdj)
      ==
      - lambdaAdj*psiAdj
    );

    psiEqnAdj.relax();
    psiEqnAdj.solve();

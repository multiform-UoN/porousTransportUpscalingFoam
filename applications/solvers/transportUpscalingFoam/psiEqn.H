
    //- Equation for the eigenfunction psi
    fvScalarMatrix psiEqn
    (
        fvm::div(phi,psi)
      - fvm::laplacian(D, psi)
      ==
        lambda*psi
    );

    psiEqn.relax();
    psiEqn.solve();

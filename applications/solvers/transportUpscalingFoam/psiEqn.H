
    //- Equation for the eigenfunction psi
    fvScalarMatrix psiEqn
    (
        fvm::div(phi,psi)
      - fvm::laplacian(D, psi)
    );

    psiEqn.relax();
    if(computeResolvent)
    {
        psiEqn.solve();
    }
    else
    {
        solve(psiEqn == lambda*psi.prevIter());
    }

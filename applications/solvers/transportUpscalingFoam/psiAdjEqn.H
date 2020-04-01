
    //- Equation for the adjoint eigenfunction psi
    fvScalarMatrix psiAdjEqn
    (
        fvm::div(phiAdj,psiAdj)
      - fvm::laplacian(D, psiAdj)
      - fvm::Sp(fvc::div(phiAdj),psiAdj)

    );

    psiAdjEqn.relax();
    if(computeResolvent)
    {
        psiAdjEqn.solve();
    }
    else
    {
        solve(psiAdjEqn == lambdaAdj*psiAdj.prevIter());
    }

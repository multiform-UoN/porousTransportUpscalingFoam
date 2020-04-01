fvVectorMatrix XEqn
(
    fvm::div(phistar, X)
//  - fvm::Sp(fvc::div(phiplus),X)  //Compensates continuity errors
  - fvm::laplacian(beta*D, X)
);

//XEqn.setReference(XRefCell, XRefValue);
XEqn.relax();
solve( XEqn == sourceX );
X.correctBoundaryConditions();

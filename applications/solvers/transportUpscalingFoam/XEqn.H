fvVectorMatrix XEqn
(
     fvm::div(phiplus, X)
   - fvm::laplacian(beta*D, X)
     ==
     sourceX
);
//XEqn.setReference(XrefCell,Xref,true);
XEqn.relax();
XEqn.solve();

#include "rescaleX.H"
dimensionedScalar smallDim("sDim",X.dimensions(),small);

//- Subtract average value
volScalarField XFieldRes
(
    mag(X -X.prevIter())/(mag(X.prevIter())+smallDim)
);

XRes = fvc::domainIntegrate(XFieldRes).value()/vol.value();

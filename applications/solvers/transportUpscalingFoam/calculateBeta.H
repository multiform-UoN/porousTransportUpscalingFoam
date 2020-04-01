//- Rescale psi and psiAdj to have a normalised beta
psi *=  (vol / fvc::domainIntegrate(psi));

psiAdj *=  (vol / fvc::domainIntegrate(psi*psiAdj));

psi.correctBoundaryConditions();
psiAdj.correctBoundaryConditions();

volScalarField beta
(
    IOobject
    (
        "beta",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    psi*psiAdj
);

dimensionedScalar epsV = fvc::domainIntegrate(beta);

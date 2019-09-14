//- Rescale psi and psiAdj to have a normalised beta
psiAdj *=  (vol / fvc::domainIntegrate(psi));

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
        IOobject::NO_WRITE
    ),
    psi*psiAdj
);

dimensionedScalar smallPsi("sp",dimless,small);
//- Calculate modified velocity Uplus
volVectorField Ustar
(
    IOobject
    (
        "Uplus",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (
        beta*(U+Up) + (D*psi*fvc::grad(psiAdj))
        - (D*psiAdj*fvc::grad(psi))
    )
);

//- Calculate modified flux
const surfaceScalarField phistar("phistar",fvc::flux(Ustar));

//- Calculate effective macroscopic velocity
const dimensionedVector Ueff ("Ueff", fvc::domainIntegrate(Ustar) / epsV);

//- Calculate source term for X equation
const volVectorField sourceX
(
    IOobject
    (
        "sourceX",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    beta*Ueff-Ustar + fvc::grad(beta*D)
);

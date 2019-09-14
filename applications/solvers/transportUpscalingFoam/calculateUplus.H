dimensionedScalar smallPsi("sp",dimless,small);
//- Calculate modified velocity Uplus
volVectorField Uplus
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
        U + (D*fvc::grad(psiAdj)/(psiAdj+smallPsi))
        - (D*fvc::grad(psi)/(psi+smallPsi))
    )
);

//- Calculate modified flux
const surfaceScalarField phiplus("phiplus",fvc::flux(beta*Uplus));

//- Calculate effective macroscopic velocity
const dimensionedVector Ueff ("Ueff", fvc::domainIntegrate(Uplus) / vol);

//- Calculate source term for X equation
const volVectorField sourceX
(
    IOobject
    (
        "sourceX",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    beta*(Ueff-Uplus) + fvc::grad(beta*D)
);

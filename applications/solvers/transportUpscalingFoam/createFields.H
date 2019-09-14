Info << "Reading field U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field p" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info << "Reading field psi" << endl;
volScalarField psi
(
    IOobject
    (
        "psi",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Create psi adjoint field" << endl;
volScalarField psiAdj
(
    IOobject
    (
        "psiAdj",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    psi
);

Info << "Reading first order corrector field X" << endl;
volVectorField X
(
    IOobject
    (
        "X",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



//- Eigenvalues
dimensionedScalar lambda("lambda",dimless/dimTime,scalar(1));
dimensionedScalar lambdaAdj("lambdaAdj",dimless/dimTime,scalar(1));
//- Total cell volume
dimensionedScalar vol("vol",dimVolume, gSum(mesh.V()));

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar D
(
    "D",
    dimViscosity,
    transportProperties.lookup("D")
);


#include "createFvOptions.H"

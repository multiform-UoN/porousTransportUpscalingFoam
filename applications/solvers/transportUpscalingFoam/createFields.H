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
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
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

//- Field for Aitken relaxation (stores increments of X)
volVectorField aitkenDelta(X);

volScalarField aitkenDeltaPsi(psi);
volScalarField aitkenDeltaPsiAdj(psiAdj);

//label XRefCell = 0;
//vector XRefValue = vector::zero;

//- Eigenvalues
dimensionedScalar lambda("lambda",dimless/dimTime,scalar(1e-10));
dimensionedScalar lambdaAdj("lambdaAdj",dimless/dimTime,scalar(1e-10));
//- Total cell volume
dimensionedScalar vol("vol",dimVolume, gSum(mesh.V()));

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

surfaceScalarField phiAdj("phiAdj",-phi);

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

volScalarField pot
(
    IOobject
    (
        "pot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimless
);

volVectorField Up
(
    IOobject
    (
        "Upot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("Up",dimVelocity,vector::zero)
);

if(transportProperties.found("inverseDistancePotential"))
{
    dictionary potentialDict
    (
        transportProperties.subDict("inverseDistancePotential")
    );

    Info<< "Reading potential properties\n" << endl;

    Switch activePot(potentialDict.lookup("active"));

    if(activePot)
    {
        dimensionedScalar refDistance
        (
            potentialDict.lookup("potRefD")
        );

        dimensionedScalar minDistance
        (
            potentialDict.lookup("potMinD")
        );

        dimensionedScalar refSpeed
        (
            potentialDict.lookup("potRefU")
        );

        const volScalarField wd(wallDist(mesh).y()/refDistance);

        pot = pow(wd+minDistance/refDistance,-1);

        Up = fvc::grad(pot*refDistance)*refSpeed;
        Up.write();

        //- Simply add to flux since (like the velocity field) it does not change in
        //  time
        phi    += fvc::flux(Up);
        phiAdj -= fvc::flux(Up);

    }
}

if(transportProperties.found("grainPotential"))
{
    dictionary potentialDict
    (
        transportProperties.subDict("grainPotential")
    );

    Info<< "Reading potential properties\n" << endl;

    Switch activePot(potentialDict.lookup("active"));

    if(activePot)
    {

        scalarField grainPots
        (
            potentialDict.lookup("grainPotentials")
        );

        vectorField grainCtrs
        (
            potentialDict.lookup("grainCenters")
        );


        forAll(grainCtrs,grain)
        {
            dimensionedVector  grainC
            (
                "grainCtr",
                dimLength,
                grainCtrs[grain]
            );

            dimensionedScalar  grainPot
            (
                "grainPot",
                dimLength,
                grainPots[grain]
            );

            pot += grainPot*pow(mag(mesh.C()-grainC),-1);
        }

        //- Smooth potential
        surfaceScalarField potS(fvc::interpolate(pot));
        pot = fvc::reconstructMag(potS*mesh.magSf());

        dimensionedScalar dim("dim",dimArea/dimTime,scalar(1));

        Up = -fvc::grad(pot)*dim;
        Up.write();

        //- Simply add to flux since (like the velocity field) it does not change in
        //  time
        phi    += fvc::flux(Up);
        phiAdj -= fvc::flux(Up);

    }
}

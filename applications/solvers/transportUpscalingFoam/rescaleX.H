{

    X -= fvc::domainIntegrate(psi*X)/vol;
    X.correctBoundaryConditions();
}

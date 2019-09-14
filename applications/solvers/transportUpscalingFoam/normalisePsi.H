{
    //- The eigenfunction psi should be normalised over the fluid volume to ensure
    //  that we obtain a proper series expansion for the average scalar field.

    psi *=
    (
        vol
        /
        fvc::domainIntegrate(mag(psi)+small)
    );

    psi.correctBoundaryConditions();
    volScalarField psiFieldRes
    (
        mag(psi
                -(
                    psi.prevIter()*(
                        vol
                        /
                        fvc::domainIntegrate(mag(psi.prevIter())))
                    )
            )/(mag(psi)+small)
    );

    psiRes = gMax(psiFieldRes);
}

{

    psiAdj *=
    (
        vol
        /
        fvc::domainIntegrate(mag(psiAdj)+small)
    );
   

    psiAdj.correctBoundaryConditions();

    volScalarField psiAdjFieldRes
    (
        mag(psiAdj -psiAdj.prevIter())/(mag(psiAdj.prevIter())+small)
    );

    psiAdjRes = gMax(psiAdjFieldRes);
}

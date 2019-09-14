    // Momentum predictor

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      - fvm::laplacian(nu, U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }

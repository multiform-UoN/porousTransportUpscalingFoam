//- Create controls for Power method
scalar lambdaRes(1.);

//- Read controls
scalar eigenTol
(
    readScalar(powerControl.dict().lookup("eigenvalueTol"))
);


//- Create controls for corrector equation

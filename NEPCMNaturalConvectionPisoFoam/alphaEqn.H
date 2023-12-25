{

    fvScalarMatrix alphaEqn
    (
        fvm::ddt(alpha)
      + fvm::div(phi, alpha)
      - fvm::laplacian(DB, alpha)
      //- fvc::laplacian((DT/T), T)
    );

    alphaEqn.relax();
    alphaEqn.solve();
}
//******************mp Equation *********************************************************** 

{
    volScalarField kappaEff("kappaEff", (turbulence->nut()/Prt) * rho * Cp + kappa);
    kappaEff.correctBoundaryConditions();

    fvScalarMatrix TEqn
    (
	fvm::ddt(rho*Cp, T)
      + rho*Cp*fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
    );

    TEqn.relax();

    TEqn.solve();

    rhok = 1.0 - beta*(T - TRef);
}

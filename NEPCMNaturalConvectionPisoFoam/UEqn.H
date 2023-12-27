
    rhoNEPCM = ((1+l)*rhoc*rhos)/(rhos+l*rhoc);
    rho = alpha*rhoNEPCM+(1-alpha)*rhof;
    DB=kB*T/(3.0*3.14*muf*ds);
    DT=0.26*(kappaf/(2.0*kappaf+kappaNEPCM))*(muf/rhof)*alpha;
    mu = (muf/pow((1.0-alpha),2.5));


    // Solve the momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    UEqn.relax();

    if (piso.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(rhok)*(g & mesh.Sf())
                  - fvc::snGrad(p)*mesh.magSf()
                )
            )
        );
    }

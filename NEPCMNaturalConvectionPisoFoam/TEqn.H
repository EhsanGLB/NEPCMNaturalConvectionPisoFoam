for (int j=0; j<7; j++ ) // added
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

    forAll(fNEPCM, i)
    {
        if (T[i] < scalar(Tf.value()-TMr.value()/2) )
        {
            fNEPCM[i] = scalar(0.0);
            NEPCM_molten[i] = scalar(0.0);
        }
        else if ( scalar(Tf.value()-TMr.value()/2) <= T[i] && T[i] <= scalar(Tf.value()+TMr.value()/2) )
        {
            fNEPCM[i] = scalar(1.0);
            NEPCM_molten[i] = scalar(0.5);
        }
        else if (T[i] > scalar(Tf.value()+TMr.value()/2) )
        {
            fNEPCM[i] = scalar(0.0);
            NEPCM_molten[i] = scalar(1.0);
        }
    }

    CpPCM = Cpc+( (3.14/2)*(hc/TMr-Cpc)*sin( (3.14/TMr)*(T-(Tf-TMr/2)) ) )*fNEPCM;
    CpNEPCM = ((CpPCM+l*Cps)*rhoc*rhos)/((rhos+l*rhoc)*rhoNEPCM);
    betaNEPCM = betac + ( (betas-betac)/2 )*(1 - l*rhos/rhoc);
    Cp = (alpha*rhoNEPCM*CpNEPCM+(1-alpha)*rhof*Cpf)/rho;
    kappa = ((kappaNEPCM+(lambda-1)*kappaf-(lambda-1)*(kappaf-kappaNEPCM)*alpha)/(kappaNEPCM+(lambda-1)*kappaf+(kappaf-kappaNEPCM)*alpha))*kappaf;
    beta = (alpha*rhoNEPCM*betaNEPCM+(1-alpha)*rhof*betaf)/rho;
}

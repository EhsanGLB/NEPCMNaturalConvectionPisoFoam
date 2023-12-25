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

rhoNEPCM = ((1+l)*rhoc*rhos)/(rhos+l*rhoc);
CpPCM = Cpc+( (3.14/2)*(hc/TMr-Cpc)*sin( (3.14/TMr)*(T-(Tf-TMr/2)) ) )*fNEPCM;
CpNEPCM = ((CpPCM+l*Cps)*rhoc*rhos)/((rhos+l*rhoc)*rhoNEPCM);
betaNEPCM = betac + ( (betas-betac)/2 )*(1 - l*rhos/rhoc);
DB=kB*T/(3.0*3.14*muf*ds);
DT=0.26*(kappaf/(2.0*kappaf+kappaNEPCM))*(muf/rhof)*alpha;
rho = alpha*rhoNEPCM+(1-alpha)*rhof;
Cp = (alpha*rhoNEPCM*CpNEPCM+(1-alpha)*rhof*Cpf)/rho;
kappa = ((kappaNEPCM+(lambda-1)*kappaf-(lambda-1)*(kappaf-kappaNEPCM)*alpha)/(kappaNEPCM+(lambda-1)*kappaf+(kappaf-kappaNEPCM)*alpha))*kappaf;
mu = (muf/pow((1.0-alpha),2.5));
beta = (alpha*rhoNEPCM*betaNEPCM+(1-alpha)*rhof*betaf)/rho;


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

    Info<< "Reading thermophysical properties\n" << endl;

    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

///////////////////////////////////////////////////////////////// NEPCM

Info<< "Reading field alpha\n" << endl;
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


//////////////- base fluid
Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh, 
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    ) 
);

dimensionedScalar rhof(transportProperties.lookup("rhof"));
dimensionedScalar Cpf(transportProperties.lookup("Cpf"));
dimensionedScalar kappaf(transportProperties.lookup("kappaf"));
dimensionedScalar muf(transportProperties.lookup("muf"));
dimensionedScalar betaf(transportProperties.lookup("betaf"));
//////////////-core
dimensionedScalar rhoc(transportProperties.lookup("rhoc"));
dimensionedScalar Cpc(transportProperties.lookup("Cpc"));
dimensionedScalar kappac(transportProperties.lookup("kappac"));
dimensionedScalar Tf(transportProperties.lookup("Tf"));
dimensionedScalar TMr(transportProperties.lookup("TMr"));
dimensionedScalar hc(transportProperties.lookup("hc"));
dimensionedScalar betac(transportProperties.lookup("betac"));
//////////////-shell
dimensionedScalar rhos(transportProperties.lookup("rhos"));
dimensionedScalar Cps(transportProperties.lookup("Cps"));
dimensionedScalar kappas(transportProperties.lookup("kappas"));
dimensionedScalar ds(transportProperties.lookup("ds"));
dimensionedScalar betas(transportProperties.lookup("betas"));
//////////////-NPCM
dimensionedScalar l(transportProperties.lookup("l"));
dimensionedScalar lambda(transportProperties.lookup("lambda"));
word alphaEqActivation(transportProperties.lookup("alphaEqActivation"));
dimensionedScalar Prt(transportProperties.lookup("Prt"));
dimensionedScalar TRef(transportProperties.lookup("TRef"));
dimensionedScalar kB(transportProperties.lookup("kB"));
////////////////////

Info<< "Reading field dc\n" << endl;
dimensionedScalar dc = ds*pow( rhos/(rhos+l*rhoc) , scalar(0.3333) );

Info<< "Reading field kappaNEPCM\n" << endl;
dimensionedScalar kappaNEPCM = ds/( (dc/kappac) + ((ds-dc)/kappas) );

Info<< "Reading field rhoNEPCM\n" << endl;
volScalarField rhoNEPCM
(
    IOobject
    (
        "rhoNEPCM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoNEPCM", dimDensity, SMALL)
);

Info<< "Reading field fNEPCM\n" << endl;
volScalarField fNEPCM
(
    IOobject
    (
        "fNEPCM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("fNEPCM", dimless, SMALL)
);

Info<< "Reading field NEPCM_molten\n" << endl;
volScalarField NEPCM_molten
(
    IOobject
    (
        "NEPCM_molten",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("NEPCM_molten", dimless, SMALL)
);

Info<< "Reading field CpPCM\n" << endl;
volScalarField CpPCM
(
    IOobject
    (
        "CpPCM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("CpPCM", dimSpecificHeatCapacity, SMALL)
);

Info<< "Reading field CpNEPCM\n" << endl;
volScalarField CpNEPCM
(
    IOobject
    (
        "CpNEPCM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("CpNEPCM", dimSpecificHeatCapacity, SMALL)
);

Info<< "Reading field betaNEPCM\n" << endl;
volScalarField betaNEPCM
(
    IOobject
    (
        "betaNEPCM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("betaNEPCM", dimless/dimTemperature, SMALL)
);

Info<< "Reading field DB\n" << endl;
volScalarField DB
(
    IOobject
    (
        "DB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DB", dimArea/dimTime, SMALL)
);

Info<< "Reading field DT\n" << endl;
volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DT", dimArea/dimTime, SMALL)
);

Info<< "Reading field beta\n" << endl;
volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("beta", dimless/dimTemperature, SMALL)
);

Info<< "Reading field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rho", dimDensity, SMALL)
);

Info<< "Reading field Cp\n" << endl;
volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Cp", dimSpecificHeatCapacity, SMALL)
);

Info<< "Reading field kappa\n" << endl;
volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kappa", dimThermalConductivity, SMALL)
);

Info<< "Reading field mu\n" << endl;
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("mu", dimMass/dimLength/dimTime, SMALL)
);

///////////////////////////////
    #include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, piso.dict(), pRefCell, pRefValue);
    mesh.schemesDict().setFluxRequired(p.name());


    // Kinematic density for buoyancy force
    volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh
        ),
        1.0 - beta*(T - TRef)
    );

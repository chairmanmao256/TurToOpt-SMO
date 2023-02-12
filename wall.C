/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DASimpleFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DASimpleFoam, 0);
addToRunTimeSelectionTable(DASolver, DASimpleFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DASimpleFoam::DASimpleFoam(
    char* argsAll,
    PyObject* pyOptions)
    : DASolver(argsAll, pyOptions),
      simplePtr_(nullptr),
      pPtr_(nullptr),
      UPtr_(nullptr),
      phiPtr_(nullptr),
      alphaPorosityPtr_(nullptr),
      alphaMaxPtr_(nullptr),
      q1Ptr_(nullptr),
      q2Ptr_(nullptr),
      omegabPtr_(nullptr),
      L_(nullptr),
      DWall_(nullptr),
      S1_(nullptr),
      S2_(nullptr),
      DB_(nullptr),
      DT_(nullptr),
      mgf_(nullptr),
      laminarTransportPtr_(nullptr),
      turbulencePtr_(nullptr),
      daTurbulenceModelPtr_(nullptr),
      daFvSourcePtr_(nullptr),
      fvSourcePtr_(nullptr),
      MRFPtr_(nullptr)
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DASimpleFoam::initSolver()
{
    /*
    Description:
        Initialize variables for DASolver
    */

    Info << "Initializing fields for DASimpleFoam" << endl;
    Time& runTime = runTimePtr_();
    fvMesh& mesh = meshPtr_();
#include "createSimpleControlPython.H"
#include "createFieldsSimple.H"
#include "createAdjointIncompressible.H"
    // initialize checkMesh
    daCheckMeshPtr_.reset(new DACheckMesh(daOptionPtr_(), runTime, mesh));

    daLinearEqnPtr_.reset(new DALinearEqn(mesh, daOptionPtr_()));
    
    this->setDAObjFuncList();

    // initialize fvSource and compute the source term
    const dictionary& allOptions = daOptionPtr_->getAllOptions();
    if (allOptions.subDict("fvSource").toc().size() != 0)
    {
        hasFvSource_ = 1;
        Info << "Initializing DASource" << endl;
        word sourceName = allOptions.subDict("fvSource").toc()[0];
        word fvSourceType = allOptions.subDict("fvSource").subDict(sourceName).getWord("type");
        daFvSourcePtr_.reset(DAFvSource::New(
            fvSourceType, mesh, daOptionPtr_(), daModelPtr_(), daIndexPtr_()));
    }
}

label DASimpleFoam::solvePrimal(
    const Vec xvVec,
    Vec wVec)
{
    /*
    Description:
        Call the primal solver to get converged state variables

    Input:
        xvVec: a vector that contains all volume mesh coordinates

    Output:
        wVec: state variable vector
    */

#include "createRefsSimple.H"

    // change the run status
    daOptionPtr_->setOption<word>("runStatus", "solvePrimal");

    // first check if we need to change the boundary conditions based on
    // the primalBC dict in DAOption. NOTE: this will overwrite whatever
    // boundary conditions defined in the "0" folder
    dictionary bcDict = daOptionPtr_->getAllOptions().subDict("primalBC");
    if (bcDict.toc().size() != 0)
    {
        Info << "Setting up primal boundary conditions based on pyOptions: " << endl;
        daFieldPtr_->setPrimalBoundaryConditions();
    }

    // call correctNut, this is equivalent to turbulence->validate();
    daTurbulenceModelPtr_->updateIntermediateVariables();

    Info << "\nStarting time loop 1: wall distance calculation\n"
         << endl;

    // deform the mesh based on the xvVec
    this->pointVec2OFMesh(xvVec);

    // check mesh quality
    label meshOK = this->checkMesh();

    if (!meshOK)
    {
        return 1;
    }

    primalMinRes_ = 1e10;
    label printInterval = daOptionPtr_->getOption<label>("printInterval");
    label printToScreen = 0;

    //solve the wall distance
    while(this->loop(runTime) && runTime.timeIndex() < 500)
    {
        printToScreen = this->isPrintTime(runTime, printInterval);

        if (printToScreen)
        {
            Info << "Time = " << runTime.timeName() << nl << endl;
        }

        // fvScalarMatrix LEqn(
        //     fvm::ddt(L) - fvm::laplacian(DT * mag(DB * (fvc::grad(L))), L) - fvm::laplacian(1e-4 * DT, L) + fvm::Sp(mgf * alphaPorosity, L) - fvc::Su(S1, S2)
        // );

        fvScalarMatrix LEqn(
            fvm::ddt(L) - fvm::laplacian(DT, L) + fvm::Sp(mgf * alphaPorosity, L) - fvc::Su(S1, S2)
        );

        // get the solver performance info such as initial
        // and final residuals
        SolverPerformance<scalar> solverL = LEqn.solve();
        //this->primalResidualControl<scalar>(solverL, printToScreen, printInterval, "L");

        if (printToScreen)
        {
            //this->calcPrimalResidualStatistics("print");
            

            Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                 << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                 << nl << endl;
        }

        DWall = (-mag(fvc::grad(L)) * DB + sqrt(2.0 * mag(L) + DB * DB * magSqr(fvc::grad(L)))) * DB;
        // DWall = (-magSqr(fvc::grad(L) * DB) + pow(mag(1.5 * L + pow(mag(DB * fvc::grad(L)),3.0)), 2.0/3.0)) * DB;

        //runTime.write();
    }

    Info << "\nStarting time loop 2: flow field calculation\n"
         << endl;

    // set the rotating wall velocity after the mesh is updated (if MRF is active)
    this->setRotingWallVelocity();

    // if the forwardModeAD is active, we need to set the seed here
#include "setForwardADSeeds.H"

    primalMinRes_ = 1e10;
    printToScreen = 0;

    const volScalarField& nut_ = mesh.thisDb().lookupObject<volScalarField>("nut");
    dimensionedScalar nuLam
    (
    "nuLam",
    dimensionSet(0,2,-1,0,0,0,0),
    scalar(1.5e-5)
    );
    
    while (this->loop(runTime)) // using simple.loop() will have seg fault in parallel
    {

        printToScreen = this->isPrintTime(runTime, printInterval);

        if (printToScreen)
        {
            Info << "Time = " << runTime.timeName() << nl << endl;
        }

        p.storePrevIter();

        // --- Pressure-velocity SIMPLE corrector
        {
#include "UEqnSimple.H"
#include "pEqnSimple.H"
        }

        laminarTransport.correct();
        daTurbulenceModelPtr_->correct();

        if (printToScreen)
        {
            daTurbulenceModelPtr_->printYPlus();

            this->printAllObjFuncs();

            Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                 << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                 << nl << endl;
        }

        volTensorField gradU_ = fvc::grad(U);

        dissipation = 0.5 *(nut_ + daTurbulenceModelPtr_->nu()) * ((gradU_ + gradU_.T()) && (gradU_ + gradU_.T()));

        magCurl = magSqr(fvc::curl(U));
        runTime.write();
    }

    this->writeAssociatedFields();

    this->calcPrimalResidualStatistics("print");

    // primal converged, assign the OpenFoam fields to the state vec wVec
    this->ofField2StateVec(wVec);

    // write the mesh to files
    mesh.write();

    Info << "End\n"
         << endl;

    return this->checkResidualTol();
}

} // End namespace Foam

// ************************************************************************* //

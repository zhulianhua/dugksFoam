#include "fieldMPIreducer.H"

fieldMPIreducer::fieldMPIreducer (
    Foam::argList& args,
    int* argc,
    char*** argv
)
{
    if(args.optionFound("dvParallel"))
    {
        dvParallel_ = true;
        MPI_Init(argc, argv);
        MPI_Type_contiguous(3, MPI_DOUBLE, &vecType_);
        MPI_Type_commit(&vecType_);
        MPI_Op_create((MPI_User_function *)vectorSum, 1, &opSumVec_);

        MPI_Comm_rank(MPI_COMM_WORLD,&rank_ );
        MPI_Comm_size(MPI_COMM_WORLD,&nproc_);
    }
    else
    {
        dvParallel_ = false;
        rank_ = 0;
        nproc_ = 1;
    }
}

fieldMPIreducer::~fieldMPIreducer()
{
    if(dvParallel_)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    }
}

void fieldMPIreducer::reduceField(volScalarField& vsf)
{
    List<scalar> vsfPart(vsf);
    MPI_Allreduce
    (
        vsfPart.data(),
        vsf.data(),
        vsf.size(),
        MPI_DOUBLE,
        MPI_SUM,
        MPI_COMM_WORLD
    );
    //forAll(vsf.boundaryField(), patchi)
        //reduceField(vsf.boundaryField()[patchi]);
}

void fieldMPIreducer::reduceField(surfaceScalarField& ssf)
{
    List<scalar> ssfPart(ssf);
    MPI_Allreduce
    (
        ssfPart.data(),
        ssf.data(),
        ssf.size(),
        MPI_DOUBLE,
        MPI_SUM,
        MPI_COMM_WORLD
    );
    forAll(ssf.boundaryField(), patchi)
        reduceField(ssf.boundaryField()[patchi]);
}

void fieldMPIreducer::reduceField(volVectorField& vvf)
{
    List<vector> vvfPart(vvf);
    MPI_Allreduce
    (
        vvfPart.data(),
        vvf.data(),
        vvf.size(),
        vecType_,
        opSumVec_,
        MPI_COMM_WORLD
    );
    forAll(vvf.boundaryField(), patchi)
        reduceField(vvf.boundaryField()[patchi]);
}

void fieldMPIreducer::reduceField(surfaceVectorField& svf)
{
    List<vector> svfPart(svf);
    MPI_Allreduce
    (
        svfPart.data(),
        svf.data(),
        svf.size(),
        vecType_,
        opSumVec_,
        MPI_COMM_WORLD
    );
    forAll(svf.boundaryField(), patchi)
        reduceField(svf.boundaryField()[patchi]);
}

void fieldMPIreducer::reduceField(scalarField& sf)
{
    List<scalar> sfPart(sf);
    MPI_Allreduce
    (
        sfPart.data(),
        sf.data(),
        sf.size(),
        MPI_DOUBLE,
        MPI_SUM,
        MPI_COMM_WORLD
    );
}

void fieldMPIreducer::reduceField(vectorField& vf)
{
    List<vector> vfPart(vf);
    MPI_Allreduce
    (
        vfPart.data(),
        vf.data(),
        vf.size(),
        vecType_,
        opSumVec_,
        MPI_COMM_WORLD
    );
}

void fieldMPIreducer::vectorSum( vector *in, vector *inout, int *len, MPI_Datatype *dptr )
{
    int i;
    for (i=0; i< *len; ++i) {
        *inout  = *in + *inout;
        in++; inout++;
    }
}

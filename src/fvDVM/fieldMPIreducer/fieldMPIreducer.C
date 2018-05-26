#include "foam_defs.h"
#include "fieldMPIreducer.H"

#if FOAM_MAJOR <= 3
    #define BOUNDARY_FIELD_REF boundaryField()
#else
    #define BOUNDARY_FIELD_REF boundaryFieldRef()
#endif


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
        MPI_Type_contiguous(9, MPI_DOUBLE, &tensorType_);
        MPI_Type_commit(&vecType_);
        MPI_Type_commit(&tensorType_);
        MPI_Op_create((MPI_User_function *)vectorSum, 1, &opSumVec_);
        MPI_Op_create((MPI_User_function *)tensorSum, 1, &opSumTensor_);

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
        reduceField(ssf.BOUNDARY_FIELD_REF[patchi]);
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
        reduceField(vvf.BOUNDARY_FIELD_REF[patchi]);
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
        reduceField(svf.BOUNDARY_FIELD_REF[patchi]);
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

void fieldMPIreducer::reduceField(tensorField& vf)
{
    List<tensor> vfPart(vf);
    MPI_Allreduce
    (
        vfPart.data(),
        vf.data(),
        vf.size(),
        tensorType_,
        opSumTensor_,
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

void fieldMPIreducer::tensorSum( tensor *in, tensor *inout, int *len, MPI_Datatype *dptr )
{
    int i;
    for (i=0; i< *len; ++i) {
        *inout  = *in + *inout;
        in++; inout++;
    }
}

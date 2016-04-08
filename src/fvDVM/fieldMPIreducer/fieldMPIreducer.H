#ifndef FIELDMPIREDUCER_H
#define FIELDMPIREDUCER_H
#define fvDVM_H
#include <mpi.h>
#include "fvCFD.H"


class fieldMPIreducer
{
private:
    //-  number of MPI process
    label nproc_;
    //-  own MPI rank
    label rank_;
    //-  user defined MPI_Operator
    MPI_Op opSumVec_;
    //-  user defined MPI_Datatype for OpenFOAM vector data
    MPI_Datatype vecType_;

    //- dvParallel ?
    bool dvParallel_;

public:
    //- Constructors
        fieldMPIreducer(Foam::argList& args, int* argc, char*** argv);
    //- Destructor
        ~fieldMPIreducer();

    //- Memper functions

        //- MPI trnasfer
            void reduceField(volScalarField& vsf);
            void reduceField(volVectorField& vvf);
            void reduceField(surfaceVectorField& svf);
            void reduceField(surfaceScalarField& ssf);
            void reduceField(Field<scalar>& sf);
            void reduceField(Field<vector>& vf);

        //- Access
            inline int nproc() const {return nproc_;};
            inline int rank() const {return rank_;};
        //- Static member function
            static void vectorSum( vector *in, vector *inout, int *len, MPI_Datatype *dptr );

};
#endif 

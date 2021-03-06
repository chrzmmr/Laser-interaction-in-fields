#include "diagonalization.h"

#include  <iostream>                       // to include cout, cin

// diagonalized the Hamiltionian for B field and v velocity (orthogonal to B) and give the eigenvectors and eigenvalues and  update all Level[n].Energy_cm
void Diagonalization_Energy(vector <Internal_state> &Level, double B, double v, SelfAdjointEigenSolver<MatrixXd> &es)
{


    // Manifold 2M Sym #=1..21
    // population v=n 2J=2Spin 2N=2L 2Ω Ecm Δ C


    int nb_levels=21; // Level.size()

//    double E0_cm[nb_levels][nb_levels] =
//    {
//        {0.,0.,0.,0},
//        {0.,40000.,0.,0},
//        {0.,0.,40000.,0},
//        {0.,0.,0.,40000.}
//    };
//    double Zeeman_cm_B[nb_levels][nb_levels] =
//    {
//        {0,0,0,0},
//        {0,0,0,0},
//        {0,0,0,0},
//        {0,0,0,0}
//    };
//    double Stark_cm_Bv[nb_levels][nb_levels] =
//    {
//        {0,0,0,0},
//        {0,0,0,0},
//        {0,0,0,0},
//        {0,0,0,0}
//    };



    /******* ORDER OF LEVELS ****************

    DEAD LEVEL M=0
    n=1   1S00 3S1-1	3S10	3S11
    n=2S  1S00	3S1-1	3S10	3S11
    n=2P  1P1-1	1P10	1P11	3P00	3P1-1	3P10	3P11	3P2-2	3P2-1	3P20	3P21	3P22

    *************************************/


    double E0_cm[nb_levels][nb_levels] =
    {
        {-10000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,-6.817589266,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,41147.81261,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,41148.66481,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,41148.66481,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,41148.66481,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,41148.29958,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,41148.29958,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,41148.29958,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,41148.0561,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,41148.23871,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,41148.23871,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,41148.23871,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,41148.3848,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,41148.3848,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,41148.3848,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,41148.3848,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,41148.3848}
    };

    double Zeeman_cm_B[nb_levels][nb_levels] =
   {
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,-0.9337307964640151,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,-0.9337307964640151,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,-0.9337307964640151,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,-0.9337307964640151,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,-0.6602473779824211,0,0,0,-0.6602473779824211,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0.5390897266891428,0,0,0,0,0,-0.7623880028197907,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.6602473779824211,0,0,0,-0.6602473779824211,0},
        {0,0,0,0,0,0,0,0,0,0,0.5390897266891428,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,-0.6602473779824211,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0.6602473779824211,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,-0.6602473779824211,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,-0.7623880028197907,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,-0.6602473779824211,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    };

    double Stark_cm_Bv[nb_levels][nb_levels] =
    {
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,1.8108016266147785e-6,0,-1.8108016266147785e-6,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,-1.0454668065750543e-6,0,1.2804301095629403e-6,0,1.8108016266147785e-6,0,-7.392566684346655e-7,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,1.2804301095629403e-6,0,1.2804301095629403e-6,0,1.2804301095629403e-6,0,-1.2804301095629403e-6,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,1.0454668065750543e-6,0,1.2804301095629403e-6,0,0,0,7.392566684346655e-7,0,-1.8108016266147785e-6},
        {0,0,0,0,0,1.8108016266147785e-6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,-1.8108016266147785e-6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,-1.0454668065750543e-6,0,1.0454668065750543e-6,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,1.2804301095629403e-6,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,1.2804301095629403e-6,0,1.2804301095629403e-6,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,1.2804301095629403e-6,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,1.8108016266147785e-6,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,1.2804301095629403e-6,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,-7.392566684346655e-7,0,7.392566684346655e-7,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,-1.2804301095629403e-6,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,-1.8108016266147785e-6,0,0,0,0,0,0,0,0,0,0,0,0}
    };

    // For small sizes, especially for sizes smaller than (roughly) 16, using fixed sizes is hugely beneficial to performance, as it allows Eigen to avoid dynamic
// So here I use dynamical size
// But I use double and not flat to increase the precision but this is double size so slower !


    MatrixXd H(nb_levels,nb_levels); // Hamiltonian Matrix

    // data Allocation
    for (int i=0; i<nb_levels; i++)
        for (int j=0; j<nb_levels; j++)
        {
            // Hamiltonian
            H(i,j)=  E0_cm[i][j] + B*Zeeman_cm_B[i][j] + B*v*Stark_cm_Bv[i][j]; // H_ij = 0_<i | H | j>_0
        }


    es.compute(H);

    for( int n = 0; n < nb_levels; n++ )
    {

        Level[n].Energy_cm = es.eigenvalues()(n);
    }


}


// diagonalized the Hamiltionian for B field and v velocity (orthogonal to B) and give the eigenvectors and eigenvalues and dipoles update all Level[n].Energy_cm
void Diagonalization_Energy_dipole(vector <Internal_state> &Level, double B, double v,  SelfAdjointEigenSolver<MatrixXd> &es,  MatrixXd d[])
{

//    int nb_levels=4; //   Level.size(); 21 in our case
//
//    double dipole[3][nb_levels][nb_levels]=
//    {
//        {
//            {0,5.4772255750516612,0.0}, // Spol=10 --> sqrt(3ù10)
//            {5.4772255750516612,0,0,0},
//            {0,0,0,0},
//            {0,0,0,0}
//        },
//        {
//            {0,0,5.4772255750516612,0},
//            {0,0,0,0},
//            {5.4772255750516612,0,0,0},
//            {0,0,0,0}
//        },
//        {
//            {0,0,0.,5.4772255750516612},
//            {0,0,0,0},
//            {0,0,0,0},
//            {5.4772255750516612,0,0,0}
//        }
//    };


    int nb_levels=21; //   Level.size(); 21 in our case

    /******* ORDER OF LEVELS ****************

    DEAD LEVEL M=0
    n=1   1S00 3S1-1	3S10	3S11
    n=2S  1S00	3S1-1	3S10	3S11
    n=2P  1P1-1	1P10	1P11	3P00	3P1-1	3P10	3P11	3P2-2	3P2-1	3P20	3P21	3P22

    *************************************/
// dipole sigma- en absorption
// pi
//sigma+
    double dipole[3][nb_levels][nb_levels]=
    {
        {
            {0,0,0,0,4.722,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,3.78688,0,0,0,0,0,0,0,0,0},
            {4.722,0,0,0,0,0,0,0,0,0,0,0,2.18635,0,2.67773,0,0,0,1.54599,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.67773,0,0,0,2.67773,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.78688},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,2.18635,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,-2.67773,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,-2.67773,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,3.78688,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,2.67773,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,1.54599,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
        },
        {
            {0,159.878,0,4.722,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {159.878,0,0,0,0,0,0,0,0,0,-3.78688,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,-2.67773,0,0,0,-2.67773,0,0,0},
            {4.722,0,0,0,0,0,0,0,0,0,0,0,2.18635,0,0,0,0,0,-3.09197,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.67773,0,0,0,-2.67773,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,2.18635,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,-2.67773,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,2.67773,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,-2.67773,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,-3.09197,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,-2.67773,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
        },
        {
            {0,0,4.722,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,3.78688,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.78688,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,-2.67773,0,0,0,2.67773,0,0,0},
            {4.722,0,0,0,0,0,0,0,0,0,0,0,2.18635,0,-2.67773,0,0,0,1.54599,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,2.18635,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,2.67773,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,2.67773,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,1.54599,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,2.67773,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,3.78688,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
        },
    };


    MatrixXd d0[3]; // d0 = matrix dipole in zero field

    d0[0] = MatrixXd(nb_levels,nb_levels);
    d0[1] = MatrixXd(nb_levels,nb_levels);
    d0[2] = MatrixXd(nb_levels,nb_levels);


// dipole matrix element
    for (int i=0; i<nb_levels; i++)
        for (int j=0; j<nb_levels; j++)
        {
            for(int n_polar = -1; n_polar <= 1; n_polar++)
            {
                d0[n_polar+1](i,j)= dipole[n_polar+1][i][j]; // d0[q+1]_ij = 0_<i | d^(q) | j>_0
            }
        }


    /*** diagonalization ***/

// diagonalized the Hamiltionian for B field and v velocity and give the eigenvectors and eigenvalues
    Diagonalization_Energy(Level, B, v, es);


    /**** calcul of the new dipoles ***/

    for(int n_polar = -1; n_polar <= 1; n_polar++)
    {
        d[n_polar+1] =  (es.eigenvectors().adjoint())*d0[n_polar+1]*(es.eigenvectors()); // The new dipole are given by d[polar] = evec^dag.d0[polar].evec
    }
}


// A partir de la matrice des dipole qui contient en première ligne les M_in et dernière colonne les M_out
int Create_dipole_Lines_from_Matrices(const char *nom_file)
{
    int nb_levels=21;

    double matrice[nb_levels+1][nb_levels+1];
    double dipole[3][nb_levels][nb_levels];

    string filename = "Data/matrice_dipole.dat";
    // ifstream file(nom_file);
    ifstream file(filename.c_str());
    if ( !file )
    {
        cerr << "Erreur d'ouverture fichier Level"  << nom_file << endl;
        return 0;
    }

    for (int i=0; i<nb_levels+1; i++)
        for (int j=0; j<nb_levels+1; j++)
        {
            matrice[i][j] =0.;
        }

    for (int i=0; i<nb_levels; i++)
        for (int j=0; j<nb_levels; j++)
        {
            for(int n_polar = -1; n_polar <= 1; n_polar++)
            {
                dipole[n_polar+1][i][j]= 0.; // d0[q+1]_ij = 0_<i | d^(q) | j>_0
            }
        }


    while (!file.eof())
    {
        for (int i=0; i<nb_levels+1; i++)
            for (int j=0; j<nb_levels+1; j++)
            {
                file >> matrice[i][j] ;
            }
    }

    for (int i=0; i<nb_levels+1; i++)
        for (int j=0; j<nb_levels+1; j++)
        {
            cout  << " i  =" << i << " j = " << j   << "    "  << matrice[i][j] << endl;
        }


    ofstream dipole_file("Data/dipoles_diagonalization.dat");
    for (int i=1; i<nb_levels+1; i++)
        for (int j=0; j<nb_levels; j++)
        {
            for(int n_polar = -1; n_polar <= 1; n_polar++)
            {
                if(matrice[i][nb_levels] == matrice[0][j] + n_polar) // M_up = m_low -1
                    dipole[n_polar+1][i-1][j]= matrice[i][j]; // d0[q+1]_ij = 0_<i | d^(q) | j>_0
            }
        }

    for(int n_polar = -1; n_polar <= 1; n_polar++)
    {
        dipole_file << "{";
        for (int i=0; i<nb_levels; i++)
        {
            dipole_file << "{";
            for (int j=0; j<nb_levels; j++)
            {

                dipole_file   << dipole[n_polar+1][i][j] << ",";
            }
            dipole_file << "}," << endl;
        }
        dipole_file << "}," << endl;
    }
    dipole_file << ",}" << endl;

    cout  << " dans le fichier final il faudra juste remplacer le ',}' par des '}' " << endl;
    file.close();

    exit(1);
    return 0;
}





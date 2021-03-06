#include "diagonalization.h"

#include  <iostream>                       // to include cout, cin

// diagonalized the Hamiltionian for B field and v velocity (orthogonal to B) and give the eigenvectors and eigenvalues and  update all Level[n].Energy_cm
void Diagonalization_Energy(vector <Internal_state> &Level, double B, double v, SelfAdjointEigenSolver<MatrixXd> &es)
{
    // Manifold 2M Sym #=1..21
    // population v=n 2J=2Spin 2N=2L 2Ω Ecm Δ C

    int nb_levels=21; // Level.size()


    /******* ORDER OF LEVELS (the n=0 and n=1 manifold, the one for spontaneous emission, should be order in Energy) ****************



    The matrix are calculated using a Mathematica code which use n S L J M_J ordering that is the most natural orderning. So we use it here.

            Dead Level     1S00     3S1-1   3S10    3S11,     1S00    3S1-1    3S10    3S11,    1P1-1    1P10    1P11    3P00    3P1-1    3P10    3P11    3P2-2    3P2-1    3P20    3P21    3P22

    number          1        2        3       4        5        6        7       8       9       10       11      12       13      14      15      16       17       18      19    20        21
    in C++ [i]      0        1        2        3       4        5        6       7       8       9        10      11       12      13      14      15       16       17      18    19        20


    Howver the diagonalization then gives an ordering in energy that is

            Dead Level     1S00     3S1-1    3S11    3S10,     1S00    3P00   3P1-1    3P10      3P11    1P1-1    1P10    1P11   3P2-2    3P2-1    3P20    3P21    3P22    3S1-1    3S10    3S11

    correspond:  1            2       3      5        4         6       13      14      15       16      10       11      12      17      18       19      20      21        7       8      9
    in C++ [i]   0            1       2      4        3         5       12      13      14       15      9        10      11      16      17       18      19      20        6       7      8

new number C++[i]0            1       2      3         4        5        6       7       8       9       10       11      12      13      14      15       16       17      18      19      20


    Indeed the block ordered in energy are:

    DEAD LEVEL M=0
    n=1   1S00

    3S1-1    3S11    3S10

    1S00

    3P00

    3P1-1    3P10    3P11

    1P1-1    1P10    1P11

    3P2-2    3P2-1    3P20    3P21    3P22

    3S1-1    3S10    3S11




    BUT in the following the order in only ordered in Energy in the ground state
    DEAD LEVEL M=0
    n=1   1S00 3S1-1    3S11    3S10
    n=2S  1S00	3S1-1	3S10	3S11
    n=2P  1P1-1	1P10	1P11	3P00	3P1-1	3P10	3P11	3P2-2	3P2-1	3P20	3P21	3P22






    *************************************/



// I add a small shift to ensure the energy orderning in the ground state
    double E0_cm[nb_levels][nb_levels] =
    {
        {-10000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,-6.817589266,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,-0.000000000001,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.000000000001,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
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
        {0,0,0,0,-0.9337307964640151,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,-0.9337307964640151,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
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
            H(i,j)=  E0_cm[i][j] + B*Zeeman_cm_B[i][j] + B*v*Stark_cm_Bv[i][j]; // H_ij = 0_<i | H | j>_0 = 0_<j | H | i>_0
//           cout << " i,j " << i <<  " " << j << "   " << E0_cm[i][j] << " " << Zeeman_cm_B[i][j] << "  "  << Stark_cm_Bv[i][j] << endl;
        }

    es.compute(H);  // calculate the new eigenvectors |i> (i start from 0)  and new eingen_Energies
    // E_i = es.eigenvalues()(i) (in incresing order)
    // es.eigenvectors()(i,j) = 0<i | j>  gives (i=line, j = column index) the new (column) vector |j> in function of the old |i>_0

    for( int n = 0; n < nb_levels; n++ )
    {
        Level[n].Energy_cm = es.eigenvalues()(n);
    }
}


// diagonalized the Hamiltionian for B field and v velocity (orthogonal to B) and give the eigenvectors and eigenvalues and dipoles update all Level[n].Energy_cm
void Diagonalization_Energy_dipole(vector <Internal_state> &Level, double B, double v,  SelfAdjointEigenSolver<MatrixXd> &es,  MatrixXd d[])
{
    int nb_levels=21; //   Level.size(); 21 in our case

    /******* ORDER OF LEVELS (ordered in energy for the stable states n=1) *****************************

    DEAD LEVEL M=0
    n=1   1S00 3S1-1    3S11    3S10
    n=2S  1S00	3S1-1	3S10	3S11
    n=2P  1P1-1	1P10	1P11	3P00	3P1-1	3P10	3P11	3P2-2	3P2-1	3P20	3P21	3P22

    ***************************************************************************************************/

// E_laser=E' + E'^dag
// for absorption (that is for E_j> E_i and a transition i --> j) the Rabi frequency comes from d_{ji} = <j|d.E'|i>
// so for a laser with polarization vector E' along e_q the dipole is <j|d.E'|i> =<j|d_(q)|i> E' with m_j = q+ m_i
// so the dipole is <j|d_q|i> that is coded here as d[q+1][j][i]
//
// q=-1 for sigma-
// q=0 for pi
// q= +1 for sigma+.
// So be CAREFUL m_i - m_j = q (polar) ONLY for E_j> E_i
//
// For (spontaneous or) stimulated emission with the same laser, that is j-->i for E_j> E_i,  the transition is <i|d.E'^dag|j> = (<j|d.E'|i>)*
// that is <i|d.E'|j> if the dipoles are REAL (THAT IS OUR CASE).

/***  <i|d_q|j> is coded here as d[q+1][i][j] (real), i = line, j = column ; that is for a i<-->j transition (with E_i> E_j) ****/

    double dipole[3][nb_levels][nb_levels]=
    {
        {
            {0, 0, 0, 4.722, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // this is d[-1][0][j] = <0|d_-1|j>
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.78688, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // this is d[-1][1][j] = <1|d_-1|j>
            {4.722, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.18635, 0, 2.67773, 0, 0, 0, 1.54599, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.78688},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.67773, 0, 0, 0, 2.67773, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 3.78688, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 2.18635, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, -2.67773, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, -2.67773, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 3.78688, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 2.67773, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 1.54599, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
        },
        {
            {0, 159.878, 0, 0, 4.722, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {159.878, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3.78688, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2.67773, 0, 0, 0, -2.67773, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.67773, 0, 0, 0, -2.67773, 0},
            {4.722, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.18635, 0, 0, 0, 0, 0, -3.09197, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, -3.78688, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 2.18635, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, -2.67773, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 2.67773, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0},
            {0, 0, -2.67773, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, -3.09197, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, -2.67773, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
        },
        {
            {0, 0, 4.722, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // this is d[+1][0][j] = <0|d_+1|j>
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 3.78688, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // this is d[+1][1][j] = <1|d_+1|j>
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.78688, 0, 0, 0, 0},
            {4.722, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.18635, 0, -2.67773, 0, 0, 0, 1.54599, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2.67773, 0, 0, 0, 2.67773, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 3.78688, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 2.18635, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 2.67773, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 2.67773, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 1.54599, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 2.67773, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 3.78688, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
        }
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
                d0[n_polar+1](i,j)= dipole[n_polar+1][i][j]; // d0[q+1]_ij = 0_<i | d_q | j>_0   (i = out and j = in)
                // cout << "polar " << n_polar << " i  "  << i  << " j " << j  << " d " <<  d0[n_polar+1](i,j) << endl;
            }
        }


    /*** diagonalization ***/

// diagonalized the Hamiltionian for B field and v velocity and give the eigenvectors and eigenvalues
    Diagonalization_Energy(Level, B, v, es);


    /**** calcul of the new dipoles ***/

    for(int n_polar = -1; n_polar <= 1; n_polar++)
    {
        d[n_polar+1] =  (es.eigenvectors().adjoint())*d0[n_polar+1]*(es.eigenvectors());

        // evec = es.eigenvectors() verifie evec(j0,j) = 0<j0 | j>  gives (j0=line, j = column index) the new (column) vector |j> in function of the old |j0>_0
        //  d0[q+1]_i0 j0 = 0_<i0 | d_q | j0>_0 . So d[q+1]_ij = <i | d_q | j> = Sum i0,j0   <i |i0>0 0<i0| d_q | j0>0 0<j0|j> =  Sum i0,j0  evec^dag (i,i0)  d_q(i0,j0)  evec(j0,j)
        // The new dipole are given by d[polar] = evec^dag.d0[polar].evec =  <i | d_q | j> with evec_j = |j> = sum_|j>_0   0_<j| j>.
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





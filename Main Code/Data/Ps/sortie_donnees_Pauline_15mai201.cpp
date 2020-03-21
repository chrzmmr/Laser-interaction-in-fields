/*
  Fichier de sortie des données moléculaires: transition, probablilitées ...
  */

#include "sortie_donnees.h"


//Sauve l'état du générateur de nombre aléatoire pour ne jamais reprendre le même
void save_random_generator(gsl_rng * r, const char *nom_file)
{
    FILE *fp;
    fp=fopen(nom_file, "wb");    // Fichier contenant les données du générateur de nombre aléatoire
    int ran_gen_result = gsl_rng_fwrite (fp, r); // This function writes the random number state of the random number generator r
    if (ran_gen_result == GSL_EFAILED)
        cout << " Problem writing to the file: Data/random_gen.txt" << endl;
    fclose(fp);
}


/* Création du fichier de sortie des données totales  */
void sortie_fichier(ofstream &  /* file2 */, Molecule /* Mol[] */)
{

    /* A mettre si on veut une numérotation automatique (ainsi que la partie finale)


    ifstream file_in("Data/file_num.txt"); // Ouvre le fichier contenant le futur numéro de fichier
    int num;
    file_in >> num  ;
    file_in.close();

    char snum[256];
    sprintf(snum,"%d",num);
    string ssnum = snum ;
    string nom_fichier = "Data/param_Ryd_" + ssnum + ".txt" ;

    ofstream file_out(nom_fichier.c_str()); // créer le fichier de donné de sortie param_Ryd_N° .txt contenant i, exc[i],  x[i], y[i], z[i], pot[i], pot_ji (plus proche voisin), theta_ji, d_ji

    */


    /* A mettre si on veut une numérotation automatique (ainsi que la partie commentée initiale )
        file_out.close();
        ofstream file_num_out("Data/file_num.txt"); // créer le fichier contenant le futur numéro de fichier
        file_num_out << ++num;
        file_num_out.close();
    */
}



void Sortie_donnee_etat_int_simple(ofstream & file_out, const vector <Molecule> &Mol, const vector <Laser> &my_laser, const double t, FitParams &params)
{
    file_out<< setprecision(8);
    for(int i=0; i<(int) Mol.size(); i++)
    {
        file_out << "   " << t << " " << i << " " << Mol[i].exc  << "   " << Mol[i].two_J/2.  << "   " << Mol[i].two_N/2. << "   " << Mol[i].two_M/2. << "   " << Mol[i].Energy0_cm << endl ;
    }
    // cout << endl << "    time t =  " << t << endl  ;
    file_out << endl;

    return;
}





void Sortie_test_debug(ofstream & file_out,  vector <Molecule> &Mol,  const vector <Internal_state> &Level, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, const int nb_mol,  FitParams &params,  DataCards &data, const int number_photons)
{

//    int number_matrix_lines = 2480;
//    fieldB.Calculate_Derivative_Matrix(number_matrix_lines, "Data/Na/MagneticField2D_derivative.dat");

    double step_r = 0.001;
    double step_z = 0.001;
    double nb_steps = 100;

    // for (double  x= -nb_steps*step_r; x < nb_steps*step_r; x+=step_r)
    //for (double  y= -nb_steps*step_r; y < nb_steps*step_r; y+=step_r)
    for (double  x= 0.003; x < 0.005; x+=step_r)

        for (double  z= -nb_steps*step_z; z < nb_steps*step_z; z+=step_z)
        {
            double y=0.;
            file_out << sqrt(x*x+y*y) << " " <<  x << " " << y << " " << z << " "  << fieldB.get_Field(Vecteur3D(x,y,z)) << " " << fieldB.get_grad_field_F2(Vecteur3D(x,y,z)) << endl;
        }
    file_out << endl;

    return;
}





void Sortie_donnee(ofstream & file_out,  vector <Molecule> &Mol,  const vector <Internal_state> &Level, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, const int nb_mol,  FitParams &params,  DataCards &data, const int number_photons)
{
    //set_pot_all_mol(Mol, fieldB, fieldE, laser, t, nb_mol, params); //Met à jour de tous les potentiels (gravité, PAS dipolaire, magnétique, electrique, ...) avec la nouvelle position pour uen sortie
// ATTENTIION THIS DOES NOT WORK FOR THE POTENTIALS


// SOrtie des paramètres scannés
    if( ((int) params.LocateParam("is_param_scanned_out")->val) == ((int) true) )
    {
        for (ParamIterator i=params.begin(); i != params.end(); ++i) // boucle sur les paramètres
        {
            Param &p = **i;
            if (p.is_scanned == true)
            {
                // file_out << p.name ;
                file_out <<  "  " << p.val_t0 << " " ;
            }
        }
    }

    file_out<< setprecision(8);

//    const int i = (int) params.LocateParam("num_niveau_etudie")->val; // numéro du niveau étudié pour faire des stats. -1 pour toutes les molécules
//
//
    Stat stat_Mol;
//    stat_molecule_un_niveau(Mol, stat_Mol, Level,  fieldB, fieldE, i, Mol.size(), params);


//    for (int i = 0; i < nb_mol; i++)
//    {
//         set_pot_mol(Mol[i], fieldB, fieldE, laser, t, params); //Met à jour de tous les potentiels (gravité, PAS dipolaire, magnétique, electrique, ...) avec la nouvelle position pour une sortie
//    }

//    double Temp_ini_z  = params.LocateParam("Temp_ini_z[0]")->val;
//    const double vitesse_therm_z = sqrt(kB*Temp_ini_z/Mol[0].get_mass()); // vitesse thermique en m/s 1/2 m vx^2 = 1/2 kB T



//    file_out  << t << " ";
//    file_out << number_photons << " ";
//    file_out << stat_Mol.Temp_3D_50 << " ";
//    file_out  << stat_Mol.population << " ";
//    file_out  << stat_Mol.Temp_1D_50.z() << " ";


//  Attention relative à la température E_pot = 3/2 k T, E_cin =3/2 kT; E tot=3 kT

//    file_out  << (stat_Mol.E_pot/kB)/mK/1.5/nb_mol << "  ";
//    file_out  << (stat_Mol.E_cin/kB)/mK/1.5/nb_mol << "  ";
//    file_out << (stat_Mol.E_pot+stat_Mol.E_cin)/kB/mK/1.5/nb_mol << "  ";


   // for (int i = 0; i < nb_mol; i++)
    //{
//        double z_init = 3.*(i-nb_mol/2.)*(params.LocateParam("size_z[0]")->val)/nb_mol;
//        double vz_init = 3.*(i-nb_mol/2.)*vitesse_therm_z/nb_mol;
//        file_out  << z_init << " ";
//        file_out  << Mol[i].get_pos().z() << " ";
//        file_out  << vz_init << " ";
//        file_out  << Mol[i].get_vel().z() << " ";
//        file_out  << Mol[i].get_pos() << " ";
        //    file_out  << Mol[i].get_vel() << " ";
//         file_out  << Euler_angles(Mol[i].get_vel()) << " ";
//        file_out  << (Mol[i].get_vel().z()-vz_init)/t/1000. << " "; // a=v/t et acc 1e3 m.s^-2
//        file_out  << 2. * (Mol[i].get_pos().z()-vz_init*t-z_init)/(t*t)/1000. << " "; // dx = v T + 1/2 a T^2
//        file_out  << Mol[i].two_M  << " ";
//        file_out  << t << " ";
//        file_out << endl;
  //  }



    /**  Stat for specific states v  **/

//    vector <Molecule> liste_Adresses_Mol_dans_niveau; //List of Molecule of first type
//    int nb_Mol_in_this_state = 0;
//    for (int i=0; i!= nb_mol; i++)
//        if (Mol[i].v == 0) // Molecules not in v=1,2,3 of X state
//        {
//            liste_Adresses_Mol_dans_niveau.push_back(Mol[i]);
//            nb_Mol_in_this_state++;
//        }
//    stat_molecule_form_list(Mol, liste_Adresses_Mol_dans_niveau, stat_Mol, Level, fieldB, fieldE,  nb_mol, params);
//    file_out << nb_Mol_in_this_state << " ";
//    file_out  << stat_Mol.Temp_1D_50.z() << " ";
//    file_out  << stat_Mol.population << " ";
//    file_out << stat_Mol.E_pot/kB/mK/1.5/nb_mol << "  ";
//    file_out << stat_Mol.E_cin/kB/mK/1.5/nb_mol << "  ";
//    file_out << (stat_Mol.E_pot+stat_Mol.E_cin)/kB/mK/1.5/nb_mol << "  ";
//    liste_Adresses_Mol_dans_niveau.clear(); // erase the vector:
//
//    cout << " t " << t << " photons = " << (double) number_photons ;
//    cout << " Epot " << stat_Mol.E_pot/kB/mK/1.5/nb_mol << "  ";
//    cout << " Ekin " <<  stat_Mol.E_cin/kB/mK/1.5/nb_mol << "  ";
//    cout << " E " << (stat_Mol.E_pot+stat_Mol.E_cin)/kB/mK/1.5/nb_mol << "  ";
//    cout << endl;

    /**  Stat for specific states of "best" molecules in the sens of position  **/
    vector <Molecule> liste_Adresses_Mol; //List of Molecule of this type
    //  double size_limite = params.LocateParam("size_x[0]")->val;
    double size_limite = 0.02;
    int nb_Mol_in_this_state = 0;
    double niveau_moyen = 0.;

    int num_manifold_not_studied = (int) params.LocateParam("num_manifold_not_studied")->val;



    for (int i=0; i!= nb_mol; i++)
    if(Mol[i].exc != num_manifold_not_studied)
        {
        if (Mol[i].get_pos().mag() < size_limite ) // Molecules  within initial size (in x)
        {
            liste_Adresses_Mol.push_back(Mol[i]);
            nb_Mol_in_this_state++;
            niveau_moyen += Mol[i].exc;
        }
        }
    niveau_moyen = niveau_moyen/nb_Mol_in_this_state;
    stat_molecule_form_list(Mol, liste_Adresses_Mol, stat_Mol, Level, fieldB, fieldE,  nb_mol, params);

    if(t!=0)
    {
    file_out  << t/1e-9 << " ";
    file_out << nb_Mol_in_this_state << " ";
    file_out << niveau_moyen << " ";
    // file_out  << stat_Mol.sigma_pos.mag() << " ";
    // file_out << (stat_Mol.E_pot+stat_Mol.E_cin)/kB/mK/1.5/nb_mol << "  ";
    file_out << stat_Mol.E_cin/kB/mK/1.5/stat_Mol.population << "  ";
    file_out << stat_Mol.E_cin/kB/mK/1.5 << "  ";
    file_out << stat_Mol.E_cin<< "  ";

     // cout << " N " <<  nb_Mol_in_this_state <<  " T " << stat_Mol.E_cin/kB/mK/1.5/nb_mol << endl;





    for (int i = 0; i < nb_mol; i++)
    {
        set_pot_mol(Mol[i], fieldB, fieldE, laser, t, params); //Met à jour de tous les potentiels (gravité, PAS dipolaire, magnétique, electrique, ...) avec la nouvelle position pour uen sortie
    }

     file_out << stat_Mol.Temp_1D_50.z() << " ";

     file_out << stat_Mol.population << " ";

     file_out << nb_mol << " ";
     }

    liste_Adresses_Mol.clear(); // erase the vector:


    file_out << endl;

    return;
}



double Sortie_donnee_hist(ofstream & file_hist,  vector <Molecule> &Mol,  const vector <Internal_state> &Level, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, const int nb_mol,  FitParams &params,  DataCards &data, const int number_photons)
{
    set_pot_all_mol(Mol, fieldB, fieldE, laser, t, nb_mol, params); //Met à jour de tous les potentiels (gravité, PAS dipolaire, magnétique, electrique, ...) avec la nouvelle position pour uen sortie
// ATTENTIION THIS DOES NOT WORK FOR THE POTENTIALS


// SOrtie des paramètres scannés
    if( ((int) params.LocateParam("is_param_scanned_out")->val) == ((int) true) )
    {
        for (ParamIterator i=params.begin(); i != params.end(); ++i) // boucle sur les paramètres
        {
            Param &p = **i;
            if (p.is_scanned == true)
            {
                // file_out << p.name ;
                file_hist <<  "  " << p.val_t0 << " " ;
            }
        }
    }




    file_hist<< setprecision(8);

    Stat stat_Mol;
    stat_molecule_un_niveau(Mol, stat_Mol, Level,  fieldB, fieldE, all, Mol.size(), params);


    for (int i = 0; i < nb_mol; i++)
    {
        set_pot_mol(Mol[i], fieldB, fieldE, laser, t, params); //Met à jour de tous les potentiels (gravité, PAS dipolaire, magnétique, electrique, ...) avec la nouvelle position pour uen sortie
    }


   Vecteur3D r;
    r = Mol[0].get_pos();
    double B = fieldB.get_Field(r).mag();

//Sortie donne pour avoir un histogramme sur les vitesses en km/s
file_hist << t << " " ;
 file_hist  << stat_Mol.population << " ";
 file_hist << endl;
for(int i = 0; i < nb_mol; i++)
{if(Mol[i].exc != 0 && Mol[i].exc != 30)


{
file_hist << Mol[i].get_vel().z()/1000<< " " ;

}}
 file_hist << endl;
double Temp =stat_Mol.Temp_1D_50.z();


    return Temp;
}


// Toutes à la suites en temps
void Sortie_donnee_pop_vJ(ofstream & file_out, const vector <Molecule> &Mol, const int nb_Mol, const double t, const int NX, const int N_two_JX, FitParams &params)
{
    if( ((int) params.LocateParam("is_param_scanned_out")->val) == ((int) true) )
    {
        for (ParamIterator i=params.begin(); i != params.end(); ++i) // boucle sur les paramètres
        {
            Param &p = **i;
            if (p.is_scanned == true)
            {
                // file_out << p.name ;
                file_out <<  "  " << p.val_t0 << " " ;
            }
        }
    }


    int **pX=new int*[NX];
    for (int vX=0; vX<NX; vX++)
    {
        pX[vX]=new int[N_two_JX];

    }

    for (int vX = 0; vX < NX; vX++)
        for (int two_JX = 0; two_JX < N_two_JX; two_JX++)
            pX[vX][two_JX] = 0;

    for (int i = 0; i < nb_Mol; i++) // Calcul des populations dans vX,JX
        if (Mol[i].exc == 0) pX[Mol[i].v][Mol[i].two_J]++;



    // cout << "    time t =  " << t << endl << endl ;
    // file_out << "    time t =  " << t << endl << endl ;

    file_out << t  << " ";

    for (int vX = 0; vX < NX; vX++)
        for (int two_JX = 0; two_JX < N_two_JX; two_JX++)
        {
            file_out << pX[vX][two_JX] << " ";
        }


    file_out <<  endl ;
    for (int vX = 0; vX < NX; vX++)
    {
        delete[] pX[vX];
    }
    delete[] pX;


    return;

}



// Sortie des populations dans l'état vX à la suite les unes des autres en temps
void Sortie_donnee_pop_v(ofstream & file_out, const vector <Molecule> &Mol, const int nb_Mol, const double t, const int NX, FitParams &params, int number_photons)
{
    // SOrtie des paramètres scannés
    if( ((int) params.LocateParam("is_param_scanned_out")->val) == ((int) true) )
    {
        for (ParamIterator i=params.begin(); i != params.end(); ++i) // boucle sur les paramètres
        {
            Param &p = **i;
            if (p.is_scanned == true)
            {
                // file_out << p.name ;
                file_out <<  "  " << p.val_t0 << " " ;
            }
        }
    }

    int *pX=new int[NX];

    for (int vX = 0; vX < NX; vX++)
        pX[vX] = 0;

    for (int i = 0; i < nb_Mol; i++) // Calcul des populations dans vX
        if (Mol[i].exc == 0) pX[Mol[i].v]++;

    cout << "    time t =  " << t << endl << endl ;
    file_out << t   << " ";
    file_out << (double) number_photons/nb_Mol << " ";


    for (int vX = 0; vX < NX; vX++)
    {
        //  cout << "  pop[vX="<< vX << "] = " << pX[vX] << endl;
        file_out <<  pX[vX] << " ";
    }

    file_out <<  endl ;

    delete[] pX;

    return;
}

void Sortie_rate(ofstream & file_rate, const  vector <double> &rate, const vector <type_codage_react> &reaction_list, const vector <Molecule> &Mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const int N_Mol, const double t,  FitParams &params)
{
    file_rate<< setprecision(12);
    file_rate << " time t = " << t << endl;

    int nb_rate = rate.size();
    for (int i = 0; i < nb_rate; i++)
    {

        file_rate  <<  rate[i];
        int n_mol= reaction_list[i].n_mol;
        int n_laser = reaction_list[i].n_laser;
        file_rate << " " << n_mol ;
        file_rate << " " << n_laser;

        Vecteur3D r;
        r = Mol[n_mol].get_pos();
        file_rate  << " " << r ;

        double B = fieldB.get_Field(r).mag();
        double E = fieldE.get_Field(r).mag();
        Internal_state Internal_state_in = Mol[n_mol] ; //  état interne de la molecule
        double Energy_in = Internal_state_in.Energy0_cm + (Internal_state_in.Energy_Shift_B_cm(B) + Internal_state_in.Energy_Shift_E_cm(E));
        Internal_state Internal_state_out = reaction_list[i].final_internal_state ; //  état interne de la molecule après la réaction
        double Energy_out = Internal_state_out.Energy0_cm + (Internal_state_out.Energy_Shift_B_cm(B) + Internal_state_out.Energy_Shift_E_cm(E));
        double Energy_transition_laser_cm = cm/laser[n_laser].get_lambda(); // Energie de la transition laser en cm^-1

//        file_rate <<  " " << B ;
//        file_rate <<  " " << fieldB.get_Field(r).z();
//        file_rate <<  " " << Energy_out- Energy_in - Energy_transition_laser_cm;
//        file_rate <<  " " << Internal_state_out.Energy0_cm;
//        file_rate << " " << Internal_state_in.Energy0_cm ;
        file_rate <<  " " << Energy_out;
        file_rate <<  " " << Energy_in ;
        file_rate <<  " " << (reaction_list[i].final_internal_state).two_M ;
        file_rate <<  " " <<  Mol[reaction_list[i].n_mol].two_M << endl;


    }

    file_rate << endl << endl;

}


void Sortie_laser_spectrum(ofstream & file_out, const vector <Laser> &laser, FitParams &params, int num_laser)
{
    Laser my_laser = laser[num_laser];
    double Energy_transition_laser_cm = cm/my_laser.get_lambda(); // Energie de la transition laser en cm^-1

    file_out<< setprecision(8);
    double intensity0=intensity_Convolution_linewidth(1., 0., 0., my_laser.get_Gamma_Laser(), my_laser.get_type_laser(),num_laser, 0., 0., params);

    for (int E_cm = 0; E_cm < 20000; E_cm++)
    {
        double I_shape = my_laser.transmission_spectrum(E_cm); // Intensité laser façonnée
        double delta =(E_cm - 0.01/my_laser.get_lambda())*Conv_Ecm_delta  ;// detuning de la transition en s^-1
        double I_laser =  intensity_Convolution_linewidth(1., delta, 0., my_laser.get_Gamma_Laser(), my_laser.get_type_laser(),num_laser,Energy_transition_laser_cm, E_cm, params)/intensity0; //
        file_out << E_cm << " " << I_shape << " " << I_laser << endl;
    }

    return;
}

// Debug. Gives state, potential, ...
void Sortie_debug(ofstream & file_rate, const  vector <double> &rate, const vector <type_codage_react> &reaction_list, const vector <Molecule> &Mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const int n_reac, const double t,  FitParams &params)
{
    // t	rate	nlas	z	B	delta	E_fin	E_in	2J_fin	2J_in	2Mfin	2M_in

    file_rate<< setprecision(12);
    //  file_rate << " time t = " ;
    file_rate << t << " ";

    int i=  n_reac;

//   file_rate  <<  " rate[" << i << "] = ";
    file_rate << rate[i] << " ";
    int n_mol= reaction_list[i].n_mol;
    int n_laser = reaction_list[i].n_laser;
//   file_rate << " nMol = " ;
    //  file_rate<< n_mol << " ";
//   file_rate << " nlas = ";
    file_rate << n_laser << " ";

    Vecteur3D r;
    r = Mol[n_mol].get_pos();
    //  file_rate << " pos = " ;
    file_rate << r.z() << " ";

    double B = fieldB.get_Field(r).mag();
    double E = fieldE.get_Field(r).mag();
    Internal_state Internal_state_in = Mol[n_mol] ; //  état interne de la molecule
    double Energy_in = Internal_state_in.Energy0_cm + (Internal_state_in.Energy_Shift_B_cm(B) + Internal_state_in.Energy_Shift_E_cm(E));
    Internal_state Internal_state_out = reaction_list[i].final_internal_state ; //  état interne de la molecule après la réaction
    double Energy_out = Internal_state_out.Energy0_cm + (Internal_state_out.Energy_Shift_B_cm(B) + Internal_state_out.Energy_Shift_E_cm(E));
    double Energy_transition_laser_cm = cm/laser[n_laser].get_lambda(); // Energie de la transition laser en cm^-1


//    file_rate << " B " ;
    file_rate << B << " ";
//   file_rate << " detuning_cm " ;
    file_rate << abs(Energy_out- Energy_in) - Energy_transition_laser_cm << " ";
//   file_rate << " Efin0 " ;
//   file_rate << Internal_state_out.Energy0_cm<< " ";
    //  file_rate << "E_in0 ";
    //  file_rate << Internal_state_in.Energy0_cm << " ";
    //  file_rate << " Efin " ;
    file_rate << Energy_out << " ";
    //  file_rate  << "E_in ";
    file_rate << Energy_in << " " ;
    //  file_rate << " 2J_fin ";
    file_rate << (reaction_list[i].final_internal_state).two_J << " ";
    //  file_rate << " 2J_in ";
    file_rate << Mol[reaction_list[i].n_mol].two_J << " ";
    file_rate << Mol[reaction_list[i].n_mol].two_M << " " ;

    file_rate << endl ;

}


// collisional cross section for charge exchange Ps Pbar.

double Cross_section_Ps(double v,  const int n)
{
    double s1= 1.32e-16;
    double s2= 1.12e-15;
    double v_electron = ALPHA * C/(2. * n);
    double kv = v/v_electron;
    double sigma = n*n*n*n*(s1/(kv * kv) + s2)/(1.+pow(kv/1.8,20)); //  Cf PRA 94, 022714 (2016) fitted

    double B = 4.5*tesla;  // MAGNETIC FIELD (PUT HERE BY HAND !!!!!)

    double E_max_field_ionization = (-QE/(4.*pi*EPSILON0 *(2.*A0)*(2.*A0)))/(16.*n*n*n*n) ; // Maximum field for field iniastion. 1/16n^4 can be replaced by 1/9 n^4  in pure electric field

    if ( v*B < E_max_field_ionization)
    {
        return sigma ;
    }
    else
    {
        return 0. ; // because the particles is field ionized
    }
}




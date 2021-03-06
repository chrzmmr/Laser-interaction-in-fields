#include "laser.h"
#include <complex>




// ------------
// Constructors
// ------------

Laser::Laser()                       // Constructeur par défaut
{
    waist_pos.set(0.,0.,0.);
    direction.set(0.,0.,1.);
    waist.set(1.,1.,0.);
    lambda = 1.;
    Gamma_Laser = 1.;
    Power = 0.;
    polarisation.set(0.,0.,0.);
    polar_angle_degree = 0.;
    type_laser = CW;
    coherent_avec_laser_num = -1;
    spectre_Ecm_attenuation.clear(); // spectre_Ecm_attenuation =  map < double, double > ();
    spectre_Ecm_attenuation.insert ( pair<double,double>(0.,1.) ); // Par défaut le spectre est non façonné. On écrit  spectre_Ecm_attenuation[0.] = 1.;

};


Laser::~Laser()
{}
;                    // Destructeur

Laser::Laser(const Laser & my_laser)               // Constructeur de (re)copie
{
    waist_pos = my_laser.waist_pos;
    direction = my_laser.direction;
    waist = my_laser.waist;
    lambda = my_laser.lambda;
    Gamma_Laser = my_laser.Gamma_Laser;
    Power = my_laser.Power;
    polarisation = my_laser.polarisation;
    polar_angle_degree = my_laser.polar_angle_degree;
    type_laser = my_laser.type_laser;
    coherent_avec_laser_num = my_laser.coherent_avec_laser_num;
    // map < double, double > spectre_Ecm_attenuation; // In order to create properly the object (cf Guarreta course p22)
    spectre_Ecm_attenuation = my_laser.spectre_Ecm_attenuation;
};

Laser & Laser::operator = (const Laser& my_laser)          // Affectation par recopie
{
    if (this != &my_laser) // On vérifie que les objets ne sont pas les mêmes !
    {
        waist_pos = my_laser.waist_pos;
        direction = my_laser.direction;
        waist = my_laser.waist;
        lambda = my_laser.lambda;
        Gamma_Laser = my_laser.Gamma_Laser;
        Power = my_laser.Power;
        polarisation = my_laser.polarisation;
        polar_angle_degree = my_laser.polar_angle_degree;
        type_laser = my_laser.type_laser;
        coherent_avec_laser_num = my_laser.coherent_avec_laser_num;
        spectre_Ecm_attenuation = my_laser.spectre_Ecm_attenuation; // Attention ne recopie pas la map seulement l'adresse.
    }
    return *this;
}



// ----------
// Comparison
// ----------

// On ne test pas le spectre (inutile car cette fonction est inutile en fait)
bool Laser::operator == (const Laser& my_laser) const
{
    return ( waist_pos == my_laser.waist_pos &&
             direction == my_laser.direction &&
             waist == my_laser.waist &&
             lambda == my_laser.lambda &&
             Gamma_Laser == my_laser.Gamma_Laser &&
             Power == my_laser.Power &&
             polarisation == my_laser.polarisation &&
             polar_angle_degree == my_laser.polar_angle_degree &&
             type_laser == my_laser.type_laser &&
             coherent_avec_laser_num == my_laser.coherent_avec_laser_num) ? true : false;
}

bool Laser::operator != (const Laser& my_laser) const
{
    return (waist_pos != my_laser.waist_pos ||
            direction != my_laser.direction ||
            waist != my_laser.waist ||
            lambda != my_laser.lambda ||
            Gamma_Laser != my_laser.Gamma_Laser ||
            Power != my_laser.Power ||
            polarisation != my_laser.polarisation ||
            polar_angle_degree != my_laser.polar_angle_degree ||
            type_laser != my_laser.type_laser ||
            coherent_avec_laser_num != my_laser.coherent_avec_laser_num) ? true : false;
}

//----------------------------------
// Surdéfinition des entrées sorties
//----------------------------------

// Sortie des paramètre sauf le spectre
ostream& operator << (ostream &flux, Laser my_laser)
{
    my_laser.write(flux);
    return(flux);
}

istream& operator >> (istream &flux, Laser & my_laser) //my_laser est modifié!
{
    my_laser.read(flux);
    return(flux);
}


void Laser::read_Spectrum(istream & flux)
{
    double energy_cm, attenuation;
    flux >> energy_cm;
    flux >> attenuation;
    spectre_Ecm_attenuation.insert ( pair<double,double>(energy_cm,attenuation) );  // spectre_Ecm_attenuation[energy_cm] = attenuation; Mais je préfère ainsi car as modifier un std::map dans une boucle for basée sur ses iterators. Même  le fait d'accéder à une clé via l'opérateur [ ] insère cette clé (avec la donnée T()) dans la map.
}

int Laser::read_Spectrum(const char *nom_file)
{
    ifstream file(nom_file);

    if ( !file || nom_file=='\0')
    {
        // cerr << "No able to open the file " << nom_file << endl;  // Better to not put because sometimes their is no file (and so if this line is here, we will have all the time a message) and we just as the defautl values
        file.close();
        return 0; // So spectre_Ecm_attenuation is unchanged and thus compose by the default file
    }

    int i=0;
    spectre_Ecm_attenuation.clear(); // To avoid to insert the default file at the begining

    while (!file.eof())
    {
        this->read_Spectrum(file);
        i++;
    }

    file.close();
    return i;
}


void Laser::write_Spectrum(ostream & flux)
{
    if (spectre_Ecm_attenuation.size() ==0)
        cout << " FICHIER VIDE " << endl;
    map < double, double >::const_iterator itr;

    if (flux == cout)
        for(itr = spectre_Ecm_attenuation.begin(); itr != spectre_Ecm_attenuation.end(); ++itr)
            cout << "Energie(cm^-1) " << (*itr).first << " Value attenuation: " << (*itr).second << endl;
    else
        for(itr = spectre_Ecm_attenuation.begin(); itr != spectre_Ecm_attenuation.end(); ++itr)
            flux << itr->first << itr->second << endl;

}


// ----------------------------------
//   FONCTIONS EN LIGNES
// ----------------------------------

// intensité au waist
double  Laser::intensity()  const
{
    double intensity = 2.*Power/(pi*(waist.x())*(waist.y()));
    // I=2P/(pi*w^2) lorsque w=wX=wY
    // Rappel w.Z=0;

    return intensity;

}


// intensité au waist prenant en compte le spectre
double  Laser::transmission_spectrum(const double energie_trans_cm)  const
{
    map < double, double >::const_iterator itr;

    itr = spectre_Ecm_attenuation.upper_bound (energie_trans_cm);  // itr pointe sur l'élément juste après energie_trans_c
    itr--;

    return (itr->second);
}




// Zone de Rayleigh Vecteur3D
// ZRZ est la zone de Railieh moyenne
Vecteur3D  Laser::Rayleigh_range()  const
{

    double ZRX=pi*waist.x()*waist.x()/lambda;
    double ZRY=pi*waist.y()*waist.y()/lambda; //ZR = pi wo^2/lambda
    double ZRZ=sqrt(ZRX*ZRX+ZRY*ZRY);
    return Vecteur3D(ZRX,ZRY,ZRZ);
}

// waist au point (X,Y,Z)
Vecteur3D  Laser::waist_size(const Vecteur3D& point)  const
{

    double ZRX=pi*waist.x()*waist.x()/lambda;
    double ZRY=pi*waist.y()*waist.y()/lambda; //ZR = pi wo^2/lambda

    double wX = (waist.x())*sqrt(1+(point.z()/ZRX)*(point.z()/ZRX)); // w(z)^2=w0^2(1+(z/zr)^2)
    double wY = (waist.y())*sqrt(1+(point.z()/ZRY)*(point.z()/ZRY)); // w(z)^2=w0^2(1+(z/zr)^2)

    return Vecteur3D(wX,wY,sqrt(ZRX*ZRX+ZRY*ZRY)); // waist selon zone de Rayleigh
}


// intensité au point (X,Y,Z) dans le repère laser centré sur le waist
double  Laser::intensity_repere_sur_waist(const Vecteur3D& point)  const
{
    double wX = (this->waist_size(point)).x();
    double wY = (this->waist_size(point)).y();
    double I0 = (this->intensity())*waist.x()*waist.y()/(wX*wY);


    return I0*exp (-2*point.x()*point.x()/(wX*wX)) * exp (-2*point.y()*point.y()/(wY*wY)); // I=I0 exp(-2r^2/w^2)
    // Plus généralement I= P * exp(-2 (x^2/wx^2))/(sqrt(2pi(wx/2)^2) *  exp(-2 (y^2/wy^2))/(sqrt(2pi(wy/2)^2)
}



// intensité au point (x,y,z)
// I.E. lié au repère du labo
// Pour le calculer on le remet dans le repère du laser
double  Laser::intensity_lab_axis(const Vecteur3D& point)  const
{
    return this->intensity_repere_sur_waist(rotation_lab_axis(point - this->get_waist_pos(),  (Euler_angles(direction)).x(), (Euler_angles(direction)).y(),  (Euler_angles(direction)).z())  );
}


// wave_vector k
// h c / \lambda = h c \sigma(m-1) = \hbar c k
Vecteur3D Laser::wave_vector()  const
{
    double k =  2*pi/lambda  ; //  k = 2pi/lambda
    return k*direction.unit();
}




// Energie de la transition laser en cm^-1
double Laser::Energy_transition_laser_cm()  const
{
    return   0.01/lambda  ;
}

//----------------------------------
// AUTRES fonctions
//--------------------------------

// I = ε0 c E^2 /2.
double champ_E(const double irradiance)
{
    return sqrt(2*irradiance/(C*EPSILON0));
}




// effectif dipole d.e_laser = sum_p d^p epsilon_p
// where the dipole transition vector d= sum_p d^p e_p is given in the local quantification axis
// and the polarisation vector e_laser= sum_p' epsilon_p'^* e_p'  is given in the laser axis
double effectif_dipole_local(const Vecteur3D& dipole, const Vecteur3D& axe_quant,  const Laser& my_laser)
{
    double dp,d0,dm;
    dm = dipole(0);
    d0 = dipole(1);
    dp = dipole(2);

    Vecteur3D Euler_angles_axe_quant = Euler_angles(axe_quant);
    Vecteur3D Euler_angles_axe_laser = Euler_angles(my_laser.get_direction());

// The link between the polar angles (theta, phi) and the Euler angle (alpha, beta, gamma) in ZXZ convention as we used them (Wikipedia)  are
    // alpha= phi +pi/2; beta = theta; gamma = -pi/2
    double phi_F = Euler_angles_axe_quant(0) - pi/2.; // polar angle for the quantization axis (along the field F)
    double theta_F = Euler_angles_axe_quant(1); // polar angle for the quantization axis (along the field F)
    double phi_k = Euler_angles_axe_laser(0) - pi/2.; // polar angle for the laser axis (along the vector k)
    double theta_k = Euler_angles_axe_laser(1); // polar angle for the laser axis (along the vector k)


    // The polarization vector = am exp(i psi) e'_-1 + ap exp(-i psi) e'_-1    is given in the list param by ap (for sigma+), am for sigma- and the polar_angle psi
    // So  epsilon_-1 = am exp(-i psi) and epsilon_+1 = ap exp(i psi)
    Vecteur3D polarization = my_laser.get_polarisation();
    double am = polarization(0);
    double ap = polarization(2);
    double psi = my_laser.get_polar_angle_degree()*pi/180.;
// TODO (dc#4#): We assume real dipole transition in -1,0,+1 polarizaion basis (sigma+, pi, sigma-). ...
//If different be careful



    complex<double> dip_eff;
    const complex<double> i(0., 1.);

    /*** We try to optimize the calcul so to reduce at maximmum the numebr of operations ****/

    double sin_theta_k = sin(theta_k);
    double sin_theta_F = sin(theta_F);
    double cos_theta_k = cos(theta_k);
    double cos_theta_F = cos(theta_F);
    double cos_phi_F_minus_phi_k = cos(phi_F-phi_k);
    double sin_phi_F_minus_phi_k = sin(phi_F-phi_k);
    double sqrt2 = sqrt(2.);

    double d1F = (dm-dp)*cos_theta_F + sqrt2*d0*sin_theta_F;
    double d2F = (dp-dm)*sin_theta_F + sqrt2*d0*cos_theta_F;

    dip_eff  = (am - exp(2.*i*psi)*ap) *
               (
                   cos_theta_k*(cos_phi_F_minus_phi_k*d1F + i*(dm+dp)*sin_phi_F_minus_phi_k ) - sin_theta_k*d2F
               )
               +
               (am + exp(2.*i*psi)*ap)*
               (
                   (dm+dp)*cos_phi_F_minus_phi_k + i*d1F*sin_phi_F_minus_phi_k
               );


    return abs(dip_eff/2.);
}


/***
Compare to the PRA 2014 we change notation now to be like Mathematica or Wikipedia:
 repère x,y,z du labo  et X,Y,Z de l'axe de quantification donné par le champ extérieur local
On utilise les angles d'Euler pour faire les rotations de repère
http://en.wikipedia.org/wiki/Euler_angles qui note (alpha,beta,gamma)


1. the first rotation is by an angle alpha (modulo 2pi, we chooose  (-pi,pi] to use atan2 ) about the z-axis (even if we are going to use acos in 0 pi !)

2. the second rotation is by an angle beta in [0,pi] about the former (new) x-axis (now x')

3. the third rotation is by an angle gamma in (modulo 2pi, we chooose  (-pi,pi] to use atan2 ) about the former z-axis (now z').
***/

// Calcul des angles d'EULER. Pour un repère donné uniquement par son vecteur OZ=direction
// Il reste donc un arbitraire pour  choisir l'angle de rotation autour de cet axe pour son repère
// Nous choissons les angles tel que le repère soit le repère polaire dont OZ est la direction et OX selon le méridien
Vecteur3D  Euler_angles( Vecteur3D direction)
{
    double norm = direction.mag();
    direction = direction/norm;
    double alpha;

    if (direction.x()==0 && direction.y()==0)
        alpha = pi/2.;
    else
        alpha = atan2(direction.x(), - direction.y()); //  alpha  = atan2(Z_1,-Z2); where Z1 is the x coordinate of Z, Z2 is the y, Z3 is the z .  atan2 (y,x)   is defined as the angle  between the positive x-axis and the ray to the point of coordinate (x,y)
    //  The other formula does nt depends on Z1 so is wrong for instance for direction = (-1,0,0) along -Ox :  cos(alpha) = -Z_2 / \sqrt{1 - Z_3^2}. acos( - direction.y() / sqrt(1.00000000000001 -direction.z()*direction.z()) )  to avoid 1-1=0

    double beta = acos(direction.z()); // ArcCos(Z3)
    return Vecteur3D(alpha,beta,-pi/2.);
    // Le lien entre les angles sphériques (polaires) (phi,theta,psi) et les angles d'Euler sont alpha=phi+pi/2; beta=theta; gamma = psi-pi/2
}


/**

On passe du repère X,Y,Z à x,y,z  R=Ar par (transpose of Z1 X2 Z3 convention)
a_(11)	=	cos(gamma)*cos(alpha) - cos(beta)*sin(alpha)*sin(gamma)
a_(12)	=	cos(gamma)*sin(alpha) + cos(beta)*cos(alpha)*sin(gamma)
a_(13)	=	sin(gamma)*sin(beta)
a_(21)	=	-sin(gamma)*cos(alpha) - cos(beta)*sin(alpha)*cos(gamma)
a_(22)	=	-sin(gamma)*sin(alpha) + cos(beta)*cos(alpha)*cos(gamma)
a_(23)	=	cos(gamma)*sin(beta)
a_(31)	=	sin(beta)*sin(alpha)
a_(32)	=	-sin(beta)*cos(alpha)
a_(33)	=	cos(beta)


**/

// Passage des coordonnées point(x,y,z) (labo) à (X,Y,Z): R=Ar donné par les angles d'Euler alpha beta et gamma
// cf http://mathworld.wolfram.com/EulerAngles.html
//Laser coordinates where point=(x,y,z) is lab coordinate
Vecteur3D  rotation_lab_axis(const Vecteur3D& point, double alpha, double beta, double gamma)
{
    Vecteur3D A1,A2,A3; //  A Rotation Matrix

    A1=Vecteur3D(cos(gamma)*cos(alpha)-cos(beta)*sin(alpha)*sin(gamma),cos(gamma)*sin(alpha)+cos(beta)*cos(alpha)*sin(gamma), sin(gamma)*sin(beta));      // A1=(a11,a12,a13) (gamma=0)
    A2=Vecteur3D(-sin(gamma)*cos(alpha)-cos(beta)*sin(alpha)*cos(gamma), -sin(gamma)*sin(alpha)+cos(beta)*cos(alpha)*cos(gamma),cos(gamma)*sin(beta) );
    A3=Vecteur3D(sin(beta)*sin(alpha),  -sin(beta)*cos(alpha), cos(beta));

    double X,Y,Z; //Laser coordinates where point=(x,y,z) is lab coordinate
    X=A1.dot(point);
    Y=A2.dot(point);
    Z=A3.dot(point);

    return Vecteur3D(X,Y,Z);
}


// Passage des coordonnées point (X,Y,Z) à labo (x,y,z) (labo).
// Le repère XYZ est donnée donné par les angles d'Euler alpha beta et gamma par rapport à xyz
Vecteur3D  rotation_axis_lab(const Vecteur3D& point, double alpha, double beta, double gamma)
{
    return rotation_lab_axis(point, -gamma, - beta, -alpha);
}

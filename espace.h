#ifndef espace_h_  ///Definition d’une bibliotheque
#define espace_h_

///Inclusions des bibliotheques
# include <stdio.h>
# include <stdlib.h>
# include <math.h> // pour qu'on puisse faire le calcul des puissance , la valeur abs

///Declaration des constantes
#define permitiviteduvide 0.00000000000885418782         // est une constante physique. Elle est notée ε0 unités F/m
#define taille 51    // la taille de notre espace
#define h 0.1         // un pas constant h de notre intervalle
#define pi 3.14159265  // constante d'Archimède



///Declaration des elements pour calculer le potentiel
typedef struct
{
    double E[2]; // Declaration d'un vecteur champ de taille 2 de type reel
    double v;    // Declaration d'une variable potentiel de type reel
    double epsilon; //Declaration d'une variable epsilon
    double rho;     //Declaration d'une charge rho de type reel ;

}Espace;

/// Prototypes des fonctions  espace
Espace** creerespace();
void voltage (Espace**grille,double v,int*x,int*y,int nbrpoint);
void Electric_fields(Espace**grille);
void charge_density_center (Espace**grille,double rho);
void charge_density_any_point( Espace**grille,int*point,double rho);
void charge_density (Espace**grille,int*x,int*y,int nbrpoints, double densite,int type);
void espace_epsilon(Espace**grille);
void epsilon_r(Espace**grille,int*y,int*x,int n,double epsilonr,int type);
void carre(double *centre,int *cordx,int *cordy,int *dimension);
void cercle(double *centre,int *cordx,int *cordy,int r);
void permitivite_relative(Espace**grille,int*x,int*y,double n,double epsilonr,int type);
#endif // espace


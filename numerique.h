#ifndef numerique_h_            ///Definition díune bibliotheque
#define numerique_h_

#include"espace.h"          //inclusion de la bibliothËque espace//


// Prototypes des fonctions//

void gauss_seidel(int dimension, double**mat,double*b,double*x, int C);
void successive_over_relaxation(int dimension, double**mat,double*b,double*x,double presicion);
void matrices(Espace**grille,double *b,double ** mat,double*x,int type);
double capacite(Espace**grille,int*x,int y,int *cordx2,int cordy2,double V1,double V2);

#endif // MÈTHODES_NUMÈRIQUES_H_INCLUDED

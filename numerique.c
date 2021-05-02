#include"numerique.h"  // Fichier d'en-tête méthode numérique

///*******************************************************************************************************************////
 /** Fonction calculant sert à resoudre l'equation du potentiel électrostatique partout à l'interieur de notre grid   ***
 // et  est une méthode itérative                                                                                     ***
 * @param dimension:la taille*taille de LA matrice (grand matrice mat)                                                          ***
 * @param b: le second membre de l'equation                                                                           ***
 * @param x: le vecteur des solutions                                                                                 ***
 * @param c :nombre d'iterations d'arrets                                                                             ***
 *@ return : la fonction ne retourne rien car elle est du type void                                                   ***
 */

void gauss_seidel(int dimension, double**mat,double*b,double*x, int C)
{
    int compteur=0;//Déclarer une variable compteur initialisé à zéro ,
    double som1=0,som2=0;// Déclarer deux variables qui calcule la somme et qui sont initialisées à zéro ;
    do
    {
        // Déclarer deux variables 'i,j ' avec la valeur initiale est 0;
        //après chaque itération i, j augmente de 1;
        // Et la boucle se terminera lorsque les valeurs sont supérieure à taille ;
        // Remrque: la valeur ++ équivaut à : "i=i+1, j=j+1";
        for(int i=0; i<dimension;i++)
        {
            for(int j=0; j<dimension; j++)
            {
                if(i+j<dimension-1)
                    som1=som1+mat[i][j]*x[j+i+1];//{calcul de la somme des coeff et leur x de l'iteration precedente}
                if(j<i)
                {
                    som2=som2+mat[i][j]*x[j];//{calcul de la sommme des coeff et leur x de l'iteration actuel}
                }

            }
            x[i]=(1./mat[i][i])*(b[i]-som1-som2); //{calcul de la solution}

            som1=0; // on initialise la somme une autre fois  //
            som2=0;
        }
        compteur=compteur+1;// {Augmentez la valeur du compteur d'arret ,après chaque itération };

    }
    while (C>compteur);// {on repète tant que le critère d'arret et supérieur au compteur};

}


///********************************************************************************************************************************
/**Fonction calculant les potentiels X en utilsant cette fois-ci la méthode de relaxation successive(SOR) est une variante de la***
/ la méthode Gauss-Seidel pour résoudre un système linéaire d'équations,ce qui entraine une convergence plus rapide /           ***
 * \param dimension:taille*taille de LA matrice (grand matrice mat)                                                                    ***
 * \param b: le second membre de l'equation                                                                                     ***
 * \param x: le vecteur des solutions                                                                                           ***
 * \param presicion en pourcentage pour indiqué la condition d'arret                                                            ***
 * @return : la fonction ne retoune rien                                                                                        ***
 */
//********************************************************************************************************************************

void successive_over_relaxation(int dimension, double**mat,double*b,double*x,double presicion)
{


    int i,k,f=0;
    double som1=0;
    double som2=0;
    double critere;
    double temp;

    double t=2*cos(pi/taille); //calcul de t pour le coefficient de ponderation
    double w=(8-sqrt(64-16*pow(t,2)))/(pow(t,2));//calcul du coefficient de ponderation

    do
    {

        critere=0;
        for(k=0; k<dimension; k++)
        {


            for(i=0; i<dimension; i++)
            {
                if(k+i<dimension-1)
                    som1=som1+x[k+1+i]*mat[k][k+1+i];//les memes sommes comme gauss seidell
                if(i<k)
                {
                    som2=som2+x[i]*mat[k][i];
                }



            }
            temp=x[k];//la solution precedente
            x[k]=(1-w)*x[k]+(w/mat[k][k])*(b[k]-som1-som2);//calcul de la solution



            if(f>500 && abs(temp-x[k])<presicion)
                
                // on a choisi 500 car c'est la valeur qui nous permet de ne pas avoir une matrice nul au milieu 
                critere=critere+1;

            else if(f<500)
                critere=-1;






            som1=0; // on réinitialise la somme une autre fois //
            som2=0;

        }
        f=f+1;//condition en plus dans le cas ou la presicion est vraiment inferieur
    }
    while ( (critere<taille) && f<3000);

}


///*******************************************************************************************************////
/** La fonction sert à remplir notre matrice avec les coefficients et remplir le vecteur b                  ***
    *  @param  space: notre espace (monde)                                                                  ***
    *  @param  b:  le second membre
    *  @return mat: la matrice qu'on veut resoudre                                                           ***
    *  @param type: choix de la fonction soit avec un diélectrique ou bine sans utilisant un diélectrique    ***
    *  @ return : la fonction ne retourne rien car elle est de type void                                     ***
*/
///******************************************************************************************************//////

void matrices(Espace**grille,double *b,double ** mat,double*x,int type)
{

    int i,j;

    if(type==0)
    {
        for(i=0; i<taille*taille; i++)
        {
            for(j=0; j<taille*taille; j++)
            {
                //les conditions sur i pour qu'on soit  sur  les extrimités ou sur un point avec un potentiel definie par l'utilisateur.
                if((i>=0 && i<taille)|| (i>(taille-1)*taille && i<(taille*taille)) ||((i+1)%(taille)==0 || i%(taille)==0) || (grille[i/taille]+i%taille)->v!=0)

                {
                    //on definit le coeff de diagonale a 1 et le rest a 0
                    if (i==j)
                        mat[i][j]=1;
                    else
                        mat[i][j]=0;

                }
                else
                {
                    //pour les autres point
                    if(j==i+1 || j==i-1 ||j==i+taille ||j==i-(taille) )
                    {
                        mat[i][j]=1;
                    }
                    else if(j==i)
                    {
                        mat[i][j]=-4;
                    }
                    else
                    {
                        mat[i][j]=0;
                    }
                }
            }
        }
    }


    if(type==1)
    {
        //le cas de dielectrique

        for(i=0; i<taille*taille; i++)
        {
            for(j=0; j<taille*taille; j++)
            {
//la meme chose comme avant
                if((i>=0 && i<taille) || (i>(taille-1)*taille && i<(taille*taille) )||((i+1)%(taille)==0 || i%(taille)==0 )|| (grille[i/taille]+i%taille)->v!=0)

                {
                    if (i==j)
                        mat[i][j]=1;
                    else
                        mat[i][j]=0;

                }
                else
                {
                    //ici les coefficients sont definie par la valeurs de epsilon r qui entourne chaque points on a arrivé ici apres l'echantillonnage decalé
                    if(j==i+1)
                    {
                        mat[i][j]=0.5*((grille[i/taille]+i%taille)->epsilon+(grille[i/taille+1]+(i)%taille)->epsilon);



                    }
                    else if(j==i-1)
                    {

                        mat[i][j]=0.5*((grille[(i)/taille]+i%taille-1)->epsilon+(grille[(i)/taille+1]+(i)%taille-1)->epsilon);
                    }
                    else if(j==i+(taille))
                    {

                        mat[i][j]=0.5*((grille[(i)/taille]+i%taille)->epsilon+(grille[(i)/taille]+(i)%taille-1)->epsilon);
                    }
                    else if(j==i-taille)

                    {

                        mat[i][j]=0.5*((grille[(i)/taille+1]+(i)%taille)->epsilon+(grille[(i)/taille+1]+(i)%taille-1)->epsilon);
                    }
                    else if(j==i)
                    {

                        mat[i][j]=-((grille[j/taille]+j%taille)->epsilon+(grille[(j)/taille]+(j)%taille-1)->epsilon+(grille[(j)/taille+1]+(j)%taille)->epsilon+(grille[(j)/taille+1]+(j)%taille-1)->epsilon);
                    }
                    else
                    {
                        mat[i][j]=0;
                    }
                }
            }
        }
    }
//initialisation du seconde membre
    for(i=0; i<taille; i++)
    {
        for(j=0; j<taille; j++)
        {
            if((grille[i]+j)->rho!=0)
            {
                //le cas ou il y a une charge
                //on a changé la valeur de epsilon 0 et on a mltiplié ici par ce qui nou manque la puissance parceque avant il y avait une erreur numerique
                //et le second membre etait pas bien defini
                b[i*taille+j]=(-((grille[i]+j)->rho)*h*h)/permitiviteduvide ;
            }
            else
            {
                //si il y a pas de charge

                b[i*taille+j]=(grille[i]+j)->v;
            }
        }

    }

//on implimente notre methode de solution
    successive_over_relaxation(2601,mat,b,x,0.005);
    // gauss_seidel(taille*taille,mat,b,x,1000);


    for(i=0; i<taille; i++)
    {
        for(j=0; j<taille; j++)
        {
            (grille[i]+j)->v=x[i*taille+j];
        }
    }


}

///***************************************************************************************************////
/** La fonction capacité sert à calculer la capacité d'un condensateur                                ****
 * \param grille: notre monde                                                                         ****
 * \param  X :limite x de la structure
 * \param cordx: 2 eme limite x de la structure                                                       ****
 * \param cordy cordy2 : les deux limites sur les y                                                   ****
 * \param V1 V2 les potentiel de la structure (il faut avoir un dipole ou une diff de potentiel)      ****
 * \return la capacité calculé                                                                         ***
 */
//***************************************************************************************************////
double capacite(Espace**grille,int*x,int y,int *cordx2,int cordy2,double V1,double V2)
{

    double Energie=0;
    int min,max;
    double integral;
    double C;
    double module_E;
    int i=0,j,k=0;

    if(y>cordy2)
    {
        min=y;
    }

    else
    {
        min=cordy2;
        y=cordy2;
    }
    if(x[1]>cordx2[1])
    {
        min=cordx2[1];
    }
    else
    {
        max=x[1];
    }

    for(i=taille;i<min+2;i++)
    {
        integral=0;
        if((i>=x[0] || i>=cordx2[0]) && (i<=x[1] || i<=cordx2[1]))
            //si on est entre la limite de la structure en x
        {
            k=k+1;
            for(j=0;j<y-1;j++)
            {
                //entre la limite en y
                //le premier premier integral
                module_E=sqrt(pow(grille[i][j].E[0],2)+pow(grille[i][j].E[1],2));
                integral=integral+h * module_E *(grille[i][j].epsilon);


             }
            integral=integral+(h/2)*(sqrt(pow(grille[i][min+1].E[0],2)+pow(grille[i][y-1].E[1],2))*(grille[i][j].epsilon));

            if(k == 1 || k == max )
            {
                //pour les valeurs de l'extrimité

                Energie=(Energie +(h/2)*integral);
            }

            else
            {
                //calcul du deuxime integral on a remarqué que les valeurs qui sont au milieu avec la methode de trapez participe de fois donc on a multiplié par deux
//et on a enlevé la somme des deux termes
                Energie=(Energie+h*integral);
            }
        }

    }
//calcul de la capacité equivalente
    C=2*Energie*permitiviteduvide/(V1-V2)*(V1-V2);

    return C;
}

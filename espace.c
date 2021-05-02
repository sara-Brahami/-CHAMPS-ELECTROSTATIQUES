#include "espace.h"
///*******************************************************************************
/** La fonction Espace sert à creer notre espace et renvoie un matrice (grille)*
	  * @param sans paramètre                                                  *
	  * @return la grille de type pointeur                                     *
	  */
///*****************************************************************************
Espace** creerespace()
{
    int i;
    //creation d'un tableu 2D


    //1-Allocation dynamique du tableu colonne
    //bloc mémoire de pointeurs vers la structure espace
    Espace **grille=calloc(taille,sizeof(Espace*));

    //2-Allocation dynamique du tableau ligne
    for(i=0;i<taille;i++)
    {
        //chaque ligne (grille[i]) est un tableau de type struct Espace
        grille[i]=calloc(taille,sizeof(Espace));
    }

    return grille;
}





/** La fonction potentiel sert à donner pour chaque point de cordonner x et y un potentiel v **
     * @param grille : chaque point dans la grille on va lui donner un potentiel             **
     * @param v: la valeur du potentiel                                                      **
     * @param x,y:les cordonnées du point suivant l'axe X et Y                               **
     * @param v: la valeur de potentiel                                                      **
     * @param nbrpoint :nombre de points que l'utilisateur a choisi                          **
     * @return la fonction ne retourne rien car elle est de type void                        **
     */
//*********************************************************************************************


void voltage (Espace**grille,double v,int*x,int *y,int nbrpoint)
{

    int i;
    for(i=0; i<nbrpoint; i++)
    {

        // chaque point dans l'espace on lui affecte son potentiel par exemple v(1,2)=2 et vont etre presque à l'interieur de
        //la grille d'ailleurs c'est pour ça qu'on ajouté taille/2 ;
        grille[(taille/2)-x[i]][(taille/2)-y[i]].v = v;

    }

}



//******************************************************************************************************************
/* La fonction champs_electrique(Electric_field) sert à cacluler le champs électrique dans des différents points**
     * @param grille : il va contenir des champs positionné par l'utilisateur                                    **
     * @return la fonction ne retourne rien car elle est de type void                                            **
 */
//*****************************************************************************************************************
void Electric_fields(Espace**grille)
{

  int i,j;
    double Etempx[taille][taille];
    double Etempy[taille][taille];

    for(i=0; i<taille; i++)
    {
        for(j=1; j<taille-1; j++)

        {

            Etempx[i][j]=-(grille[i]+j)->v+(grille[i]+j-1)->v;
            Etempy[i][j]=-(grille[j]+i)->v+(grille[j-1]+i)->v;
        }

    }

//on fait une moyenne de deux du E calculer (une moyenne en y pour Ex et une moyenne en x pour Ey)

    for(i=1; i<taille-2; i++)
    {
        for(j=1; j<taille-2; j++)
        {




            (grille[i]+j)->E[0]=0.5*(Etempx[i][j+1]+Etempx[i][j]);
            (grille[i]+j)->E[1]=0.5*(Etempy[j+1][i]+Etempy[j][i]);
        }

    }

}


//********************************************************************************
/** La fonction " densité de charge " sert à affecter au centre de notre grille**
une densité de charge initiale  rho },                                         **
    * @param grille : la grille de dimension taille*taillle                    **
    * @param rho: la valeur de la densité de charge initiale                   **
    * @return la fonction ne retourne rien car elle est de type void           **
*/
//*******************************************************************************

void charge_density_center (Espace**grille,double rho)
{

   grille[taille/2][taille/2].rho = rho;


}


//*****************************************************************************************************
/** La fonction "charge_density_any_point " sert à affecter une densité de charge rho               **
à n'importe quel endroit dans la grille , en donnant les cordonnées  par l'utilisateur              **
    * @param grille : c'est l'espace                                                                **
    * @param point : un vecteur de 2 de deux cases qui contient les cordonnées d'un point quelconque**
    * @param rho: la valeur de la densité de charge donnée par l'utilisateur                        **
    * @return : la fonction ne return rien elle est de type void                                    **
*/
//*****************************************************************************************************

void charge_density_any_point( Espace**grille,int*point,double rho)
{
      //on place la chrage dans un un point choisi par l'utilisateur et décalé par exemple rho(pointx+centre,pointy+centre)
    (grille[(taille/2)+point[1]][(taille/2)+point[0]]).rho=rho;

}



//*******************************************************************************************************
/**  La fonction sert à  affecter  la densité  de charge donnée par l'utilsateur avec leur coordonnée ,**
 et nous on va les placer dans la grille.                                                              **
   * @param grille : notre espace                                                                      **
   * @param densite: la valeur de rho donnée par l'utilisateur                                         **
   * @param x,y: se sont des vecteurs qui contiennent les cordonnées de la densité de charge           **
   * @param nbrpoints :nombre de points                                                                **
   * @param type: condition pour l'option de remplissage                                               **
 */
 //******************************************************************************************************

void charge_density (Espace**grille,int*x,int*y,int nbrpoints, double densite,int type)
{
    int i=0,j,pos1,pos2;
    for(i=0;i<nbrpoints;i++)
    {
        grille[((taille-1)/2)-x[i]][((taille-1)/2)+y[i]].rho=densite;


    }

    if(type==1)
    {
        //on vas parcourir nontre monde (matrice)
        for(i=0;i<taille;i++)
        {
            pos1=0;
            pos2=0;
            for(j=0; j<taille; j++)
            {
                if(grille[i][j].rho!=0)
                {
                    if(pos1==0 && pos2==0)
                    {
                        //si on trouve une valeur de rho non nul on sauvegarde la position
                        pos1=j;
                    }
                    else if(pos2==0)

                    {
                        ///si on trouve une autre valeur on compare entre les deux rho si ils sont les memes donc c'est la meme forme on sauvegarde la pos2
                        if(grille[i][j].rho==grille[i][pos1].rho)
                        {
                            pos2=j;
                        }
                        else
                        {
                            ///si non on reinitialise pos1 a cette nouvelle valeur et on continue notre recherche pour un point sur la meme ligne avec le meme rho
                            pos1=j;
                        }
                    }

                }
                if(pos1!=0 && pos2!=0)
                {
                    ///si on trouve deux point qui realise les conditions en haut en remplit les valeus de rho de tout les point entre c'est deux
                    for(int k=pos1+1; k<pos2; k++)
                    {
                        grille[i][j].rho=grille[i][pos1].rho;
                    }
                    pos1=0;
                    pos2=0;
                }
            }///on fait ca pour toute les colonnes

        }
    }




}


///*********************************************************************
/** La fonction  on vas definir notre epsilon r appellé epsilon  a 1
 */
///**********************************************************************
void espace_epsilon(Espace**grille)
{
    int i,j;
    for(i=0; i<taille; i++)
    {
        for(j=0; j<taille; j++)
        {
            grille[i][j].epsilon=1;
        }
    }
}
///***************************************************************************************************************************************************************
/** \brief cette fonction affecte un epsilon r a un ensemble de point et si on veut elle remplie ce qu'il y a entre c'est point si ils sont fermé avec epsilon r **
 * \param space : notre espace
 * \param r: la valeur de epsilon r
 * \param cordx,xordy:les cordonnées du point
 * \param n :nombre de points
 * \param type: condition pour l'option de remplissage
 */
//****************************************************************************************************************************************************************

void permitivite_relative(Espace**grille,int*x,int*y,double n,double epsilonr,int type)
{
    //elle fonctionne de plus ou moins comme la fonction precedente
    int i,j=0,pos1,pos2,k;

    for(j=0;j<n;j++)
    {

   grille[(taille-1)/2-x[j]][(taille-1)/2+y[j]].epsilon=epsilonr;

    }

    if(type==1)
    {
        for(i=0; i<taille; i++)
        {
            pos1=0;
            pos2=0;
            for(j=0; j<taille; j++)
            {
                if(grille[i][j].epsilon!=1)
                {
                    if(pos1==0 && pos2==0)
                    {
                        pos1=j;
                    }
                    else if(pos2==0)

                    {
                        if(grille[i][j].epsilon==grille[i][pos1].epsilon)
                        {
                            pos2=j;
                        }
                        else
                        {
                            pos1=j;
                        }
                    }

                }
                if(pos1!=0 && pos2!=0)
                {
                    for(k=pos1+1; k<pos2; k++)
                    {
                        grille[i][k].epsilon=grille[i][pos1].epsilon;
                    }
                    pos1=0;
                    pos2=0;
                }
            }

        }
    }
}


///***************************************************************************************************************************////
/** La fonction Cercle sert à calculer les cordonnées d'un point sachant qu'on connait sont rayon et son centre                **
    *@ parma centre : le centre du cercle qui est vecteur contenant deux cases la case 1 contient les coord x et case 2 coord y **
    *@ param cordx : les coordonnées du point qu'on cherche à calculer
    *@ param r     : le rayon R
    @ return        : la fonctio  ne return rien
 */
 ///****************************************************************************************************************************///


void cercle(double *centre,int *cordx,int *cordy,int r)
{
    int i= 0;
    for(i=0;i<2*r;i++)
//on parametrise notre circle par c'est deux equations et on vas parcourir un des axe on calculent les cordonnées du deuxieme axe
//on fait ca sur x et y parceque on est pas en continue et donc si les cordonées calculés sont proche avec une diffirence moins de h
//il vont donner le meme point donc on evite ca par un parcourir les deux axes
    {
        cordx[i]=centre[0]+r-i;
        cordy[i]=sqrt(r*r-(cordx[i]-centre[0])*(cordx[i]-centre[0]))+centre[1];
        cordx[i+2*r]=centre[0]+r-i;
        cordy[i+2*r]=-(cordy[i]-2*centre[1]);

        cordy[i+4*r]=centre[1]+r-i;
        cordx[i+4*r]=(sqrt(r*r-(cordy[i+4*r]-centre[1])*(cordy[i+4*r]-centre[1]))+centre[0]);
        cordy[i+6*r]=centre[1]+r-i;
        cordx[i+6*r]=-(cordx[i+4*r]-2*centre[0]);



    }
}
void carre(double *centre,int *cordx,int *cordy,int *dimension)
{
    int i,j,k;
    int cond=1;
    //cette condition c'est pour assurer que les quatres bares de notre rectangle sont toujours dans le meme quadron
    if(centre[0]*centre[1]<=0)
        cond=-1;

    for(i=0; i<2; i++)
    {
        for(j=0; j<dimension[0]; j++)
        {
            cordx[j+i*dimension[0]]=(centre[0]+(-dimension[0]/2+j));
            cordy[j+i*dimension[0]]=(centre[1]+(-(dimension[1]/2)+i*(dimension[1]-1)));

        }


    }

    k=2*j;


    for(i=0; i<2; i++)
    {
        for(j=0; j<dimension[1]; j++)
        {
            cordy[k+j+i*dimension[1]]=cond*(centre[0]+(-dimension[1]/2+j));
            cordx[k+j+i*dimension[1]]=cond*(centre[1]+(-(dimension[0]/2)+i*(dimension[0]-1)));

        }


    }


}


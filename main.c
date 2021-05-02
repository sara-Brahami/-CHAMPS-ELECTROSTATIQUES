
/**************************************
*                                     *
*            L3 EEA                   *

*
*                                     *
*         potentiel et champs         *
*                                     *
*     Réalisé en 2020/2021            *
*         Sorbonne université(UPMC)   *
**************************************/

//les bibliothèques

#include <stdio.h>
#include <stdlib.h>

#include"numerique.h"      // Fichier d'en-tête pour les méthodes numériques( gausse seidel ...)
#include"espace.h"                    // Fichier d'en-tête pour espace


#define r 8 // le petit rayon
#define r1 18   //le grand rayon


//*****************************************************************
//****                ********************                     ****
//******            **                    **                  *****
//*********        ***PROGRAMME PRINCIPALE***            **********
//******            **                    **                  *****
//****                ********************                     ****
//*****************************************************************


int main()
{

     //Décalaration des variables
    int i,j;
    // allocation de matrice a///
double** a=calloc(taille*taille,sizeof(double*));

    for(i=0; i<taille*taille; i++)
    {
        a[i]=calloc(taille*taille,sizeof(double));

    }
    double *b=calloc(taille*taille,sizeof(double));
    double *X=calloc(taille*taille,sizeof(double));


    //le condensateur est constitué de deux armature en parallèle l'armature y1 et armature y2 d'une taille 20 
    int armature_y1[20]= {8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8};
    int armature_y2[20]= {-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9};
    int armature_x[20]=  {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9};
    double centre[2]= {0,0};   // Déclaration de la varaible centre 
    double epsilonr,rayon,RHO;   
    int choix;
    int nbrelement;
    double dim[2];
    int dim_calc[2];
    int *cordy, *cordx;
    double pointx,pointy;
    int pointsxd;
    int pointyd;
    int nbrpoints;
    int nbrcercle;
    double Vc;
    int* cordonxc_grand=calloc(8*r1,sizeof(int));
    int* cordonyc_grand=calloc(8*r1,sizeof(int));

    int* cordonxc_petit=calloc(8*r,sizeof(int));
    int* cordonyc_petit=calloc(8*r,sizeof(int));

//appel à des fonctions 
    Espace** mat=creerespace(); 
    espace_epsilon(mat);

     ///*************************************************************************************************************************

    printf("********Bienvenue dans notre project Calcul du champ et du potentiel electrostatique dans le vide et dans un dielectrique**********\n");
     ///**************************************************************************************************************************



      ///****************************************************************************************************************************
                   //*le choix des types d'affichage *//
      ///***************************************************************************************************************************


      do {
          printf("si vous souhaitez avoir  :\n1*des cas standard deja definie (condensateur plan ou condesateur circulaire) Merci de taper 1!\n2*si vous voulez des cas ou vous choisissiez tous les parametres du simulation merci de taper 2!\n");
          if (scanf("%d", &choix)!= 1) {
          printf("Erreur fatale à la saisie");
          exit(1);
          }
          if((choix!=1)&&(choix!=2))
          {
          printf("Merci de resayer une autre fois car vous n avez pas respecter les consignes!!! \n");
          }

       } while ((choix !=1)&&(choix!=2)) ;
    switch(choix)
    {
    case 1:
    {
        do
        {
               printf("si vous souhaitez \n1.1 **un condensateur plan** Merci de taper 3\n1.2**un condestateur circulaire** Merci de taper 4\n");
          if (scanf("%d", &choix) != 1) {
          printf("Erreur fatale à la saisie");
          exit(1);
          }
          if((choix!=3)&&(choix!=4))
          {
          printf("Merci de resayer une autre fois car vous n'avez pas respecter les consignes!!! \n");
          }
        }
        while(choix!=3 && choix!=4);
        switch(choix)
        {
            ///********cas condensateur plan *******///
        case 3:
        {
        do
        {
            printf("vous avez choisi un condensateur plan\n");
            printf("si vous souhaitez avoir:\n1.1.1**avec un dielectrique dans le condesateur plan**Merci de taper 5\n1.1.2**Sans dielectrique**Merci de taper 6\n");
            scanf("%d",&choix);

            if((choix!=5)&&(choix!=6))
            {
             printf("Merci de resayer une autre fois car vous n'avez pas respecte les consignes!!! \n");
            }
        }while((choix!=5) && (choix!=6));
            switch(choix)
            {
                ///*****condensteur avec un dielectrique****///
            case 5:
                printf("vous avez choisi *un condensateur avec un dielectrique*\n");
                printf("votre condesateur a les dimensions:\nl=%f\n d=%f\n",17*h,15*h);
                printf("choisissiez les dimensions du dielectrique\n");

                printf("Merci d'entrer sa longeur\n");
                scanf("%lf",&dim[0]);


                printf("Merci d'entrer sa largeur\n");
                scanf("%lf",&dim[1]);

                dim_calc[0]=dim[0]/h;
                dim_calc[1]=dim[1]/h;

                printf("Merci d'entrer une valeur de la permitivite relative epsilonr=\n");
                scanf("%lf",&epsilonr);

                printf("Merci d'entrer la valeur du potentiel Vc =\n");
                scanf("%lf",&Vc);

                int* cordonx=calloc(2*(dim_calc[1]+dim_calc[0]),sizeof(int));
                int* cordony=calloc(2*(dim_calc[1]+dim_calc[0]),sizeof(int));

                carre(centre,cordonx,cordony,dim_calc);

                permitivite_relative(mat,cordony,cordonx,2*(dim_calc[1]+dim_calc[0]),epsilonr,1);

                voltage(mat,-Vc/2,armature_x,armature_y1,20);
                voltage(mat,-Vc/2,armature_x,armature_y2,20);
                break;

                ///******condensateur sans dielectrique *******////
            case 6:
                printf("vous avez choisi un condensateur plan sans dielctrique \n");
                printf("votre condesateur a les dimensions suivantes:\nl=%lf\nd=%lf\n",17*h,15*h);
                printf("Merci d'entrer la valeur du potentiel Uc =\n");
                scanf("%lf",&Vc);
                voltage(mat,-Vc/2,armature_x,armature_y1,20);
                voltage(mat,-Vc/2,armature_x,armature_y2,20);
                break;

            default :
                printf("I dont know what you entered");
                break;
            }



        }

        break;
        ///**********cas condensateur circulaire( avec dielectrique ou sans dielectrique) ********///
        case 4:
        {
            do
            {
            printf("vous avez choisi un condensateur circulaire\n");
            printf("si vous souhaitez avoir du:\n2.1*dielectrique dans le condensateur* Merci de taper 7 \n2.2 sans dielctrique Merci du taper 8\n");

             scanf("%d",&choix);
             if((choix!=7)&&(choix!=8))
             {
                 printf("Merci de resayer une autre fois car vous n'avez pas respecte les consignes!!! \n");
             }
            }while(choix!=7 && choix!=8);

            switch(choix)
            {
               ///****condensateur avec dielectrique ****////
            case 7:
                printf("votre condesateur a les dimensions:\nGrand rayon=%lf\nPetit Rayon=%lf\n",r1*h,r*h);
                printf("choisissiez les dimensions de votre dielectrique\nGrande rayon: ");
                scanf("%lf",&dim[0]);
                printf("Petit rayon: ");
                scanf("%lf",&dim[1]);
                dim_calc[0]=dim[0]/h;
                dim_calc[1]=dim[1]/h;
                printf("Merci de bien vouloir entrer la valeur de la permitivite relative ");
                scanf("%lf",&epsilonr);
                printf("Merci de bien vouloir entrer la valeur de la tension Vc");
                scanf("%lf",&Vc);


                int* cordonxgrand=calloc(8*dim_calc[1],sizeof(int));
                int* cordonygrand=calloc(8*dim_calc[1],sizeof(int));

                int* cordonxpetit=calloc(8*dim_calc[0],sizeof(int));
                int* cordonypetit=calloc(8*dim_calc[0],sizeof(int));
                //petit cercle //
                cercle(centre,cordonxc_petit,cordonyc_petit,r);
                //grand cercle //

                cercle(centre,cordonxc_grand,cordonyc_grand,r1);
                cercle(centre,cordonxgrand,cordonygrand,dim_calc[1]);
                cercle(centre,cordonxpetit,cordonypetit,dim_calc[0]);
                permitivite_relative(mat,cordonxgrand,cordonygrand,8*(dim_calc[1]),epsilonr,1);
                permitivite_relative(mat,cordonxpetit,cordonypetit,8*(dim_calc[0]),epsilonr,1);

                voltage(mat,Vc/2,cordonxc_grand,cordonyc_grand,8*r1);
                voltage(mat,-Vc/2,cordonxc_petit,cordonyc_petit,8*r);
                free(cordonxpetit);
                free(cordonxgrand);
                free(cordonxc_grand);
                free(cordonxc_petit);

                free(cordonxpetit);
                free(cordonxgrand);
                free(cordonyc_grand);
                free(cordonyc_petit);
                break;
                ///******condensateur sans dielectrique ******///
            case 8:
                printf("vous avez choisi un condensateur circulaire sans dielctrique a l'interieur");
                printf("votre condesateur a les dimensions:\nGrand rayon=%lf\nPetit Rayon=%lf\n",r1*h,r*h);
                printf("Merci d'entrer la valeur de la tension Vc:");
                scanf("%lf",&Vc);

                cercle(centre,cordonxc_petit,cordonyc_petit,r);
                cercle(centre,cordonxc_grand,cordonyc_grand,r1);

                voltage(mat,Vc/2,cordonxc_grand,cordonyc_grand,8*r1);
                voltage(mat,-Vc/2,cordonxc_petit,cordonyc_petit,8*r);

                free(cordonxc_grand);
                free(cordonxc_petit);


                free(cordonyc_grand);
                free(cordonyc_petit);
                break;
             default :
                printf("I dont know what you entered");
                break;

            }
        }
        break;
        default :
            printf("I dont know what you entered");
        break;

        }
    }
    break;
     ///***********avec densité de charge ou bien sans densité de charge ********///
    case 2:
    {      printf("vous etes dans le cas ou c'est a vous de choisir le type de simulation\n");
          do
          {
            printf("si vous souhaitez avoir \n1*Une source de tension*Merci de taper 1!\n2*Une distribution de charges*Merci de taper 2!\n");

            if(scanf("%d",&choix)!=1)
            {
                printf("Erreur fatale à la saisie\n");
                exit(1);
            }
            if((choix!=1)&&(choix!=2))
            {
             printf("Merci de resayer une autre fois car vous n'avez pas respecté les consignes!!! \n");
            }


          }
          while(choix!=1 && choix!=2);
        switch(choix)
        {
            ///*****************************************************************///
                           ///sans densité de charge ///
           ///*****************************************************************///
        case 1:
              printf("vous etes dans le cas tension \n");
            do
            {
                printf("on a 2 sources diponible si vous souhaitez :\n1*des cercles,Merci de taper 1,\n2 *des  points,Merci de taper 2:\n");
                if(scanf("%d",&choix)!=1)
                {
                    printf("Erreur fatale à la saisie\n");
                    exit(1);
                }
                if((choix!=1)&&(choix!=2))
                {
                  printf("Merci de resayer une autre fois car vous n'avez pas respecté les consignes!!! \n");
                }
            }while((choix!=1)&&(choix!=2));

                switch(choix)
                {
                   ///*********************************************///
                          ///cas des cercles ///
                  ///*********************************************///
                case 1:
                {

                    printf("cas:cercle\n");
                    printf("\n Merci d'entrer le nombre de cercle que vous souhaitiez utiliser");
                    scanf("%d",&nbrcercle);
                    printf("\n prenez en consideration que la simulation elle se fait dans l'espace lmiter par\n|xmax|=%lf\n|ymax|=%lf",taille*h/2,taille*h/2);

                    for(i=1; i<=nbrcercle; i++)
                    {
                        printf("\n Merci d entrer les coordonnees du centre de cercle numero %d\n",i);
                        printf("\n veuillez entrer les coordonnees suivant l'axe x ");
                        scanf("%lf",&centre[1]);

                        printf("\n veuillez entrer les coordonnées suivant l'axe y");
                        scanf("%lf",&centre[0]);


                        centre[0]= -centre[0]/h;
                        centre[1]= centre[1]/h;

                        printf("\nveuillez entrer  le rayon cercle numero %d: ",i);
                        scanf("%lf",&rayon);

                        rayon=rayon/h;

                        printf("\n veuillez entrer le potentiel du circle numero %d: ",i);
                        scanf("%lf",&Vc);

                        cordx=calloc(8*rayon,sizeof(int));
                        cordy=calloc(8*rayon,sizeof(int));


                        cercle(centre,cordx,cordy,rayon);
                        voltage(mat,Vc,cordx,cordy,8*rayon);

                        free(cordx);
                        free(cordy);

                    }

                }
                break;
                ///*****************************///
                    ///cas des points ///
                ///*****************************///
                case 2:
                {
                    printf(" Merci d'entrer le nombre de point que vous souhaitiez utiliser");
                    scanf("%d",&nbrpoints);
                    printf("prenez en consideration que la simulation elle se fait dans l'espace lmité par\n|xmax|=%lf\n|ymax|=%lf",taille*h/2,taille*h/2);

                    for(i=1; i<=nbrpoints; i++)
                    {
                        printf("\nveuillez entrer les coordonnées suivant l'axe x du point %d=",i);
                        scanf("%lf",&pointx);
                        printf("\nveuillez entrer les coordonnées suivant l'axe y du point %d=",i);
                        scanf("%lf",&pointy);
                        pointsxd=-(int)pointy/h;
                        pointyd=(int)pointx/h;

                        printf(" \nveuillez entrer le potentiel Vc du point %d =",i);
                        scanf("%lf",&Vc);

                        voltage(mat,Vc,&pointsxd,&pointyd,1);
                    }
                }
                break;
                default :
                   printf("I dont know what you entered");
                break;

                }

          ///*******************************************************************////
                        /// cas ou on a la densité de charge rho//
         ///*******************************************************************////
        case 2:
        {
            do
            {
              do
              {
               printf("\n aimeriez vous avoir une densite de charge avec des cercles ou des points \n");
               printf(" \n si vous deseriez \n1 *des cercles Merci de taper 1!\n2 *des points Merci de taper 2! ");
               if(scanf("%d",&choix)!=1)
               {
                   printf("Erreur Fatale à la saisie\n");
                   exit(1);
               }
               if((choix!=1)&&(choix!=2))
               {
                   printf("merci de retaper une autre fois car vous n avez pas respecter les consignes\n");
               }
               }while((choix!=1)&&(choix!=2));
                switch(choix)
                {
                 ///***cas  des cercle avec une densité de charge ****.////
                case 1:
                {


                    printf("\n veuillez entrez le nombre de cercle que vous voulez utiliser:\n: ");
                    scanf("%d",&nbrcercle);
                    printf("\nprenez en consideration que la simulation elle se fait dans l'espace lmite par\n|xmax|=%lf\n|ymax|=%lf\n",taille*h/2,taille*h/2);
                    \
                    for(i=1; i<=nbrcercle; i++)
                    {
                        printf("\nMerci d'entrer les cordonnees du cecle numero%d:\n",i);
                        printf("\nveuillez entrer les coordonnees suivant l'axe X=\n");
                        scanf("%lf",&centre[0]);
                        printf("\nveuillez entrer les coordonnees suivant l'axe Y=\n");
                        scanf("%lf",&centre[1]);
                        centre[0]=-(int)centre[0]/h;
                        centre[1]=(int)centre[1]/h;

                        printf("\nveuillez entrer le rayon du circle numero %d",i);
                        scanf("%lf",&rayon);
                        rayon=(int)rayon/h;

                        printf("\nveuillez entrer la densite de charge du cercle numero %d",i);
                        scanf("%lf",&RHO);

                    ///****** choix entre une densité surfacique ou bien lineique*******///
                       printf("\n aimeriez vous ajouter une densite surfacique ou bien lineique sur votre cercle");
                       do
                        {
                        printf("\n si vous souhaitez :\n1 *une densite surfacique* Merci de taper 1!\n2 *une densite lineique* Merci de taper 2!\n");
                        if (scanf("%d", &choix) != 1) {
                        printf("Erreur fatale à la saisie\n");
                        exit(1);
                        }
                        if((choix!=1)&&(choix!=2))
                        {
                        printf("Merci de resayer une autre fois car vous n'avez pas respecté les consignes!!! \n");
                        }
                        }while((choix!=1)&&(choix!=2));

                        cordx=calloc(8*rayon,sizeof(int));
                        cordy=calloc(8*rayon,sizeof(int));


                        cercle(centre,cordx,cordy,rayon);
                        charge_density (mat,cordx,cordy,8*rayon,RHO,choix);



                          	//libération de notre tableau //
                        free(cordx);
                        free(cordy);

                    }

                }
                break;


                case 2://point
                {
                    printf(" Merci d'entrez le nombre de points que vous voulez utiliser:\n");
                    scanf("%d",&nbrelement);
                    printf("prenez en consideration que la simulation elle se fait dans l'espace lmité par\n|xmax|=%lf\n|ymax|=%lf\n",taille*h/2,taille*h/2);

                    for(i=1; i<=nbrelement; i++)
                    {
                        printf("\nentrer les cordonnes  du point %d \n",i);
                        printf("cordx: \n");
                        scanf("%lf",&pointx);
                        printf("cordy: \n");
                        scanf("%lf",&pointy);
                        pointsxd=-(int)pointy/h;
                        pointyd=(int)pointx/h;
                        printf("entrer le potentiel du point %d : \n",i);
                        scanf("%lf",&Vc);
                        voltage(mat,Vc,&pointsxd,&pointyd,1);

                    }
                }
                break;

                }


                printf("Est-ce que vous voulez ajoutez d'autres structures:\n");
                printf("si oui merci de tapez 1 sinon tapez 2\n");

                scanf("%d",&choix);


            }
            while(choix==1);

        }
        break;
        }
        printf("aimeriez vous ajouter des zones avec des differents dielectriques\n");
         do
         {
           printf("si vous le souhaitez merci de taper 1 sinon taper 0\n");
           if(scanf("%d",&choix)!=1)
            {
            printf("erreur fatale a la saisie \n");
            }
           if((choix!=1)&&(choix!=0))
            {
            printf("Merci de ressayer une autre fois\n");
            }
          }while((choix!=1)&&(choix!=0));
        if(choix==1)
        {
                printf("on a comme formes des cercles \n");

              do
              {


                    printf("Merci d'entrer le nombre de cercle que vous souhaitiez utilisez !\n");
                    scanf("%d",&nbrcercle);
                    printf("prenez en consideration que la simulation elle se fait dans l'espace lmité par\n|xmax|=%lf\n|ymax|=%lf",taille*h/2,taille*h/2);

                    for(i=1; i<=nbrcercle; i++)
                    {
                        printf("Merci d'entrer les coorconnées du centre du cercle numéro %d=\n",i);
                        printf("veuillez entrer les coordonnées du cercle suivant l'axe x=\n");
                        scanf("%lf",&centre[0]);

                        printf("veuillez entrer les coordonnées du cercle suivant l'axe y=\n");
                        scanf("%lf",&centre[1]);

                        centre[0]=-(int)centre[0]/h;
                        centre[1]=(int)centre[1]/h;

                        printf("Merci d'entrer le rayon du circle numero %d=\n",i);
                        scanf("%lf",&rayon);
                        rayon=(int)rayon/h;

                        cordx=calloc(8*rayon,sizeof(int));
                        cordy=calloc(8*rayon,sizeof(int))
                        ;
                        printf("Merci d'entrer la permitivité relative du cercle numéro %d=\n",i);
                        scanf("%lf",&epsilonr);


                        printf("aimeriez vous d'ajouter un epsilon r sur le contour ou sur toute la surface de votre cercle:\n");
                        do
                        {
                        printf("si vous souhaitez de l'ajouter sur \n 1*la surface* merci de taper 1! ,\n2-contour Merci de taper 2\n");
                         if(scanf("%d",&choix)!=1);
                         {
                             printf("Erreur fatale à la saisie\n");
                             exit(1);
                         }
                         if((choix!=1)||(choix!=2))
                         {
                             printf("erreur ! merci de ressayer une autre fois\n");
                         }
                        }while ((choix!=1)||(choix!=2));

                        cercle(centre,cordx,cordy,rayon);
                        epsilon_r(mat,cordx,cordy,8*rayon,epsilonr,choix);

                        free(cordx);
                        free(cordy);

                    }


                printf("est-ce que vous voulez d'autres zone \n1*Oui* Merci de taper 1\n2*Non* Merci de taper 2\n ");
                if(scanf("%d",&choix)!=1);
                {
                    printf("erreur fatale à la saisie");
                    exit(1);
                }
                if((choix!=1)&&(choix!=2))
                {
                    printf("erreur!merci de resayer une autre fois\n");
                }
            }while(choix==1);

        }
        else
        {
            printf("vous avez decidez de ne pas ajouter des zones avec de dielectriques\n");
        }
        printf("fin");
    }
    break;


}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



printf("si tu souhaite visualiser l affichage sur MATLAB le potentiel et du champ tu pouurais y acceder au fichier *champtension*\n);
printf("good bye");
// appel à la mtrice pour remplir mat
matrices(mat,b,a,X,1);

    FILE*p=fopen("tension.txt","w");
//on sauvegarde les valeurs dans des fichiers pour l'affichage
Electric_fields(mat);


    for(i=0; i<taille; i++)
    {
        for(j=0; j<taille; j++)
        {
            fprintf(p,"%lf",(mat[i]+j)->v);
        }
        fprintf(p, "\n");
    }
    fclose(p);


 FILE*pointeur=fopen("champs.txt","w");

    for(int k=0; k<2; k++)
    {

        for(i=0; i<taille; i++)

        {
            for(j=0; j<taille; j++)
            {
                fprintf(pointeur,"%lf,",(mat[i]+j)->E[k]);
            }
            fprintf(pointeur,"\n");
        }
    }
    fclose(pointeur);

//liberation de la mémoire
   // libérer secon memebre
    free(b);
    // on commence par libérer les colonne
    for(i=0; i<taille*taille; i++)
    {
        free(a[i]);
    }
    // liberer les lignes
    free(a);

    free(X);

    for(i=0; i<taille; i++)
    {
        free(mat[i]);

    }
    free(mat);




//liberation de l'espace
return 0;

}


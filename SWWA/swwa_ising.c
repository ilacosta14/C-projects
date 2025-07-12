// ************************************************************************************************
// SWENDSEN WANG multicluster code for the Ising Model
// see https://opus4.kobv.de/opus4-zib/files/4218/wende_steinke.pdf for reference
// ************************************************************************************************


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rng.h"

Uint RNG_seed;
Uint L; //volume of the lattice

int i,j;
int **s;   // spin "field"
double* two_pt;
Scalar g;  //Temperature (for our case T=g)
Scalar action,spin,sq_spin;
Scalar sum_spin=0,sum_sq_spin=0;
Scalar magn=0,susc=0;
int N=0; //number of configs

/** restart the RNG with the current seed **/
void RNG_increment( Uint increment ) {

        RNG_seed += increment;
}

/** restart the RNG with the current seed **/
void RNG_restart(void) {

        rng_set_seed( RNG_seed );
}
void initializelattice(){
	s=(int**)calloc(L,sizeof(int*));
    for(j=0;j<L;j++){
        s[j]=(int*)calloc(L,sizeof(int));
        }
    for(i=0;i<L;i++){
        for(j=0;j<L;j++){

            Scalar r=rng_get();
            if(r<0.5){
                s[i][j]=1;
            }
            else s[i][j]=-1; 
            // s[i][j]=1;
           }
    }
}

int **ibond,**jbond;    //lattice bonds
Scalar prob;            //active bonds probability
int **cluster;          //cluster labels
int *biglabel;        // biggest label on the label matrix
int *snewchosen;        // has the spin value been chosen?
int *snew;              //random spin value in each cluster

void initializevariables(){
    //allocate arrays for bonds in x and y directions
	ibond=(int**)calloc(L,sizeof(int*));
    for(j=0;j<L;j++){
        ibond[j]=(int*)calloc(L,sizeof(int));
    }
    jbond=(int**)calloc(L,sizeof(int*));
    for(j=0;j<L;j++){
        jbond[j]=(int*)calloc(L,sizeof(int));
    }
    
    //compute active bond probability
    prob=1.0-exp(-4.0/g);

    biglabel=(int*)calloc(L*L,sizeof(int));
    snewchosen=(int*)calloc(L*L,sizeof(int));
    snew=(int*)calloc(L*L,sizeof(int));
    
    //allocate 2D array for spin cluster labels
    cluster=(int**)calloc(L,sizeof(int*));
    for(j=0;j<L;j++){
        cluster[j]=(int*)calloc(L,sizeof(int));

    //allocate the two pt vector
    two_pt=(double*)calloc(L,sizeof(double*));
    }
}

void bonds(){
    //visit all the spins in the lattice
    for (i=0;i<L;i++){ 
	   for (j=0;j<L;j++){

          ibond[i][j]=jbond[i][j]=FALSE; //sets all the active bonds to be FALSE by default
		  int inext = i == L-1 ? 0 : i+1; //NN on the right
          
          Scalar r = rng_get();

		  if (s[i][j]==s[inext][j] && r<prob){ // Assigns an active bond on the right if the condition is satisfied
			ibond[i][j]=TRUE;
		  }

		  int jnext = j == L-1 ? 0 : j+1; //same as before but looking down
          
          r=rng_get();
		  
          if (s[i][j]==s[i][jnext] && r<prob){
			jbond[i][j]=TRUE;
		  }
	   }
    }
}

int properlabel(int label) {
    while (biglabel[label] != label)
        label = biglabel[label];
    return label;
}

void labelclusters(){ //this function is to create the labels matrix

	int label=0;
	for (int i=0;i<L;i++){
		for (int j= 0; j<L; j++){
			int bonds=0; //sets bonds=0 at every loop
			int bondi[4], bondj[4];
            //create two arrays that save the active bond for the site [i][j] and counts the number of active bonds per site
			if (i>0 && ibond[i-1][j]){ //if there is an active bond between i-1 and i (on the left)
				bondi[bonds]=i-1;              //saves the bonds of the point [i][j]
				bondj[bonds++]=j;
			}
			if (i==L-1 && ibond[i][j]) {
                bondi[bonds] = 0;                         // same but on the boundary because of periodic boundary conditions
                bondj[bonds++] = j;
            }

            if (j>0 && jbond[i][j-1]){
				bondi[bonds]=i;              //same as before but looking up
				bondj[bonds++]=j-1;
			}
			if (i==L-1 && jbond[i][j]) {
                bondi[bonds] = i;                         
                bondj[bonds++] = 0;
            }


            //At this point, I assign the labels, starting with the bonds
            // The following part of the code uses the Hoshen-Kopelman algorithm for labelling, see https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html

            if (bonds==0){ // if the point [i][j] has no active bond, then assign a cluster label to this site and increase the label
        	   cluster[i][j]=label;
        	   biglabel[label] = label;
        	   ++label;
        
            }

            else{
        	   int minlabel= label, plabel;
                for(int b=0;b<bonds;b++){
        	       plabel=properlabel(cluster[bondi[b]][bondj[b]]); //check all the neighbors
        	       if (minlabel>plabel){
        		      minlabel=plabel;
        	       }
                }

        	cluster[i][j]=minlabel; // assign the minimum label to all the elements of the cluster
        	   for(int b=0;b<bonds;b++){
        		  plabel=cluster[bondi[b]][bondj[b]];
                  biglabel[plabel]=minlabel;
        		  cluster[bondi[b]][bondj[b]] = minlabel;
        	
        	    }

            }
        }

	}
}

// Once everything is labeled, I flip the spin

void flipclusterspin(){
	for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++) {

            // random new cluster spins values have not been set
            int n = i * L + j;
            snewchosen[n] = FALSE;

            // replace all labels by their proper values
            cluster[i][j] = properlabel(cluster[i][j]);
        }   
    }
    int flips = 0;    // count number of spins that are flipped
    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++) {
            // find the now proper label of the cluster
            int label = cluster[i][j];

            // choose a random new spin value for cluster
            // only if this has not already been done
            if (!snewchosen[label]) {  
                Scalar r =rng_get();  
                snew[label] = r<0.5 ? +1 : -1;
                snewchosen[label] = TRUE;
            }

            // re-set the spin value and count number of flips
            if (s[i][j] != snew[label]) {
                s[i][j] = snew[label];
                ++flips;
            }
        }
    }
}


void MCstep(){
    bonds();
    labelclusters();
    flipclusterspin();
}

void measureObservables() {
    char file_name1[128],file_name2[128],file_name3[128],file_name4[128];
    FILE *fp1,*fp2,*fp3,*fp4;
    action=spin=sq_spin=0;
    snprintf(file_name1,sizeof(file_name1),"/home/ilaria/Dottorato/Coding/SWWA/%ldx%ld/action_%0.1f.dat",L,L,g);
    snprintf(file_name2,sizeof(file_name2),"/home/ilaria/Dottorato/Coding/SWWA/%ldx%ld/spin_%0.1f.dat",L,L,g);
    snprintf(file_name3,sizeof(file_name3),"/home/ilaria/Dottorato/Coding/SWWA/%ldx%ld/sq_spin_%0.1f.dat",L,L,g);
    snprintf(file_name4,sizeof(file_name4),"/home/ilaria/Dottorato/Coding/SWWA/%ldx%ld/two_pt_%0.1f.dat",L,L,g);
    // compute the action density and print it on file
    // printf("File Name: %s\n", file_name1);
    fp1 = fopen(file_name1, "a");
    if (fp1 == NULL) {
        fprintf(stderr, "Error opening file for writing.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++) {
            int inext = i == L-1 ? 0 : i+1;
            int jnext = j == L-1 ? 0 : j+1;
            action += s[i][j]*(s[inext][j] + s[i][jnext]);
        }
   }
   // action=action/(L*L);
   // printf("Value of g: %lf\n", g);
   action=-2/(g)*(action/(L*L));
   fprintf(fp1, "%lf\n", action);
   fclose(fp1);
   // compute the spin density and print it on file
   fp2=fopen(file_name2,"a");
   for (int i = 0; i < L; ++i){
       for (int j = 0; j < L; ++j){
           spin+=s[i][j];
       }
   }
   spin=spin/(L*L);
   fprintf(fp2, "%lf\n", spin);
   fclose(fp2);
   // compute the susceptibility density and print it on file
   fp3=fopen(file_name3,"a");
   sq_spin=spin*spin;
   fprintf(fp3, "%lf\n", sq_spin);
   fclose(fp3);
   // compute two pt function and print it on file
   fp4=fopen(file_name4,"a");
   for (int i = 0; i < L; ++i)
   {
       two_pt[i]=0;
       for (int j = 0; j < L; ++j)
       {
           two_pt[i]+=s[i][j]*s[0][0];
       }
       two_pt[i]=two_pt[i]/(L);
       fprintf(fp4, " %lf", two_pt[i]);
   }
   fprintf(fp4, "\n");
   fclose(fp4);
   // save the sum of spins, square of spins and number of configs for the avg
   sum_spin+=spin;
   sum_sq_spin+=sq_spin;
   ++N;
}
void compute_averages(){
    magn=sum_spin/N;
    susc=sum_sq_spin/N;
}


int main (){
    int steps;
    
    printf("Insert: seed   L  g   MC steps:\n");
    scanf("%ld %ld %lf %d",&RNG_seed,&L,&g,&steps);
    RNG_restart();
    printf( "L = %ld, g = %lf\n",L,g);
    // printf("Insert the size of the lattice\n");
    initializelattice();
    initializevariables();
    // g= g/2;
    int Thermsteps= steps/1000;
    printf("Performing thermalization steps...\n");
    for (int k = 0; k<Thermsteps;k++){
        RNG_increment(+2);
        RNG_restart();
    	MCstep();
 	
    }
    printf("Done\nPerforming MC steps and computing the observables...\n");
    for (int k = 0; k<steps;k++){
        RNG_increment(+2);
        MCstep();
        measureObservables();
    
    }
// printf ("Output config:\n");
//     for(i=0;i<L;i++){
//         for(j=0;j<L;j++){
//             printf("%d      ",s[i][j]);
//         }
//         printf("\n");
//     }
//     printf("\n"); 

    printf("Done\n");
    printf("Computing the averages:\n");
    compute_averages();
    printf("magn=%0.12f\n",magn);
    printf("susc=%0.12f\n",susc);
    printf("Done\n");
    return 0;
}
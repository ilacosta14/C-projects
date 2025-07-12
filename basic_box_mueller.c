#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES // for C
#include <math.h>
#include <time.h>

int main()
{
	double v[5000],u[5000],p[10000];
	char file_name[32];
    FILE * fpointer;
    sprintf(file_name,"box_muller.dat");
    fpointer= fopen(file_name,"w");

	for (int i = 0; i < 5000; ++i)
	{
		
		srandom(time(NULL)^random());
        v[i]=(double)rand()/RAND_MAX;
        u[i]=(double)rand()/RAND_MAX;
    	p[i]=sqrt(-2*log(1-v[i]))*cos(2*M_PI*(1-u[i]));
    	p[5000+i]=sqrt(-2*log(1-v[i]))*sin(2*M_PI*(1-u[i]));
	}
for(int i = 0; i < 10000; i++){

    fprintf(fpointer, "%f\n ", p[i]); 
    }

    fclose (fpointer);
	return 0;
}
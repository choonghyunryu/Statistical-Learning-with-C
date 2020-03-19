/* HMEAN.C */

#include <stdio.h>
#include <stdlib.h>

#define MAX 500
#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))

int n;
double total,x[MAX];

void main(int argc,char *argv[])
{
  FILE *stream;
  double *px,num,h_temp=0;

  px=x;

  if (argc<=1) {
    puts("Usage : HMEAN Datafile");
    exit(EXIT_FAILURE);
  }
  stream=fopen(argv[1],"rt");
  if (stream==NULL) printerrmsg("File not found!!!");

  while (fscanf(stream,"%lf\n",&num)!=EOF) {
    ++n;
    *(px+n)=num;
    if (*(px+n)==0) printerrmsg("Hamonic Mean has not 0 data.");
    else
      h_temp+=1./(*(px+n));
    total+=*(px+n);
  }

  printf("N               : %d \n"
         "Total           : %f\n"
         "Harmonic  Mean  : %f\n",n,total,n/h_temp);
}


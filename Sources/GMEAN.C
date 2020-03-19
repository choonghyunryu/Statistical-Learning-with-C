/* GMEAN.C */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 500
#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))

int n;
double total,x[MAX];

void main(int argc,char *argv[])
{
  FILE *stream;
  double *px,num,g_temp=0;
  px=x;
  if (argc<=1) {
    puts("Usage : GMEAN Datafile");
    exit(EXIT_FAILURE);
  }
  stream=fopen(argv[1],"rt");
  if (stream==NULL) printerrmsg("File not found!!!");

  while (fscanf(stream,"%lf\n",&num)!=EOF) {
    ++n;
    *(px+n)=num;
    if (*(px+n)<0) printerrmsg("Geometric Mean has not minus data.");
    else {
      if (*(px+n)==0) ; /*¸aža· ˆt·¡ 0·¡¡e Skip */
      else
        g_temp+=log(*(px+n));
    }
    total+=*(px+n);
  }
  printf("N               : %d \n"
         "Total           : %f\n"
         "Geometric  Mean : %f\n",n,total,exp(g_temp/n));
}


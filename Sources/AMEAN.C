/* AMEAN.C */

#include <stdio.h>
#include <stdlib.h> /* defined EXIT_FAILURE 1 */

#define MAX 500 /* ÂA” ¸aža ·³b® */
#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))

int n;
double total,x[MAX];

void main(int argc,char *argv[])
{
  FILE *stream;
  double *px,num;

  px=x;

  if (argc<=1) {
    puts("Usage : AMEAN Datafile");
    exit(EXIT_FAILURE);
  }
  stream=fopen(argv[1],"rt");
  if (stream==NULL) printerrmsg("File not found!!!");

  while (fscanf(stream,"%lf\n",&num)!=EOF) { /* EOFŒa»¡  ¸ažaŸi ·ª·q */
    ++n; /* ¸aža· ˆ•® ¬eÂ‰ */
    *(px+n)=num;
    total+=*(px+n); /* ˆbˆb ¸aža· Ðs ‰¬e */
  }

  printf("N               : %d\n"
         "Total           : %f\n"
         "Arithmetic Mean : %f\n",n,total,total/n);
}


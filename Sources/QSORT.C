/* QSORT.C */

#include <stdio.h>
#include <stdlib.h>

#define MAX 500  /* ÀáŸ¡Ði ¸aža· ÂA” ˆ•® */
#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))

void sort(int left,int right);

double *px,x[MAX];  /* ¸ažaŸi ¸á¸wÐi ¤µi */

void main(int argc,char *argv[])
{
  FILE *stream;
  int n=0,i;
  double num;

  px=x;

  if (argc<=1) {
    puts("Usage : QSORT Datafile");
    exit(EXIT_FAILURE);
  }
  stream=fopen(argv[1],"rt");                         /* Ìa·© µ¡Ïe */
  if (stream==NULL) printerrmsg("File not found!!!"); /* ¸ažaÌa·©·¡ ¸i¡µ ·³b */
                                                      /* –A´ö·i ˜            */ 
  while (fscanf(stream,"%lf\n",&num)!=EOF) {          /* Ìa·©· {Œa»¡ ¤e¥¢   */
    ++n;                                              /* ¸aža· ˆ•® Äa¶…Ëa   */
    *(px+n)=num;                                      /* ¤µiµA ¸aža ”·³     */
  }

  sort(1,n);

  for (i=1;i<=n;i++) {                                /* ‰i‰ÁÃ¡Ÿi ÑÁ¡eµA Â‰b */
    if (i%5==1) printf("\n");
    printf("%f  ",*(px+i));
  }
}

void sort(int left,int right) /* Quick sort */
{
  int i=left,j=right;
  double temp,mid=x[(left+right)/2];

  for (;;) {
    while (*(px+i) < mid && i<right)
      i++;
    while (mid < *(px+j) && j>left)
      j--;
    if (i<=j) {
      temp=*(px+i); /* temp“e i¤å¼µÁ j¤å¼ ¶¥­¡Ÿi Ã¡ÑÅÐa‹¡ ¶áÐe ·±¯¡¸á¸w ¥e® */
      *(px+i)=*(px+j);
      *(px+j)=temp;
      i++;
      j--;
    }
    else break;
  }

  if (left<j)      /* jŸi ‹¡º…·a¡ ¹ÁÃb·i Æâ­¡Ë· */
    sort(left,j);
  if (i<right)     /* iŸi ‹¡º…·a¡ ¶Ãb·i Æâ­¡Ë· */
    sort(i,right);
}


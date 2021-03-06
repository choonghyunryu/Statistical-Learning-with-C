/* QSORT.C */

#include <stdio.h>
#include <stdlib.h>

#define MAX 500  /* �១�i �a�a�� �A�� ���� */
#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))

void sort(int left,int right);

double *px,x[MAX];  /* �a�a�i ��w�i ���i */

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
  stream=fopen(argv[1],"rt");                         /* �a�� ���e */
  if (stream==NULL) printerrmsg("File not found!!!"); /* �a�a�a���� �i�� ���b */
                                                      /* �A���i ��            */ 
  while (fscanf(stream,"%lf\n",&num)!=EOF) {          /* �a���� �{�a�� �e��   */
    ++n;                                              /* �a�a�� ���� �a���a   */
    *(px+n)=num;                                      /* ���i�A �a�a ����     */
  }

  sort(1,n);

  for (i=1;i<=n;i++) {                                /* �i��á�i ���e�A �b */
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
      temp=*(px+i); /* temp�e i�弁�� j�弁 �����i á���a�� ���e ������w �e�� */
      *(px+i)=*(px+j);
      *(px+j)=temp;
      i++;
      j--;
    }
    else break;
  }

  if (left<j)      /* j�i �����a�� ���b�i �⭡˷ */
    sort(left,j);
  if (i<right)     /* i�i �����a�� ���b�i �⭡˷ */
    sort(i,right);
}


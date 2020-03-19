/* Program : SINUM.C ver 0.1     */
/* Author  : Ryu choong hyun     */
/* Date    : 94.8.9.	         */
/* Note    : Simple Index Number */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define FACTOR 50

void get_simple_index(double basic_term);

int row;
double data[FACTOR],time[FACTOR],simple_index[FACTOR];

void main(int argc,char *argv[])
{
  FILE *stream;
  int i,j;
  double num,basic;

  clrscr();

  if (argc<=1) {
    puts("Usage : sinum datafile");
    exit(EXIT_FAILURE);
  }

  stream=fopen(argv[1],"rt");
  if (stream==NULL) printerrmsg("File not found !!");

  fscanf(stream,"%d\n",&row);

  for (i=0;i<row;i++)
    for (j=0;j<2;j++) {
      fscanf(stream,"%lf\n",&num);
      if (j==0) time[i]=num;
      else data[i]=num;
    }

  printf("\t\t *** Simple Index Number ***\n\n");
  printf("\t\t       Basic Term = ");
  scanf("%lf",&basic);

  get_simple_index(basic);
  printf("\n\t  Time       Data     Simple-Index\n");
  for (i=0;i<row;i++)
    printf("\t%7.f   %10.5f   %7.2lf\n",time[i],data[i],simple_index[i]);

  getch();
}

void get_simple_index(double basic_term)
{
  int i;
  double p0;

  for (i=0;i<row;i++)
    if (time[i]==basic_term) {
      p0=data[i];
      break;
    }

  for (i=0;i<row;i++)
    simple_index[i]=(data[i]/p0)*100;
}


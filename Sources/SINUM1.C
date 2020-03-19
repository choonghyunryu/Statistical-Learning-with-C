/* Program : SINUM1.C ver 0.1                      */
/* Author  : Ryu choong hyun                       */
/* Date    : 94.8.10.	                           */
/* Note    : Weighted Aggregate Price Index Number */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define FACTOR 5
#define GOODS 5+1

void get_laspeyres(double basic_term);
void get_paasche(double basic_term);
void get_fisher(double basic_term);
void get_me(double basic_term);

int row,col,goods;
double p[FACTOR][GOODS],q[FACTOR][GOODS],time[FACTOR],index[FACTOR];

void main(int argc,char *argv[])
{
  FILE *stream;
  char opt;
  int i,j;
  double num,basic;

  clrscr();

  if (argc<=1) {
    puts("Usage : sinum1 datafile");
    exit(EXIT_FAILURE);
  }

  stream=fopen(argv[1],"rt");
  if (stream==NULL) printerrmsg("File not found !!");

  fscanf(stream,"%d\n",&row);
  fscanf(stream,"%d\n",&col);
  fscanf(stream,"%d\n",&goods);

  for (i=0;i<row;i++)
    for (j=0;j<col;j++) {
      fscanf(stream,"%lf\n",&num);
      if (j==0) time[i]=num;
      else if (j<=goods) p[i][j-1]=num;
      else q[i][j-goods-1]=num;
    }


  printf("\t     *** Weighted Aggregate Price Index Number ***\n\n");
  printf("\t 1.Laspeyres  2.Paasche  3.Fisher  4.Marshall-Edgeworth\n\n");

  do {
  opt=getch();
  } while (!(opt=='1' || opt=='2' || opt=='3' || opt=='4'));

  printf("\t\t\t   Basic Term = ");
  scanf("%lf",&basic);

  switch (opt) {
    case '1' : get_laspeyres(basic);
	       printf("\n\t\t\t   Time  Laspeyres\n");
	       break;
    case '2' : get_paasche(basic);
	       printf("\n\t\t\t   Time  Paasche\n");
	       break;
    case '3' : get_fisher(basic);
	       printf("\n\t\t\t   Time  Fisher\n");
	       break;
    case '4' : get_fisher(basic);
	       printf("\n\t\t\t   Time    M-E\n");
	       break;
    default  : break;
  }

  for (i=0;i<row;i++) printf("\t\t\t%7.f  %7.2lf\n",time[i],index[i]);

  getch();
}

void get_laspeyres(double basic_term)
{
  int i,j,no;
  double sum_pqn,sum_pq0=0;

  for (i=0;i<row;i++)
    if (time[i]==basic_term) {
      no=i;
      break;
    }

  for (i=0;i<goods;i++) sum_pq0+=p[no][i]*q[no][i];

  for (i=0;i<row;i++) {
    sum_pqn=0;
    for (j=0;j<goods;j++) sum_pqn+=p[i][j]*q[no][j];
    index[i]=(sum_pqn/sum_pq0)*100;
  }
}

void get_paasche(double basic_term)
{
  int i,j,no;
  double sum_pqn,sum_pq0;

  for (i=0;i<row;i++)
    if (time[i]==basic_term) {
      no=i;
      break;
    }

  for (i=0;i<row;i++) {
    sum_pqn=0;
    sum_pq0=0;
    for (j=0;j<goods;j++) {
      sum_pqn+=p[i][j]*q[i][j];
      sum_pq0+=p[no][j]*q[i][j];
    }
    index[i]=(sum_pqn/sum_pq0)*100;
  }
}

void get_fisher(double basic_term)
{
  int i,j,no;
  double sum_p0q0=0,sum_pnq0,sum_p0qn,sum_pnqn;

  for (i=0;i<row;i++)
    if (time[i]==basic_term) {
      no=i;
      break;
    }

  for (i=0;i<goods;i++) sum_p0q0+=p[no][i]*q[no][i];

  for (i=0;i<row;i++) {
    sum_p0qn=0;
    sum_pnq0=0;
    sum_pnqn=0;
    for (j=0;j<goods;j++) {
      sum_p0qn+=p[no][j]*q[i][j];
      sum_pnq0+=p[i][j]*q[no][j];
      sum_pnqn+=p[i][j]*q[i][j];
    }
    index[i]=sqrt((sum_pnq0/sum_p0q0)*(sum_pnqn/sum_p0qn))*100;
  }
}

void get_me(double basic_term)
{
  int i,j,no;
  double sum_pnq,sum_p0q;

  for (i=0;i<row;i++)
    if (time[i]==basic_term) {
      no=i;
      break;
    }

  for (i=0;i<row;i++) {
    sum_p0q=0;
    sum_pnq=0;
    for (j=0;j<goods;j++) {
      sum_p0q+=p[no][j]*(q[no][j]+q[i][j]);
      sum_pnq+=p[i][j]*(q[no][j]+q[i][j]);
    }
    index[i]=(sum_pnq/sum_p0q)*100;
  }
}


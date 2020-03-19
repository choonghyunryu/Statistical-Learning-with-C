/* Program : CAURND.C ver 0.1     */
/* Author  : Ryu choong hyun      */
/* Date    : 94.8.21.	  	  */
/* Note    : Cauchy Random Number */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <time.h>

#define PI M_PI
#define N 1000

double rnd(void);
double caurnd(void);

int no;
double x[N];

void main(void)
{
  int i;
  double xx,sum=0,sum_s=0,mean,var,dev;

  clrscr();

  randomize();

  printf("\t\t** Cauchy Random Numbers **\n\n");
  printf("\t\tNumber of Random-Number = ");
  scanf("%d",&no);

  for (i=0;i<no;i++) {
    x[i]=caurnd();
    sum+=x[i];
  }

  mean=sum/no;

  for (i=0;i<no;i++) {
    sum_s+=pow(x[i]-mean,2);
  }

  var=sum_s/no;
  dev=sqrt(var);

  printf("\n\t\t      < Random-Number >\n");
  for (i=0;i<no;i++)  printf("\t\t\t%10.6f\n",x[i]);
  printf("\n\t\t    Mean  = %8.4f\n",mean);
  printf("\t\t    Var   = %8.4f\n",var);
  printf("\t\t    S.D.  = %8.4f\n\n",dev);

  getch();
}

double rnd(void)
{
  return (1.0/(RAND_MAX+1.0))*rand();
}

double caurnd(void)
{
  return tan(PI*rnd());
}


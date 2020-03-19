/* Program : EXPRND.C ver 0.1          */
/* Author  : Ryu choong hyun           */
/* Date    : 94.8.21.	  	       */
/* Note    : Exponential Random Number */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <time.h>

#define N 1000

double rnd(void);
double exprnd(double lambda);

int no;
double x[N];

void main(void)
{
  int i;
  double lambda,sum=0,sum_s=0,mu,sigma2,sigma,mean,var,dev;

  clrscr();

  randomize();

  printf("\t\t** Exponential Random Numbers **\n\n");
  printf("\t\tNumber of Random-Number = ");
  scanf("%d",&no);
  printf("\t\tExp(lambda) => lambda   = ");
  scanf("%lf",&lambda);

  for (i=0;i<no;i++) {
    x[i]=exprnd(lambda);
    sum+=x[i];
  }

  mu=1./lambda;
  mean=sum/no;

  for (i=0;i<no;i++) {
    sum_s+=pow(x[i]-mean,2);
  }

  sigma2=1./pow(lambda,2);
  sigma=sqrt(sigma2);
  var=sum_s/no;
  dev=sqrt(var);

  printf("\n\t\t      < Random-Number >\n");
  for (i=0;i<no;i++)  printf("\t\t\t%10.6f\n",x[i]);
  printf("\n\t\t Mean  = %lf(%lf)\n",mean,mu);
  printf("\t\t Var   = %lf(%lf)\n",var,sigma2);
  printf("\t\t S.D.  = %lf(%lf)\n\n",dev,sigma);

  getch();
}

double rnd(void)
{
  return (1.0/(RAND_MAX+1.0))*rand();
}

double exprnd(double lambda)
{
  return -1/lambda*log(1-rnd());
}


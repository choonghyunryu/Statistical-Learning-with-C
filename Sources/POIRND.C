/* Program : POIRND.C ver 0.1      */
/* Author  : Ryu choong hyun       */
/* Date    : 94.8.21.	  	   */
/* Note    : Poisson Random Number */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <time.h>

#define N 1000

double rnd(void);
int poirnd(double lambda);
void probability(void);

int no,x[N],max,min=1000;
double pro[N];

void main(void)
{
  int i,n;
  double lambda,sum=0,sum_s=0,mean,var,dev,r;

  clrscr();

  randomize();

  printf("\t\t** Poisson Random Numbers **\n\n");
  printf("\t\tNumber of Random-Number   = ");
  scanf("%d",&no);
  printf("\t\tPoisson(lambda) => lambda = ");
  scanf("%lf",&lambda);

  for (i=0;i<no;i++) {
    x[i]=poirnd(lambda);
    sum+=x[i];
  }

  mean=sum/no;

  for (i=0;i<no;i++) {
    sum_s+=pow(x[i]-mean,2);
    if (min>x[i]) min=x[i];
    if (max<x[i]) max=x[i];
  }

  var=sum_s/no;
  dev=sqrt(var);

  printf("\n\t\t      < Random-Number >\n");
  for (i=0;i<no;i++) {
    if (i%10==0) printf("\t %4d",x[i]);
    else if (i%10==9) printf("%4d\n",x[i]);
    else printf("%4d",x[i]);
  }
  printf("\n\t\t Mean  = %lf(%lf)\n",mean,lambda);
  printf("\t\t Var   = %lf(%lf)\n",var,lambda);
  printf("\t\t S.D.  = %lf(%lf)\n\n",dev,sqrt(lambda));

  probability();
  for (i=min;i<=max;i++) printf("\t\t\t%2d   %5.3f\n",i,pro[i]);

  getch();
}

double rnd(void)
{
  return (1.0/(RAND_MAX+1.0))*rand();
}

int poirnd(double lambda)
{
  int i,sum=0;
  double l,p=1;

  l=exp(-lambda);

  for (;;) {
    p*=rnd();
    if (p>=l) sum++;
    else return sum;
  }
}

void probability(void)
{
  int i,j,temp[N]={0,};

  for (i=0;i<=no;i++) {
    j=x[i];
    temp[j]++;
  }

  for (i=0;i<=no;i++) pro[i]=(double)temp[i]/no;
}

/* Program : BINRND.C ver 0.1       */
/* Author  : Ryu choong hyun        */
/* Date    : 94.8.21.	  	    */
/* Note    : Binomial Random Number */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <time.h>

#define N 1000

double rnd(void);
int binrnd(int n,double p);
void probability(void);

int no,x[N];
double pro[N];

void main(void)
{
  int i,n;
  double p,sum=0,sum_s=0,mu,sigma2,mean,var,dev,r;

  clrscr();

  randomize();

  printf("\t\t** Binomial Random Numbers **\n\n");
  printf("\t\tNumber of Random-Number = ");
  scanf("%d",&no);
  printf("\t\tBinomial(n,p) => n      = ");
  scanf("%d",&n);
  printf("\t\tBinomial(n,p) => p      = ");
  scanf("%lf",&p);

  for (i=0;i<no;i++) {
    x[i]=binrnd(n,p);
    sum+=x[i];
  }

  mean=sum/no;

  for (i=0;i<no;i++) {
    sum_s+=pow(x[i]-mean,2);
  }

  mu=n*p;
  sigma2=n*p*(1-p);

  var=sum_s/no;
  dev=sqrt(var);

  printf("\n\t\t      < Random-Number >\n");
  for (i=0;i<no;i++) {
    if (i%10==0) printf("\t %4d",x[i]);
    else if (i%10==9) printf("%4d\n",x[i]);
    else printf("%4d",x[i]);
  }
  printf("\n\t\t Mean  = %lf(%lf)\n",mean,mu);
  printf("\t\t Var   = %lf(%lf)\n",var,sigma2);
  printf("\t\t S.D.  = %lf(%lf)\n\n",dev,sqrt(sigma2));

  probability();
  for (i=0;i<=n;i++) printf("\t\t\t%2d   %5.3f\n",i,pro[i]);

  getch();
}

double rnd(void)
{
  return (1.0/(RAND_MAX+1.0))*rand();
}

int binrnd(int n,double p)
{
  int i,sum=0;

  for (i=0;i<n;i++) if (rnd()<p) sum++;

  return sum;
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

/* Program : UNIRND.C ver 0.1      */
/* Author  : Ryu choong hyun       */
/* Date    : 94.8.20.		   */
/* Note    : Uniform Random Number */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <time.h>

#define N 1000

double rnd(void);
double get_r(void);
void probability(void);

int n;
double x[N],pro[10];

void main(void)
{
  int i;
  double sum=0,sum_s=0,mean,dev,r;

  clrscr();

  randomize();

  printf("\t\t** Uniform Random Numbers (0 to 1) **\n\n");
  printf("\t\t   Number of Random-Number = ");
  scanf("%d",&n);

  for (i=0;i<n;i++) {
    x[i]=rnd();
    sum+=x[i];
  }

  mean=sum/n;

  for (i=0;i<n;i++) {
    sum_s+=pow(x[i]-mean,2);
  }

  dev=sqrt(sum_s/n);
  r=get_r();

  printf("\n\t\t\t  Random-Number\n");
  for (i=0;i<n;i++) printf("\t\t\t    %lf\n",x[i]);
  printf("\n\t\t\t Mean  = %lf\n",mean);
  printf("\t\t      (%lf , %lf)\n",.5-2/sqrt(12*n),.5+2/sqrt(12*n));
  printf("\t\t\t S.D.  = %lf\n",dev);
  printf("\t\t\t 1st R = %lf\n",r);
  printf("\t\t      (%lf , %lf)\n\n",-1./(n-1)-2*sqrt(n*(n-3)/(n+1))/(n-1),
				     -1./(n-1)+2*sqrt(n*(n-3)/(n+1))/(n-1));
  probability();
  printf("\t\t\t Interval     P\n");
  for (i=0;i<10;i++) printf("\t\t\t[%3.1f ,%3.1f)  %5.3f\n",
			     i/10.,(i+1.)/10.,pro[i]);

  getch();
}

double rnd(void)
{
  return (1.0/(RAND_MAX+1.0))*(rand()+0.5);
}

double get_r(void)
{
  int i;
  double sum=0,sum_mult=0,sum_s=0;

  for (i=0;i<n;i++) {
    sum+=x[i];
    sum_s+=pow(x[i],2);
    if (i!=n-1) sum_mult+=x[i]*x[i+1];
  }
  sum_mult+=x[0]*x[n-1];

  return (n*sum_mult-pow(sum,2))/(n*sum_s-pow(sum,2));
}

void probability(void)
{
  int i,j,temp[10]={0,};

  for (i=0;i<n;i++) {
    j=x[i]*10;
    temp[j]++;
  }

  for (i=0;i<10;i++) pro[i]=(double)temp[i]/n;
}

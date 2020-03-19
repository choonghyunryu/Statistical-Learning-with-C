/* Program : BETARND.C ver 0.1   */
/* Author  : Ryu choong hyun     */
/* Date    : 94.8.21.	  	 */
/* Note    : Beta Random Number  */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <time.h>

#define N 1000

double rnd(void);
double betarnd(double alpha,double beta);

int no;
double x[N];

void main(void)
{
  int i;
  double alpha,beta,sum=0,sum_s=0,mu,sigma2,sigma,mean,var,dev;

  clrscr();

  randomize();

  printf("\t\t** Beta Random Numbers **\n\n");
  printf("\t\tNumber of Random-Number = ");
  scanf("%d",&no);
  printf("\t\tBeta(\xe0,\xe1) => \xe0          = ");
  scanf("%lf",&alpha);
  printf("\t\tBeta(\xe0,\xe1) => \xe1          = ");
  scanf("%lf",&beta);

  for (i=0;i<no;i++) {
    x[i]=betarnd(alpha,beta);
    sum+=x[i];
  }

  mu=alpha/(alpha+beta);
  mean=sum/no;

  for (i=0;i<no;i++) {
    sum_s+=pow(x[i]-mean,2);
  }

  sigma2=alpha*beta/(pow(alpha+beta,2)*(alpha+beta+1));
  sigma=sqrt(sigma2);
  var=sum_s/no;
  dev=sqrt(var);

  printf("\n\t\t    < Random-Number >\n");
  for (i=0;i<no;i++)  printf("\t\t\t%lf\n",x[i]);
  printf("\n\t\t Mean  = %lf(%lf)\n",mean,mu);
  printf("\t\t Var   = %lf(%lf)\n",var,sigma2);
  printf("\t\t S.D.  = %lf(%lf)\n\n",dev,sigma);

  getch();
}

double rnd(void)
{
  return (1.0/(RAND_MAX+1.0))*rand();
}

double betarnd(double alpha,double beta)
{
  double u1,u2,v1,v2,w;

  do {
    u1=rnd();
    u2=rnd();
    v1=pow(u1,1/alpha);
    v2=pow(u2,1/beta);
    w=v1+v2;
    if (w<=1) return v1/w;
  } while (1);
}

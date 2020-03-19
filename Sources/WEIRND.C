/* Program : WEIPRND.C ver 0.1     */
/* Author  : Ryu choong hyun       */
/* Date    : 94.8.21.	  	   */
/* Note    : Weibull Random Number */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <time.h>

#define N 1000

double rnd(void);
double weirnd(double lambda,double alpha);
double gammaf(double x);

int no;
double x[N];

void main(void)
{
  int i;
  double lambda,alpha,sum=0,sum_s=0,mu,sigma2,sigma,mean,var,dev;

  clrscr();

  randomize();

  printf("\t\t** Weibull Random Numbers **\n\n");
  printf("\t\tNumber of Random-Number = ");
  scanf("%d",&no);
  printf("\t\tExp(lambda,\xe0) => lambda = ");
  scanf("%lf",&lambda);
  printf("\t\tExp(lambda,\xe0) => \xe0      = ");
  scanf("%lf",&alpha);

  for (i=0;i<no;i++) {
    x[i]=weirnd(lambda,alpha);
    sum+=x[i];
  }

  mu=pow(lambda,-1/alpha)*gammaf(1/alpha+1);
  mean=sum/no;

  for (i=0;i<no;i++) {
    sum_s+=pow(x[i]-mean,2);
  }

  sigma2=pow(lambda,-2/alpha)*(gammaf(2/alpha+1)-pow(gammaf(1/alpha+1),2));
  sigma=sqrt(sigma2);
  var=sum_s/no;
  dev=sqrt(var);

  printf("\n\t\t      < Random-Number >\n");
  for (i=0;i<no;i++)  printf("\t\t\t%10.4f\n",x[i]);
  printf("\n\t\t Mean  = %10.2f(%10.2f)\n",mean,mu);
  printf("\t\t Var   = %10.2f(%10.2f)\n",var,sigma2);
  printf("\t\t S.D.  = %10.2f(%10.2f)\n\n",dev,sigma);

  getch();
}

double rnd(void)
{
  return (1.0/(RAND_MAX+1.0))*rand();
}

double weirnd(double lambda,double alpha)
{
  return pow(-log(1-rnd())/lambda,1/alpha);
}

double gammaf(double x)
{
  double a,b,d,p,q;

  b=x;
  a=1/(x*(x+1));

  while(1) {
    b-=1;
    if (b<0) break;
    a*=b+2;
  }

  d=(int)x;
  b=x-d;
  p=(((-.4530104765151624390265e2*b-.26102849506164577995484e3)*b-
       .18177991911564340957213e4)*b-.495476897267726670260117e4)*b-
       .1357265011103685653360092e5;
  q=((((((b-.1616160629773462081477e2)*b+.9377691732243930278552e2)*b-
	    .1247947256227383271411e3)*b-.93117850455934846630217e3)*b+
	    .344069924134645411261773e4)*b+.78353488005595515120714e3)*b-
	    .1357265011103690542985686e5;
  return a*p/q;
}

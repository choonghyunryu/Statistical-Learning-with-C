/* Program : FRND.C ver 0.1   */
/* Author  : Ryu choong hyun  */
/* Date    : 94.8.22.	      */
/* Note    : F Random Number  */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <time.h>

#define E M_E
#define N 1000

double rnd(void);
double gammarnd(double alpha);
double exprnd(double lambda);
double chirnd(double df);
double frnd(double u,double v);

int no;
double x[N];

void main(void)
{
  int i;
  double df1,df2,sum=0,sum_s=0,mean,var,dev,mu,sigma2,sigma;

  clrscr();

  randomize();

  printf("\t\t** F Random Numbers **\n\n");
  printf("\t\tNumber of Random-Number = ");
  scanf("%d",&no);
  printf("\t\tF(u,v) => u             = ");
  scanf("%lf",&df1);
  printf("\t\tF(u,v) => v             = ");
  scanf("%lf",&df2);

  for (i=0;i<no;i++) {
    x[i]=frnd(df1,df2);
    sum+=x[i];
  }

  if (df2!=2) mu=df2/(df2-2);
  mean=sum/no;

  for (i=0;i<no;i++) {
    sum_s+=pow(x[i]-mean,2);
  }

  if (df2>4) {
    sigma2=2*pow(df2,2)*(df1+df2-2)/(df1*pow(df2-2,2)*(df2-4));
    sigma=sqrt(sigma2);
  }
  var=sum_s/no;
  dev=sqrt(var);

  printf("\n\t\t    < Random-Number >\n");
  for (i=0;i<no;i++)  printf("\t\t\t%9.5f\n",x[i]);
  if (df2!=2) printf("\n\t\t Mean  = %9.5f(%9.5f)\n",mean,mu);
  else printf("\n\t\t Mean  = %9.5f(*******)\n",mean);
  if (df2>4) {
    printf("\t\t Var   = %9.5f(%9.5f)\n",var,sigma2);
    printf("\t\t S.D.  = %9.5f(%9.5f)\n\n",dev,sigma);
  }
  else {
    printf("\t\t Var   = %9.5f(*******)\n",var);
    printf("\t\t S.D.  = %9.5f(*******)\n\n",dev);
  }

  getch();
}

double rnd(void)
{
  return (1.0/(RAND_MAX+1.0))*rand();
}

double gammarnd(double alpha)
{
  double u,b,p,x,y,temp;

  if (alpha<1) {
    b=(E+alpha)/E;
    p=b*rnd();
    do {
      if (p>1) {
	x=-log((b-p)/alpha);
	temp=pow(x,alpha-1);
      }
      else {
	x=pow(p,1/alpha);
	temp=exp(-x);
      }
    } while (rnd()>=temp);
  }
  else if (alpha==1) x=exprnd(1);
  else {
    temp=sqrt(2*alpha-1);
    do {
      do {
	do {
	  x=1-rnd();
	  y=2*rnd()-1;
	} while (pow(x,2)+pow(y,2)>1);
	y/=x;
	x=temp*y+alpha-1;
      } while (x<=0);
      u=(alpha-1)*log(x/(alpha-1))-temp*y;
    } while (u<-50 || rnd()>(1+pow(y,2))*exp(u));
  }

  return x;
}

double exprnd(double lambda)
{
  return -1/lambda*log(1-rnd());
}

double chirnd(double df)
{
  return 2*gammarnd(df/2);
}

double frnd(double u,double v)
{
  return (chirnd(u)*v)/(chirnd(v)*u);
}
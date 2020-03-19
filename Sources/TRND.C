/* Program : TRND.C ver 0.1    */
/* Author  : Ryu choong hyun   */
/* Date    : 94.8.22.	       */
/* Note    : T Random Number   */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <time.h>

#define E M_E
#define PI M_PI
#define N 1000

double rnd(void);
double exprnd(double lambda);
double gammarnd(double alpha);
double chirnd(double df);
double nornd2(void);
double trnd(double df);

int no;
double x[N];

void main(void)
{
  int i;
  double df,sum=0,sum_s=0,mean,var,dev,sigma2,sigma;

  clrscr();

  randomize();

  printf("\t\t    ** T Random Numbers **\n\n");
  printf("\t\tNumber of Random-Number = ");
  scanf("%d",&no);
  printf("\t\tX\xfd(df) => df            = ");
  scanf("%lf",&df);

  for (i=0;i<no;i++) {
    x[i]=trnd(df);
    sum+=x[i];
  }

  mean=sum/no;

  for (i=0;i<no;i++) {
    sum_s+=pow(x[i]-mean,2);
  }

  if (df!=2) {
    sigma2=df/(df-2);
    sigma=sqrt(sigma2);
  }
  var=sum_s/no;
  dev=sqrt(var);

  printf("\n\t\t    < Random-Number >\n");
  for (i=0;i<no;i++)  printf("\t\t\t%8.4f\n",x[i]);
  printf("\n\t\t Mean  = %8.4f(%8.4f)\n",mean,0.);
  if (df!=2) {
    printf("\t\t Var   = %8.4f(%8.4f)\n",var,sigma2);
    printf("\t\t S.D.  = %8.4f(%8.4f)\n\n",dev,sigma);
  }
  else {
    printf("\t\t Var   = %8.4f(*****)\n",var);
    printf("\t\t S.D.  = %8.4f(*****)\n\n",dev);
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

double nornd2(void)
{
  static int i=0;
  static double u,temp;

  if (i==0) {
    i++;
    temp=sqrt(-2*log(rnd()));
    u=2*PI*rnd();
    return temp*cos(u);
  }
  else {
    i=0;
    return temp*sin(u);
  }
}

double trnd(double df)
{
  return nornd2()/sqrt(chirnd(df)/df);
}
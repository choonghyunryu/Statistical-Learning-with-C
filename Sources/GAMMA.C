/* Program : GAMMA.C ver 0.1                  */
/* Author  : Ryu choong hyun                  */
/* Date    : 94.7.9.		      	      */
/* Note    : Gamma Probability Distribution   */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>

#define PI M_PI

double nor_a(double z);
double gammaf(double x);
double gamma_p(double alpha,double t);
double logamma(double y);

double alpha,beta;

void main(void)
{
  int i;
  char c;
  double a,b,mean,var,dev,value;

  clrscr();
  printf("\n\t ** Gamma Probability Distribution **\n");
  printf("\n\t             Alpha     = ");
  scanf("%lf",&alpha);
  printf("\n\t             Beta      = ");
  scanf("%lf",&beta);
  mean=alpha*beta;
  var=mean*beta;
  dev=sqrt(var);

  printf("\n\t               << Interval >>\n");
  printf("\t        1. a<x<b   2. x>a   3. x<a\n");
  do
    switch (c=getch()) {
      case '1' : printf("\n\t             a    = ");
		 scanf("%lf",&a);
		 printf("\n\t             b    = ");
		 scanf("%lf",&b);
		 value=gamma_p(alpha,b*beta)-
		 gamma_p(alpha,a*beta);
		 break;
      case '2' : printf("\n\t             a    = ");
		 scanf("%lf",&a);
		 value=1-gamma_p(alpha,a*beta);
		 break;
      case '3' : printf("\n\t             a    = ");
		 scanf("%lf",&a);
		 value=gamma_p(alpha,a*beta);
		 break;
      default  : break;
    } while (!(c=='1' || c=='2' || c=='3'));

  printf("\n\n\t\t Mean        = %lf\n\n",mean);
  printf("\t\t Devation    = %lf\n\n",dev);
  printf("\t\t Variance    = %lf\n\n",var);
  printf("\t\t Probability = %lf\n\n",value);
  getch();
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

double gamma_p(double alpha,double t)
{
  double oflo,tol,elimit,arg,gamma,a,b,c,an,rn,pn1,pn2,pn3,pn4,pn5,pn6;

  oflo=1.0e30;
  tol=1.0e-7;
  elimit=-88;

  if (t==0) return 0;

  if (alpha>1000)
    return nor_a(3*sqrt(alpha)*(pow(t/alpha,1/3.)+1/(9*alpha)-1));

  if (t>1.0e6) return 1;

  if (t<=1 || t<alpha) {
    arg=alpha*log(t)-t-logamma(alpha+1);
    c=1;
    gamma=1;
    a=alpha;

    do {
      a+=1;
      c*=(t/a);
      gamma+=c;
    } while (c>tol);

    arg+=log(gamma);

    if (arg>=elimit) return exp(arg);
  }
  else {
    arg=alpha*log(t)-t-logamma(alpha);
    a=1-alpha;
    b=a+t+1;
    c=0;
    pn1=1;
    pn2=t;
    pn3=t+1;
    pn4=t*b;
    gamma=pn3/pn4;

    while (1) {
      a+=1;
      b+=2;
      c+=1;
      an=a*c;
      pn5=b*pn3-an*pn1;
      pn6=b*pn4-an*pn2;

      if (fabs(pn6)>0) {
	rn=pn5/pn6;
	if (fabs(gamma-rn)<=min(tol,tol*rn)) {
	  arg+=log(gamma);
	  if (arg>=elimit) return 1-exp(arg);
	}
	gamma=rn;
      }

      pn1=pn3;
      pn2=pn4;
      pn3=pn5;
      pn4=pn6;

      if (fabs(pn5)>=oflo) {
	pn1/=oflo;
	pn2/=oflo;
	pn3/=oflo;
	pn4/=oflo;
      }
    }
  }
  return gamma;
}

double logamma(double y)
{
  int i;
  double x,a,b,c,n,lgamma;

  if (y-8>=0) {
    x=y;
    n=1;
  }
  else {
    x=y+8;
    n=-1;
  }

  c=1/pow(x,2);
  a=(x-.5)*log(x)-x+.918938533204673;
  b=((((.000766345188*c-.00059409561052)*c+.0007936431104845)*c-
	.00277777775657725)*c+.0833333333333169234)/x;

  lgamma=a+b;
  if (n>=0) return lgamma;

  x-=1;
  a=x;

  for (i=1;i<=7;i++) a*=x-i;

  return lgamma-log(a);
}

double nor_a(double z)
{
  double d[6],p;

  d[0]=.0498673470;
  d[1]=.0211410061;
  d[2]=.0032776263;
  d[3]=.0000380036;
  d[4]=.0000488906;
  d[5]=.0000053830;

  p=1-pow(1+d[0]*z+d[1]*pow(z,2)+d[2]*pow(z,3)+d[3]*pow(z,4)+d[4]*
	  pow(z,5)+d[5]*pow(z,6),-16)/2;
  return 2*p-1;
}



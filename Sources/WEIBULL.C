/* WEIBULL.C  Ver 0.1        */
/* Author : Ryu Choong Hyun  */
/* Date   : 94.7.11.         */
/* Weibull Probability &     */
/* Percentage Point          */

#include <stdio.h>
#include <conio.h>
#include <math.h>

double weibull_p(double lambda,double alpha,double x);
double weibull_point(double lambda,double alpha,double point);
double gammaf(double x);

void main(void)
{
  char opt,c;
  double a,b,lambda,alpha,value,point,mean,var;

  clrscr();
  printf("\n\t    ** Weibull  Distribution **\n");
  printf("\n\t1. Probability  2. Percentage Point\n");

  do
  switch (opt=getch()) {
    case '1' :
      printf("\n\t1. P(a<X<b)  2. P(X>a)  3. P(X<a)\n");
      do
      switch (c=getch()) {
      case '1' : printf("\n\t   Lambda = ");
		 scanf("%lf",&lambda);
		 printf("\n\t   Alpha  = ");
		 scanf("%lf",&alpha);
		 printf("\n\t   From   = ");
		 scanf("%lf",&a);
		 printf("\n\t   To     = ");
		 scanf("%lf",&b);
		 value=weibull_p(lambda,alpha,b)-weibull_p(lambda,alpha,a);
		 printf("\n\n\t  ** Exponential  Distribution **\n\n");
		 printf("\t    P(%f \x3c X \x3c %f) = %f\n\n",a,b,value);
		 break;
      case '2' : printf("\n\t   Lambda = ");
		 scanf("%lf",&lambda);
		 printf("\n\t    Alpha = ");
		 scanf("%lf",&alpha);
		 printf("\n\t    a    = ");
		 scanf("%lf",&a);
		 value=1-weibull_p(lambda,alpha,a);
		 printf("\n\n\t  ** Weibull  Distribution **\n\n");
		 printf("\t    P(X \x3e %f) = %lf\n\n",a,value);
		 break;
      case '3' : printf("\n\t   Lambda = ");
		 scanf("%lf",&lambda);
		 printf("\n\t    Alpha = ");
		 scanf("%lf",&alpha);
		 printf("\n\t    a    = ");
		 scanf("%lf",&a);
		 value=weibull_p(lambda,alpha,a);
		 printf("\n\n\t  ** Weibull  Distribution **\n\n");
		 printf("\t    P(X \x3c %f) = %lf\n\n",a,value);
		 break;
      default :  break;
    } while (!(c=='1' || c=='2' || c=='3'));
    break;
    case '2' : printf("\n\t1. Lower Percentage Point"
		      "  2. Upper Percentage Point\n\n");
	       do
	       switch (opt=getch()) {
		 case '1' :
		 case '2' : printf("\n\t   Lambda = ");
			    scanf("%lf",&lambda);
			    printf("\n\t    Alpha = ");
			    scanf("%lf",&alpha);
			    printf("\n\tPercentage Point : ");
			    scanf("%lf",&point);
			    break;
		 default  : break;
	       } while (!(opt=='1' || opt=='2'));

	       printf("\n\n\t** Weibull  Distribution **\n\n");
	       if (opt=='1') {
		 a=weibull_point(lambda,alpha,point);
		 printf("\t %5.2f Percentage = %f\n\n",point,a);
	       }
	       else {
		 a=weibull_point(lambda,alpha,1-point);
		 printf("\t %5.2f Percentage = %f\n\n",point,a);
	       }
	       break;
    default  : break;
  } while (!(opt=='1' || opt=='2'));

  mean=pow(lambda,-1/alpha)*gammaf(1/alpha+1);
  var=pow(lambda,-2/alpha)*(gammaf(2/alpha+1)-pow(gammaf(1/alpha+1),2));
  printf("\t Mean = %f Variance = %f",mean,var);
  getch();
}

double weibull_p(double lambda,double alpha,double x)
{
  return 1-exp(-lambda*pow(x,alpha));
}

double weibull_point(double lambda,double alpha,double point)
{
  return pow(-log(1-point)/lambda,1/alpha);
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

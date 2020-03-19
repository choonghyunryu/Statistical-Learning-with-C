/* Program : UNIFORM.C ver 0.1    */
/* Author  : Ryu choong hyun      */
/* Date    : 94.7.9.		  */
/* Note    : Uniform Distribtion  */

#include <stdio.h>
#include <conio.h>
#include <math.h>

double uni_p(double x);
double uni_point(double point);

char job;
double a,b;

void main(void)
{
  char opt,c;
  double i,ax,bx,mean,var,aval,bval,value,point;

  clrscr();
  printf("\n\t ** Uniform Distribution **\n");
  printf("\n\t    First Interval  = ");
  scanf("%lf",&a);
  printf("\n\t    Last  Interval  = ");
  scanf("%lf",&b);
  printf("\n\t1. Probability  2. Percentage Point\n");

  mean=(a+b)/2;
  var=(b-a)*(b-a)/12;

  do
  switch (job=getch()) {
    case '1' :
      printf("\n\t1. P(a<X<b)  2. P(X>a)  3. P(X<a)\n");
      do
      switch (c=getch()) {
      case '1' : printf("\n\t    From = ");
		 scanf("%lf",&ax);
		 printf("\n\t    To   = ");
		 scanf("%lf",&bx);
		 value=(ax<a) ? uni_p(bx) : (bx>b) ? 1-uni_p(ax) :
					       uni_p(bx)-uni_p(ax);         ;
		 clrscr();
		 printf("\t      ** Uniform Probability Distribution **\n\n");
		 printf("\t\tMean = %-10.4lf  Var = %-10.4lf\n",mean,var);
		 printf("\t\tP(%f \x3c X \x3c %f) = %f\n\n",ax,bx,value);
		 break;
      case '2' : printf("\n\t    a    = ");
		 scanf("%lf",&ax);
		 bx=b;
		 value=1-uni_p(ax);
		 clrscr();
		 printf("\t      ** Uniform Probability Distribution **\n\n");
		 printf("\t\tMean = %-10.4lf   Var = %-10.4lf\n",mean,var);
		 printf("\t\tP(X \x3e %f) = %lf\n\n",ax,value);
		 break;
      case '3' : printf("\n\t    a    = ");
		 scanf("%lf",&ax);
		 bx=a;
		 value=uni_p(ax);
		 clrscr();
		 printf("\t      ** Uniform Probability Distribution **\n\n");
		 printf("\t\tMean = %-10.4lf   Var = %-10.4lf\n",mean,var);
		 printf("\t\tP(X \x3c %f) = %lf\n\n",ax,value);
		 break;
      default :  break;
    } while (!(c=='1' || c=='2' || c=='3'));
    break;
    case '2' : printf("\n\t1. Lower Percentage Point"
		      "  2. Upper Percentage Point\n\n");
	       do
	       switch (opt=getch()) {
		 case '1' :
		 case '2' : printf("\n\tPercentage Point : ");
			    scanf("%lf",&point);
			    break;
		 default  : break;
	       } while (!(opt=='1' || opt=='2'));

	       clrscr();
	       printf("\t      ** Uniform Probability Distribution **\n\n");
	       printf("\t\tMean = %-10.4lf  Var = %-10.4lf\n",mean,var);
	       if (opt=='1') {
		 ax=uni_point(point);
		 printf("\t\t%5.2f Percentage = %f\n\n",point,ax);
	       }
	       else {
		 ax=b-(uni_point(point)-a);
		 printf("\t\t%5.2f Percentage = %f\n\n",point,ax);
	       }
	       break;
    default  : break;
  } while (!(job=='1' || job=='2'));

  getch();
}

double uni_p(double x)
{
  return (b>=x && x>=a ) ? 1/(b-a)*(x-a) : (x<=a) ? 0 : 1 ;
}

double uni_point(double point)
{
  return point*(b-a)+a;
}


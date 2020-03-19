/* CAUCHY.C  Ver 0.1        */
/* Author : Ryu Choong Hyun */
/* Date   : 94.7.11.        */
/* Cauchy Probability &     */
/* Percentage Point         */

#include <stdio.h>
#include <conio.h>
#include <math.h>

#define PI M_PI

double cauchy_p(double x);
double cauchy_point(double point);

void main(void)
{
  char opt,c;
  double a,b,value,point;

  clrscr();
  printf("\n\t ** Cauchy  Distribution **\n");
  printf("\n\t1. Probability  2. Percentage Point\n");

  do
  switch (opt=getch()) {
    case '1' :
      printf("\n\t1. P(a<X<b)  2. P(X>a)  3. P(X<a)\n");
      do
      switch (c=getch()) {
      case '1' : printf("\n\t    From = ");
		 scanf("%lf",&a);
		 printf("\n\t    To   = ");
		 scanf("%lf",&b);
		 value=cauchy_p(b)-cauchy_p(a);
		 printf("\n\n\t      ** Cauchy  Distribution **\n\n");
		 printf("\t\tP(%f \x3c X \x3c %f) = %f\n\n",a,b,value);
		 break;
      case '2' : printf("\n\t    a    = ");
		 scanf("%lf",&a);
		 value=1-cauchy_p(a);
		 printf("\n\n\t      ** Cauchy  Distribution **\n\n");
		 printf("\t\tP(X \x3e %f) = %lf\n\n",a,value);
		 break;
      case '3' : printf("\n\t    a    = ");
		 scanf("%lf",&a);
		 value=cauchy_p(a);
		 printf("\n\n\t      ** Cauchy  Distribution **\n\n");
		 printf("\t\tP(X \x3c %f) = %lf\n\n",a,value);
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

	       printf("\n\n\t      ** Cauchy  Distribution **\n\n");
	       if (opt=='1') {
		 a=cauchy_point(point);
		 printf("\t\t%5.2f Percentage = %f\n\n",point,a);
	       }
	       else {
		 a=cauchy_point(1-point);
		 printf("\t\t%5.2f Percentage = %f\n\n",point,a);
	       }
	       break;
    default  : break;
  } while (!(opt=='1' || opt=='2'));

  getch();
}

double cauchy_p(double x)
{
  if (x==0) return .5;
  return (x<0) ? .5-1/PI*atan(fabs(x)) : .5+1/PI*atan(x); 
}

double cauchy_point(double point)
{
  if (point==.5) return 0; 
  return (point<.5) ? -tan(PI*(.5-point)) : tan(PI*(point-.5));
}

/* EXPDIS.C  Ver 0.1         */
/* Author : Ryu Choong Hyun  */
/* Date   : 94.7.11.         */
/* Exponential Probability & */
/* Percentage Point          */

#include <stdio.h>
#include <conio.h>
#include <math.h>

double exp_p(double lamda,double x);
double exp_point(double lambda,double point);

void main(void)
{
  char opt,c;
  double a,b,lambda,value,point;

  clrscr();
  printf("\n\t ** Exponential  Distribution **\n");
  printf("\n\t1. Probability  2. Percentage Point\n");

  do
  switch (opt=getch()) {
    case '1' :
      printf("\n\t1. P(a<X<b)  2. P(X>a)  3. P(X<a)\n");
      do
      switch (c=getch()) {
      case '1' : printf("\n\t    Lambda = ");
		 scanf("%lf",&lambda);
		 printf("\n\t    From  = ");
		 scanf("%lf",&a);
		 printf("\n\t    To    = ");
		 scanf("%lf",&b);
		 value=exp_p(lambda,b)-exp_p(lambda,a);
		 printf("\n\n\t  ** Exponential  Distribution **\n\n");
		 printf("\t    P(%f \x3c X \x3c %f) = %f\n\n",a,b,value);
		 break;
      case '2' : printf("\n\t    Lambda = ");
		 scanf("%lf",&lambda);
		 printf("\n\t    a    = ");
		 scanf("%lf",&a);
		 value=1-exp_p(lambda,a);
		 printf("\n\n\t  ** Exponential  Distribution **\n\n");
		 printf("\t    P(X \x3e %f) = %lf\n\n",a,value);
		 break;
      case '3' : printf("\n\t    Lambda = ");
		 scanf("%lf",&lambda);
		 printf("\n\t    a    = ");
		 scanf("%lf",&a);
		 value=exp_p(lambda,a);
		 printf("\n\n\t  ** Exponential  Distribution **\n\n");
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
		 case '2' : printf("\n\t    Lambda = ");
			    scanf("%lf",&lambda);
			    printf("\n\tPercentage Point : ");
			    scanf("%lf",&point);
			    break;
		 default  : break;
	       } while (!(opt=='1' || opt=='2'));

	       printf("\n\n\t** Exponential  Distribution **\n\n");
	       if (opt=='1') {
		 a=exp_point(lambda,point);
		 printf("\t %5.2f Percentage = %f\n\n",point,a);
	       }
	       else {
		 a=exp_point(lambda,1-point);
		 printf("\t %5.2f Percentage = %f\n\n",point,a);
	       }
	       break;
    default  : break;
  } while (!(opt=='1' || opt=='2'));

  printf("\t Mean = %f Variance = %f",1/lambda,1/pow(lambda,2));
  getch();
}

double exp_p(double lambda,double x)
{
  return 1-exp(-lambda*x);
}

double exp_point(double lambda,double point)
{
  return -log(1-point)/lambda;
}
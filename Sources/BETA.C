/* Program : BETA.C ver 0.1                   */
/* Author  : Ryu choong hyun                  */
/* Date    : 94.7.15.		      	      */
/* Note    : Beta Probability Distribution &  */
/*           Percentage Point                 */

#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>

#define PI M_PI

double logamma(double y);
double gammaf(double x);
void beta_p(double alpha,double beta,double x,double *p,double *d);
double beta_point(double alpha,double beta,double point);
double nor_point(double point);

void main(void)
{
  char opt,c;
  double alpha,beta,point,a,b,pa,pb,mean,var,d1,d2;

  clrscr();
  printf("\n\t    ** Beta  Distribution **\n");
  printf("\n\t1. Probability  2. Percentage Point\n");

  do
  switch (opt=getch()) {
    case '1' :
      printf("\n\t1. P(a<X<b)  2. P(X>a)  3. P(X<a)\n");
      do
      switch (c=getch()) {
      case '1' : printf("\n\t    Alpha = ");
		 scanf("%lf",&alpha);
		 printf("\n\t    Beta  = ");
		 scanf("%lf",&beta);
		 printf("\n\t    From  = ");
		 scanf("%lf",&a);
		 printf("\n\t    To    = ");
		 scanf("%lf",&b);
		 beta_p(alpha,beta,a,&pa,&d1);
		 beta_p(alpha,beta,b,&pb,&d2);
		 printf("\n\n\t  ** Beta  Distribution **\n\n");
		 printf("\t    P(%f \x3c X \x3c %f) = %f\n\n",a,b,pb-pa);
		 break;
      case '2' : printf("\n\t   Alpha = ");
		 scanf("%lf",&alpha);
		 printf("\n\t   Beta  = ");
		 scanf("%lf",&beta);
		 printf("\n\t    a    = ");
		 scanf("%lf",&a);
		 beta_p(alpha,beta,a,&pa,&d1);
		 printf("\n\n\t  ** Beta  Distribution **\n\n");
		 printf("\t    P(X \x3e %f) = %lf\n\n",a,1-pa);
		 break;
      case '3' : printf("\n\t    Alpha = ");
		 scanf("%lf",&alpha);
		 printf("\n\t    Beta  = ");
		 scanf("%lf",&beta);
		 printf("\n\t     a    = ");
		 scanf("%lf",&a);
		 beta_p(alpha,beta,a,&pa,&d1);
		 printf("\n\n\t  ** Beta  Distribution **\n\n");
		 printf("\t    P(X \x3c %f) = %lf\n\n",a,pa);
		 break;
      default :  break;
    } while (!(c=='1' || c=='2' || c=='3'));
    break;
    case '2' : printf("\n\t1. Lower Percentage Point"
		      "  2. Upper Percentage Point\n\n");
	       do
	       switch (opt=getch()) {
		 case '1' :
		 case '2' : printf("\n\t    Alpha = ");
			    scanf("%lf",&alpha);
			    printf("\n\t    Beta  = ");
			    scanf("%lf",&beta);
			    printf("\n\tPercentage Point : ");
			    scanf("%lf",&point);
			    break;
		 default  : break;
	       } while (!(opt=='1' || opt=='2'));

	       printf("\n\n\t** Beta  Distribution **\n\n");
	       if (opt=='1') {
		 a=beta_point(alpha*2,beta*2,point);
		 printf("\t %5.2f Percentage = %f\n\n",point,a);
	       }
	       else {
		 a=beta_point(alpha*2,beta*2,1-point);
		 printf("\t %5.2f Percentage = %f\n\n",point,a);
	       }
	       break;
    default  : break;
  } while (!(opt=='1' || opt=='2'));

  mean=alpha/(alpha+beta);
  var=alpha*beta/(pow(alpha+beta,2)*(alpha+beta+1));
  printf("\t Mean = %f Variance = %f",mean,var);
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

void beta_p(double alpha,double beta,double x,double *p,double *d)
{
  int i;
  double a,b,xx,sum,e;

  sum=1;
  e=1;
  a=alpha;
  b=beta;
  xx=x;

  if (a<=0 || b<=0) {
    fprintf(stderr,"Alpha and Beta is non-positive !!!");
    exit(1);
  }

  if (x>=.5000000012345) {
    a=beta;
    b=alpha;
    xx=1-x;
  }

  *d=exp((a-1)*log(xx)+(b-1)*log(1-xx)+logamma(a+b)-logamma(a)-logamma(b));

  for (i=1;i<=61;i++) {
    e*=(a+b+i-1)*xx/(a+i);
    sum+=e;
  }

  if (e<0.1e-13) {
    *p=*d*sum*xx*(1-xx)/a;
    if (x>=0.5000000012345) *p=1-*p;
  }
  else {
    fprintf(stderr,"Not Converged !!!");
    exit(1);
  }
}

double beta_point(double alpha,double beta,double point)
{
  double a,b,c,d,f1,f2,q,n,x,m,r,eps;

  m=0;
  eps=1.0e-7;
  f1=beta;
  f2=alpha;

  if (beta-1<=0) {
    q=point*.5;
    n=alpha;

    if (n-2>0) {
      a=n;
      c=nor_point(1-q);
      b=c*c;
      x=c+(b+1)*c/(4*a)+
	((5*b+16)*b+3)*c/(96*pow(a,2))+
	(((3*b+19)*b+17)*b-15)*c/(384*pow(a,3))+
	((((79*b+776)*b+1482)*b-1920)*b-945)*c/(92160*pow(a,4))+
	(((((27*b+339)*b+930)*b-1782)*b-765)*b+17955)*c/(368640*pow(a,5));
    }
    else if (n-2<0) x=sin(PI*(.5-q))/cos(PI*(.5-q));
    else {
      a=pow(1-2*q,2);
      x=sqrt(2*a/(1-a));
    }
    x*=x;

    if (alpha==1 && beta!=1) x=1/x;
    x=f2/(f2+f1*fabs(x));

    if (n<=2 && (f1==1 || f2==1)) return x;

    do {
	beta_p(alpha/2,beta/2,x,&q,&d);
	a=d*((.5*f2-1)/x+(.5*f1-1)/(1-x));
	b=pow(d,2)-2*(q-point);
	r=(b<=0) ? -d/a : 2*(q-point)/(-d-sqrt(b));
	x+=r;
	if (++m-10>0) break;
    } while (fabs(r)-eps>0);
  }
  else {
    if (alpha-1>0) {
      x=nor_point(1-point);
      a=1-(2/(9*f1));
      b=1-(2/(9*f2));
      c=pow(b,2)-(1-b)*pow(x,2);
      d=sqrt(fabs(pow(a,2)*pow(b,2)-c*(pow(a,2)-(1-a)*pow(x,2))));

      if (point>.5) d=-d;

      x=pow((a*b+d)/c,3);
    }
    else {
      q=(1-point)*.5;
      n=beta;

      if (n-2>0) {
	a=n;
	c=nor_point(1-q);
	b=c*c;
	x=c+(b+1)*c/(4*a)+
	  ((5*b+16)*b+3)*c/(96*pow(a,2))+
	  (((3*b+19)*b+17)*b-15)*c/(384*pow(a,3))+
	  ((((79*b+776)*b+1482)*b-1920)*b-945)*c/(92160*pow(a,4))+
	  (((((27*b+339)*b+930)*b-1782)*b-765)*b+17955)*c/(368640*pow(a,5));
      }
      else if (n-2<0) x=sin(PI*(.5-q))/cos(PI*(.5-q));
      else {
	a=pow(1-2*q,2);
	x=sqrt(2*a/(1-a));
      }
      x*=x;

      if (alpha==1 && beta!=1) x=1/x;
    }
    x=f2/(f2+f1*fabs(x));

    if (n<=2 && (f1==1 || f2==1)) return x;

    do {
      beta_p(alpha/2,beta/2,x,&q,&d);
      a=d*((.5*f2-1)/x+(.5*f1-1)/(1-x));
      b=pow(d,2)-2*(q-point);
      r=(b<=0) ? -d/a : 2*(q-point)/(-d-sqrt(b));
      x+=r;
      if (++m-10>0) break;
    } while (fabs(r)-eps>0);
  }
  return x;
}

double nor_point(double point)
{
  double q,r,px;
  double a[8]={3.3871328727963666080e0,
	       1.3314166789178437745e2,
	       1.9715909503065514427e3,
	       1.3731693765509461125e4,
	       4.5921953931549871457e4,
	       6.7265770927008700853e4,
	       3.3430575583588128105e4,
	       2.5090809287301226727e3};
  double b[7]={4.2313330701600911252e1,
	       6.8718700749205790830e2,
	       5.3941960214247511077e3,
	       2.1213794301586595867e4,
	       3.9307895800092710610e4,
	       2.8729085735721942674e4,
	       5.2264952788528545610e3};
  double c[8]={1.42343711074968357734e0,
	       4.63033784615654529590e0,
	       5.76949722146069140550e0,
	       3.64784832476320460504e0,
	       1.27045825245236838258e0,
	       2.41780725177450611770e-1,
	       2.27238449892691845833e-2,
	       7.74545014278341407640e-4};
  double d[7]={2.05319162663775882187e0,
	       1.67638483018380384940e0,
	       6.89767334985100004550e-1,
	       1.48103976427480074590e-1,
	       1.51986665636164571966e-2,
	       5.47593808499534494600e-4,
	       1.05075007164441684324e-9};
  double e[8]={6.65790464350110377720e0,
	       5.46378491116411436990e0,
	       1.78482653991729133580e0,
	       2.96560571828504891230e-1,
	       2.65321895265761230930e-2,
	       1.24266094738807843860e-3,
	       2.71155556874348757815e-5,
	       2.01033439929228813265e-7};
  double f[7]={5.99832206555887937690e-1,
	       1.36929880922735805310e-1,
	       1.48753612908506148525e-2,
	       7.86869131145613259100e-4,
	       1.84631831751005468180e-5,
	       1.42151175831644588870e-7,
	       2.04426310338993978564e-15};

  q=point-.5;

  if (fabs(q)<=.425) {
    r=.425*.425-q*q;
    px=(((((((a[7]*r+a[6])*r+a[5])*r+a[4])*r+a[3])*r+a[2])*r+a[1])*r+a[0])/
       (((((((b[6]*r+b[5])*r+b[4])*r+b[3])*r+b[2])*r+b[1])*r+b[0])*r+1);
    return (q*px);
  }
  else {
    r=(q<0.) ? point : 1-point;
    if (r<=0) return 0;
    r=sqrt(-log(r));

    if (r<=5.0) {
      r=r-1.6;
      px=(((((((c[7]*r+c[6])*r+c[5])*r+c[4])*r+c[3])*r+c[2])*r+c[1])*r+c[0])/
	 (((((((d[6]*r+d[5])*r+d[4])*r+d[3])*r+d[2])*r+d[1])*r+d[0])*r+1);
    }
    else {
      r=r-5.0;
      px=(((((((e[7]*r+e[6])*r+e[5])*r+e[4])*r+e[3])*r+e[2])*r+e[1])*r+e[0])/
	 (((((((f[6]*r+f[5])*r+f[4])*r+f[3])*r+f[2])*r+f[1])*r+f[0])*r+1);
    }
  }
  return (q<0) ? -px : px;
}


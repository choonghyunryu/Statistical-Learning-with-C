/* Program : ESTDMEAN.C ver 0.1                   */
/* Author  : Ryu choong hyun                      */
/* Date    : 94.7.18.		                  */
/* Note    : Interval Estimate of Difference Mean */

#include <stdio.h>
#include <conio.h>
#include <math.h>

#define PI M_PI

double nor_point(double point);
double t_q(int df,double t);
double t_point(int df,double point);
double t_sub(int df,double p);

void est_dmean1(int n1,int n2,double mean1,double mean2,
		double var1,double var2,double alpha,double *l,double *u);
void est_dmean2(int n1,int n2,double mean1,double mean2,
		double s1,double s2,double alpha,double *l,double *u);

void main(void)
{
  char opt;
  int n1,n2;
  double mean1,mean2,var1,var2,alpha,lower,upper;

  clrscr();
  printf("  *** Inverval Estimation of Difference Between \xe6 ***\n\n");
  printf("\t      1.Known \xe5\xfd  2.Unknown \xe5\xfd\n\n");

  do
  switch (opt=getch()) {
    case '1': printf("\t     Number of Sample(N1)  = ");
	      scanf("%d",&n1);
	      printf("\t     Mean of Sample(Mean1) = ");
	      scanf("%lf",&mean1);
	      printf("\t     Variance1             = ");
	      scanf("%lf",&var1);
	      printf("\t     Number of Sample(N2)  = ");
	      scanf("%d",&n2);
	      printf("\t     Mean of Sample(Mean2) = ");
	      scanf("%lf",&mean2);
	      printf("\t     Variance2             = ");
	      scanf("%lf",&var2);
	      printf("\t     Condidence Level = ");
	      scanf("%lf",&alpha);
	      est_dmean1(n1,n2,mean1,mean2,var1,var2,alpha,&lower,&upper);
	      break;
    case '2': printf("\t     Number of Sample(N1)  = ");
	      scanf("%d",&n1);
	      printf("\t     Mean of Sample(Mean1) = ");
	      scanf("%lf",&mean1);
	      printf("\t     Sample's Variance1    = ");
	      scanf("%lf",&var1);
	      printf("\t     Number of Sample(N2)  = ");
	      scanf("%d",&n2);
	      printf("\t     Mean of Sample(Mean2) = ");
	      scanf("%lf",&mean2);
	      printf("\t     Sample's Variance2    = ");
	      scanf("%lf",&var2);
	      printf("\t     Condidence Level  = ");
	      scanf("%lf",&alpha);
	      est_dmean2(n1,n2,mean1,mean2,var1,var2,alpha,&lower,&upper);
    default : break;
  } while (!(opt=='1' || opt=='2'));

  printf("\n\n\tN1    = %5d\n",n1);
  printf("\tN2    = %5d\n",n2);
  printf("\tMEAN1 = %lf\n",mean1);
  printf("\tMEAN2 = %lf\n",mean2);

  if (opt=='2') {
    printf("\tS1\xfd   = %lf\n",var1);
    printf("\tS2\xfd   = %lf\n",var2);
  }
  else {
    printf("\t\xe5\xfd    = %lf\n",var1);
    printf("\t\xe5\xfd    = %lf\n",var2);
  }
  printf("\n\tCondidence Level = %5.2f\n",alpha*100);
  printf("\tInterval of Estimation's Lower = %lf\n",lower);
  printf("\tInterval of Estimation's Upper = %lf\n",upper);
  getch();
}

void est_dmean1(int n1,int n2,double mean1,double mean2,double var1,
	       double var2,double alpha,double *l,double *u)
{
  double z;

  z=nor_point((1-alpha)/2);

  *l=(mean1-mean2)-fabs(z)*sqrt(var1/n1+var2/n2);
  *u=(mean1-mean2)+fabs(z)*sqrt(var1/n1+var2/n2);
}

void est_dmean2(int n1,int n2,double mean1,double mean2,double s1,double s2,
		double alpha,double *l,double *u)
{
  double t;

  t=(n1>=30 || n2>=30) ? nor_point((1-alpha)/2) :
    (alpha<=.5) ? -t_point(n1+n2-2,alpha/2) : t_point(n1+n2-2,(1-alpha)/2);

  *l=(mean1-mean2)-fabs(t)*sqrt(s1/n1+s2/n2);
  *u=(mean1-mean2)+fabs(t)*sqrt(s1/n1+s2/n2);
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

double t_point(int df,double point)
{
  double u,y[5];

  if(df>5 && (point>5.e-3 || df>13)) {
  u=nor_point(1-point);
  y[0]=(pow(u,3)+u)/pow(2,2);
  y[1]=(5*pow(u,5)+16*pow(u,3)+3*u)/(pow(2,5)*3);
  y[2]=(3*pow(u,7)+19*pow(u,5)+17*pow(u,3)-15*u)/(pow(2,7)*3);
  y[3]=(79*pow(u,9)+776*pow(u,7)+1482*pow(u,5)-1920*pow(u,3)-945*u)/
       (pow(2,11)*pow(3,2)*5);
  y[4]=(27*pow(u,11)+339*pow(u,9)+930*pow(u,7)-1782*pow(u,5)-765*pow(u,3)+
	17955*u)/(pow(2,13)*pow(3,2)*5);

  return u+y[0]/df+y[1]/pow(df,2)+y[2]/pow(df,3)+y[3]/pow(df,4)+
	 y[4]/pow(df,5);
  }
  else return t_sub(df,point);
}

double t_sub(int df,double p)
{
  double  eps,q,point=0,s=10000;

  if (df==1 && p<0.01 && p>=0.001) eps = 1.0e-4;
  else if(df==2 && p< 0.0001)      eps = 1.0e-4;
  else if(df==1 && p< 0.001)       eps = 1.0e-2;
  else                             eps = 1.0e-5;

  while (1) {
    point+=s;
    if (s<=eps) return point;
    q=t_q(df,point)-p;
    if (q==0) return point;
    if (q<0) {
      point-=s;
      s/=10;
    }
  }
}

double t_q(int df,double t)
{
  int i;
  double theta,s,c,sum=0,temp=1,p,tx;

  tx=t;
  t=(t<0) ? -t : t;

  theta=atan(t/sqrt(df));
  c=sqrt(df/(df+t*t));
  s=sqrt(1-df/(df+t*t));

  for (i=df%2+4;i<=df;i+=2) {
    temp*=(i-3.)/(i-2.);
    sum+=temp*s*pow(c,i-2);
  }

  if (df & 1) {
    if (df==1) p=1-(.5+theta/PI);
    else p=1-(.5+(theta+s*c+sum)/PI);
  }
  else p=1-(.5+(s+sum)/2);

  return (tx<0) ? 1-p : p;
}

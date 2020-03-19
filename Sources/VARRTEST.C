/* Program : VARRTEST.C                        */
/* Author  : Ryu choong hyun                   */
/* Date    : 94.7.21.		               */
/* Note    : Statistical Hypothesis Testing    */
/*           Ratio Between Population Variance */

#include <stdio.h>
#include <conio.h>
#include <math.h>

#define PI M_PI

int vr_one_test(char opt,int n1,int n2,double s1,double s2,
		double alpha,double *cv);
int vr_two_test(int n1,int n2,double s1,double s2,double alpha,
		double *cv1,double *cv2);

double nor_point(double z);
double f_q(int df1,int df2,double f);
double f_point(int df1,int df2,double p);

void main(void)
{
  char opt;
  int n1,n2,h0;
  double s1,s2,alpha,cv1,cv2;

  clrscr();
  printf("\t   ** Statistical Hypothesis Testing Ratio \xe5\xfd **\n\n");
  printf("\t\t1. H0 : \xe5\xfd%d \xf3 \xe5\xfd%d"
	     "   H1 : \xe5\xfd%d > \xe5\xfd%d\n",1,2,1,2);
  printf("\t\t2. H0 : \xe5\xfd%d \xf2 \xe5\xfd%d"
	     "   H1 : \xe5\xfd%d < \xe5\xfd%d\n",1,2,1,2);
  printf("\t\t3. H0 : \xe5\xfd%d = \xe5\xfd%d"
	     "   H1 : \xe5\xfd%d \xd8 \xe5\xfd%d\n\n",1,2,1,2);
  do
  switch (opt=getch()) {
    case '1':
    case '2':
    case '3': printf("\t\tN1                 = ");
	      scanf("%d",&n1);
	      printf("\t\tN2                 = ");
	      scanf("%d",&n2);
	      printf("\t\t\S\xfd%d                = ",1);
	      scanf("%lf",&s1);
	      printf("\t\t\S\xfd%d                = ",2);
	      scanf("%lf",&s2);
	      printf("\t\tSignificance Level = ");
	      scanf("%lf",&alpha);
	      break;
    default : break;
  } while (!(opt=='1' || opt=='2' || opt=='3'));

  if (opt=='3') {
    h0=vr_two_test(n1,n2,s1,s2,alpha,&cv1,&cv2);
    printf("\n\t\tLower Cirtical Value = %lf\n",cv1);
    printf("\t\tUpper Cirtical Value = %lf\n",cv2);
  }
  else {
    h0=vr_one_test(opt,n1,n2,s1,s2,alpha,&cv1);
    printf("\n\t\tCirtical Value = %lf\n",cv1);
  }

  if (h0==1) printf("\t\tNull Hypothesis(H0) is Accept\n");
  else printf("\t\tNull Hypothesis(H0) is Reject\n");
  getch();
}

int vr_one_test(char opt,int n1,int n2,double s1,double s2,
		double alpha,double *cv)
{
  switch (opt) {
    case '1': *cv=f_point(n1-1,n2-1,alpha);
	      break;
    case '2': *cv=f_point(n1-1,n2-1,1-alpha);
	      break;
    default : break;
  }

  if (opt=='1') return (s1/s2<=*cv) ? 1 : 0;
  else return (s1/s2>=*cv) ? 1 : 0;
}

int vr_two_test(int n1,int n2,double s1,double s2,
		double alpha,double *cv1,double *cv2)
{
  *cv1=f_point(n1-1,n2-1,1-alpha/2);
  *cv2=f_point(n1-1,n2-1,alpha/2);

  return (*cv1<s1/s2 && s1/s2<*cv2) ? 1 : 0;
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

double f_q(int df1,int df2,double f)
{
  int i;
  double theta,s,c,sum=0,temp=1;

  theta=atan(sqrt(df1*f/df2));
  c=sqrt(1/(1+df1*f/df2));
  s=sqrt(1-c*c);

  if ((df1*df2)%2==0 && df1 & 1) return 1-f_q(df2,df1,1/f);

  if ((df1*df2)%2==0) {
    for (i=4;i<=df1;i+=2) {
      temp*=(df2+i-4.)/(i-2.);
      sum+=temp*pow(s,i-2);
    }
    return pow(c,df2)*(1+sum);
  }
  else {
    for (i=3;i<=df2;i+=2) {
      if (i != 3) temp*=((i-3.)/(i-2.));
      sum+=temp*s*pow(c,i-2);
    }
    for (i=3;i<=df1;i+=2) {
      if (df2==1 && i==3) temp=1;
      else temp*=(i+df2-4.)/(i-2.);
      sum-=temp*pow(c,df2)*pow(s,i-2);
    }
    return 1-2*(theta+sum)/PI;
  }
}

double f_point(int df1,int df2,double p)
{
  double gamma,h,delta,u,z,t1,t2,t3,t4,eps,s,fw,qe;

  if (df1<=10 || df2<=10) {
    eps=1.0e-5;
    if (df2==1) eps=1.0e-4;
    s=1000;
    fw=0;

    while (1) {
      fw+=s;
      if (s<=eps) return fw;
      qe=f_q(df1,df2,fw)-p;
      if (qe==0) return fw;
      if (qe<0) {
	fw-=s;
	s/=10;
      }
    }
  }
  else {
    gamma=1./df1+1./df2;
    delta=1./df1-1./df2;
    h=2/gamma;
    u=nor_point(1-p);

    t1=(u*u+2)*delta/6.;
    t2=((pow(u,3)+3*u)/(12*h)+(pow(u,3)+11*u)/144.*h*delta*delta)/sqrt(h);
    t3=(pow(u,4)+9*u*u+8)*delta/(60*h)-(3*pow(u,4)+7*u*u-16)*h*pow(delta,3)
       /6480.;
    t4=((pow(u,5)+20*pow(u,3)+15*u)/(480*h*h)+(pow(u,5)+44*pow(u,3)+183*u)*
	 pow(delta,4)/2880.+(9*pow(u,5)-284*pow(u,3)-1513*u)*h*h*pow(delta,4)
	 /622080.)/sqrt(h);
    z=u/sqrt(h)-t1+t2-t3+t4;

    return exp(2*z);
  }
}

/* Program : PTEST.C                                  */
/* Author  : Ryu choong hyun                          */
/* Date    : 94.7.23.		                      */
/* Note    : Statistical Hypothesis Testing           */
/*           Difference Between Population Proportion */

#include <stdio.h>
#include <conio.h>
#include <math.h>

int pd_one_test(char opt,int n1,int n2,int x1,int x2,double alpha,double *cv);
int pd_two_test(int n1,int n2,int x1,int x2,double alpha,
		double *cv1,double *cv2);

double nor_point(double z);

void main(void)
{
  char opt;
  int n1,n2,x1,x2,h0;
  double alpha,cv1,cv2;

  clrscr();
  printf("    ** Statistical Hypothesis Testing Difference Between p **\n\n");
  printf("\t\t1. H0 : p1 \xf3 p2  H1 : p1 > p2\n");
  printf("\t\t2. H0 : p1 \xf2 p2  H1 : p1 < p2\n");
  printf("\t\t3. H0 : p1 = p2  H1 : p1 \xd8 p2\n\n");
  do
  switch (opt=getch()) {
    case '1':
    case '2':
    case '3': printf("\t\tN1                 = ");
	      scanf("%d",&n1);
	      printf("\t\tN2                 = ");
	      scanf("%d",&n2);
	      printf("\t\t\X1                 = ");
	      scanf("%d",&x1);
	      printf("\t\t\X2                 = ");
	      scanf("%d",&x2);
	      printf("\t\tSignificance Level = ");
	      scanf("%lf",&alpha);
	      break;
    default : break;
  } while (!(opt=='1' || opt=='2' || opt=='3'));

  if (opt=='3') {
    h0=pd_two_test(n1,n2,x1,x2,alpha,&cv1,&cv2);
    printf("\n\t\tLower Cirtical Value = %lf\n",cv1);
    printf("\t\tUpper Cirtical Value = %lf\n",cv2);
  }
  else {
    h0=pd_one_test(opt,n1,n2,x1,x2,alpha,&cv1);
    printf("\n\t\tCirtical Value = %lf\n",cv1);
  }

  if (h0==1) printf("\t\tNull Hypothesis(H0) is Accept\n");
  else printf("\t\tNull Hypothesis(H0) is Reject\n");
  getch();
}

int pd_one_test(char opt,int n1,int n2,int x1,int x2,double alpha,double *cv)
{
  double p1,p2;

  p1=x1/(double)n1;
  p2=x2/(double)n2;

  switch (opt) {
    case '1': *cv=nor_point(1-alpha)*sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2);
	      break;
    case '2': *cv=nor_point(alpha)*sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2);
	      break;
    default : break;
  }

  if (opt=='1') return (p1-p2<=*cv) ? 1 : 0;
  else return (p1-p2>=*cv) ? 1 : 0;
}

int pd_two_test(int n1,int n2,int x1,int x2,double alpha,
	       double *cv1,double *cv2)
{
  double p,p1,p2;

  p1=x1/(double)n1;
  p2=x2/(double)n2;
  p=(x1+x2)/(double)(n1+n2);

  *cv1=nor_point(alpha/2)*sqrt(p*(1-p)*(1./n1+1./n2));
  *cv2=nor_point(1-alpha/2)*sqrt(p*(1-p)*(1./n1+1./n2));

  return (*cv1<p1-p2 && p1-p2<*cv2) ? 1 : 0;
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

/* Program : VARTEST.C                       */
/* Author  : Ryu choong hyun                 */
/* Date    : 94.7.21.		             */
/* Note    : Statistical Hypothesis Testing  */
/*           Population Variance             */

#include <stdio.h>
#include <conio.h>
#include <math.h>

int v_one_test(char opt,int n,double var,double s,double alpha,double *cv);
int v_two_test(int n,double var,double s,double alpha,
	       double *cv1,double *cv2);

double nor_point(double z);
double chi_point(int df, double point);

void main(void)
{
  char opt;
  int n,h0;
  double var,s,alpha,cv1,cv2;

  clrscr();
  printf("\t     ** Statistical Hypothesis Testing \xe5\xfd **\n\n");
  printf("\t\t1. H0 : \xe5\xfd \xf3 \xe5\xfd%d"
	     "   H1 : \xe5\xfd > \xe5\xfd%d\n",0,0);
  printf("\t\t2. H0 : \xe5\xfd \xf2 \xe5\xfd%d"
	     "   H1 : \xe5\xfd < \xe5\xfd%d\n",0,0);
  printf("\t\t3. H0 : \xe5\xfd = \xe5\xfd%d"
	     "   H1 : \xe5\xfd \xd8 \xe5\xfd%d\n\n",0,0);
  do
  switch (opt=getch()) {
    case '1':
    case '2':
    case '3': printf("\t\t\xe5\xfd                 = ");
	      scanf("%lf",&var);
	      printf("\t\tN                  = ");
	      scanf("%d",&n);
	      printf("\t\t\S\xfd                 = ");
	      scanf("%lf",&s);
	      printf("\t\tSignificance Level = ");
	      scanf("%lf",&alpha);
	      break;
    default : break;
  } while (!(opt=='1' || opt=='2' || opt=='3'));

  if (opt=='3') {
    h0=v_two_test(n,var,s,alpha,&cv1,&cv2);
    printf("\n\t\tLower Cirtical Value = %lf\n",cv1);
    printf("\t\tUpper Cirtical Value = %lf\n",cv2);
  }
  else {
    h0=v_one_test(opt,n,var,s,alpha,&cv1);
    printf("\n\t\tCirtical Value = %lf\n",cv1);
  }

  if (h0==1) printf("\t\tNull Hypothesis(H0) is Accept\n");
  else printf("\t\tNull Hypothesis(H0) is Reject\n");
  getch();
}

int v_one_test(char opt,int n,double var,double s,double alpha,double *cv)
{
  switch (opt) {
    case '1': *cv=(n-1)*var/chi_point(n-1,1-alpha);
	      break;
    case '2': *cv=(n-1)*var/chi_point(n-1,alpha);
	      break;
    default : break;
  }

  if (opt=='1') return (s<=*cv) ? 1 : 0;
  else return (s>=*cv) ? 1 : 0;
}

int v_two_test(int n,double var,double s,double alpha,
	       double *cv1,double *cv2)
{
  *cv1=(n-1)*var/chi_point(n-1,alpha/2);
  *cv2=(n-1)*var/chi_point(n-1,1-alpha/2);

  return (*cv1<s && s<*cv2) ? 1 : 0;
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

double chi_point(int df, double point)
{
  int i;
  double u,c[11],temp=0;

  u=nor_point(1-point);

  switch (df) {
    case 1  : return pow(nor_point(point*.5),2);
    case 2  : return -2*log(point);
    default :
  c[0]=1;
  c[1]=u;
  c[2]=(pow(u,2)-1)/3;
  c[3]=(pow(u,3)-7*u)/(pow(2,2)*pow(3,2));
  c[4]=(-3*pow(u,4)-7*pow(u,2)+16)/(2*pow(3,4)*5);
  c[5]=(9*pow(u,5)+256*pow(u,3)-433*u)/(pow(2,5)*pow(3,5)*5);
  c[6]=(12*pow(u,6)-243*pow(u,4)-923*pow(u,2)+1472)/(pow(2,3)*pow(3,6)*5*7);
  c[7]=(-3753*pow(u,7)-4353*pow(u,5)+289517*pow(u,3)+289717*u)/
       (pow(2,7)*pow(3,8)*pow(5,2)*7);
  c[8]=(270*pow(u,8)+4614*pow(u,6)-9513*pow(u,4)-104989*pow(u,2)+35968)/
       (pow(2,4)*pow(3,9)*pow(5,2)*7);
  c[9]=(-5139*pow(u,9)-547848*pow(u,7)-2742210*pow(u,5)+7016224*pow(u,3)+
	37501325*u)/(pow(2,11)*pow(3,10)*pow(5,2)*7);
  c[10]=(-364176*pow(u,10)+6208146*pow(u,8)+125735778*pow(u,6)+303753831*
	 pow(u,4)-672186949*pow(u,2)-2432820224)/
	(pow(2,7)*pow(3,13)*pow(5,3)*7*11);
  c[11]=(199112985*pow(u,11)+1885396761*pow(u,9)-31857434154*pow(u,7)-
	 287542736226*pow(u,5)-556030221167*pow(u,3)+487855454729*u)/
	(pow(2,13)*pow(3,14)*pow(5,3)*pow(7,2)*11);

  for (i=0;i<=11;i++)
    temp+=c[i]*pow(2./df,i/2.);
  return df*temp;
  }
}

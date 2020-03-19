/* Program : SPEARMAN.C ver 0.1               */
/* Author  : Ryu choong hyun                  */
/* Date    : 94.8.20.	                      */
/* Note    : Spearman's Rank Correlation Test */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define PI M_PI
#define N 50

void sort(int th,int left,int right);
void sort1(int th,int left,int right);
void rank(int th);
void get_txy(void);
double get_rs(char opt,double sum_d2);
int spearman_test(double t,double alpha,char opt);
double nor_point(double point);
double t_q(int df,double t);
double t_point(int df,double point);
double t_sub(int df,double p);

int n;
typedef struct{
  int id;
  double no;
  double x;
} data;

data x[N][2];
double d[N],tx,ty;

void main(int argc,char *argv[])
{
  FILE *stream;
  char opt;
  int i,j,h0;
  double num,rs,s,t,sum_ds=0,alpha;

  clrscr();

  if (argc<=1) {
    puts("Usage : spearman datafile");
    exit(EXIT_FAILURE);
  }

  stream=fopen(argv[1],"rt");
  if (stream==NULL) printerrmsg("File not found !!");

  fscanf(stream,"%d\n",&n);

  for (i=0;i<n;i++)
    for (j=0;j<2;j++) {
      fscanf(stream,"%lf\n",&num);
      x[i][j].x=num;
      x[i][j].id=i+1;
    }

  for (i=0;i<2;i++) {
    sort(i,0,n-1);
    rank(i);
  }

  get_txy();

  for (i=0;i<2;i++) sort1(i,0,n-1);

  for (i=0;i<n;i++) {
    d[i]=x[i][0].no-x[i][1].no;
    sum_ds+=pow(d[i],2);
  }

  rs=(tx==0 && ty==0) ? get_rs('0',sum_ds) : get_rs('1',sum_ds);
  s=sqrt((1-pow(rs,2))/(n-2));
  t=rs/s;

  printf("\t    ** Spearman's Rank Correlation Test **\n\n");
  printf("\t\t1. H0 : Rho = 0  H1 : Rho > 0\n");
  printf("\t\t2. H0 : Rho = 0  H1 : Rho < 0\n");
  printf("\t\t3. H0 : Rho = 0  H1 : Rho \xd8 0\n");
  do
  switch (opt=getch()) {
    case '1':
    case '2':
    case '3': printf("\n\t\t   Significance Level = ");
	      scanf("%lf",&alpha);
	      break;
    default : break;
  } while (!(opt=='1' || opt=='2' || opt=='3'));

  printf("\tObj      X      Rank1      Y      Rank2    Diff    Diff\xfd\n");
  printf("\t--------------------------------------------------------\n");
  for (i=0;i<n;i++)
    printf("\t%2d  %8.2f  %6.1f  %8.2f  %6.1f  %7.1f  %7.2f\n",
	   i+1,x[i][0].x,x[i][0].no,x[i][1].x,x[i][1].no,d[i],pow(d[i],2));
  printf("\t---------------------------------------------------------\n");
  printf("\tSUM\t\t\t\t\t\t%8.2f\n\n",sum_ds);

  printf("\t\tNumber of Sample       = %d\n",n);
  printf("\t\tSpearman's Correlation = %lf\n",rs);
  printf("\t\tT Statistic            = %lf\n",t);

  if (n>=30) {
    h0=(opt=='3') ? spearman_test(t,alpha,'0') : spearman_test(t,alpha,opt);
    if (h0==1) printf("\n\t\tNull Hypothesis(H0) is Accept\n");
    else printf("\n\t\tNull Hypothesis(H0) is Reject\n");
  }

  getch();
}

void sort(int th,int left,int right)
{
  int i=left,j=right;
  double temp,mid=x[((left+right)/2)][th].x;

  do {
    while (x[i][th].x < mid && i<right)
      i++;
    while (mid < x[j][th].x && j>left)
      j--;
    if (i<=j) {
      temp=x[i][th].x;
      x[i][th].x=x[j][th].x;
      x[j][th].x=temp;

      temp=x[i][th].id;
      x[i][th].id=x[j][th].id;
      x[j][th].id=temp;
      i++;
      j--;
    }
  } while (i<=j);
  if (left<j)
    sort(th,left,j);
  if (i<right)
    sort(th,i,right);
}

void sort1(int th,int left,int right)
{
  int i=left,j=right;
  double temp,mid=x[((left+right)/2)][th].id;

  do {
    while (x[i][th].id < mid && i<right)
      i++;
    while (mid < x[j][th].id && j>left)
      j--;
    if (i<=j) {
      temp=x[i][th].id;
      x[i][th].id=x[j][th].id;
      x[j][th].id=temp;

      temp=x[i][th].x;
      x[i][th].x=x[j][th].x;
      x[j][th].x=temp;

      temp=x[i][th].no;
      x[i][th].no=x[j][th].no;
      x[j][th].no=temp;
      i++;
      j--;
    }
  } while (i<=j);
  if (left<j)
    sort1(th,left,j);
  if (i<right)
    sort1(th,i,right);
}

void rank(int th)
{
  int i,j,k,same,sum;

  for (i=0;i<n;i++) {
    sum=0;
    same=1;
    if (x[i][th].x!=x[i+1][th].x) x[i][th].no=i+1;
    else {
      for (j=i;;j++) {
	if (x[j+1][th].x!=x[j+2][th].x) {
	  same++;
	  for (k=i;k<i+same;k++) sum+=k+1;
	  for (k=i;k<i+same;k++) x[k][th].no=(double)sum/same;
	  i=i+same-1;
	  break;
	}
	same++;
      }
    }
  }
}

void get_txy(void)
{
  int i,t=1;
  double temp=0;

  for (i=0;i<n;i++) {
    if (temp==x[i][0].no) t++;
    else {
      tx+=pow(t,3)-t;
      t=1;
    }
    temp=x[i][0].no;
  }

  for (i=0;i<n;i++) {
    if (temp==x[i][1].no) t++;
    else {
      ty+=pow(t,3)-t;
      t=1;
    }
    temp=x[i][1].no;
  }
}

double get_rs(char opt,double sum_d2)
{
  return (opt=='0') ? 1-6*sum_d2/(pow(n,3)-n) :
		      ((pow(n,3)-n)-6*sum_d2-(tx+ty)/2)/
		      sqrt(pow(pow(n,3)-n,2)-(tx+ty)*(pow(n,3)-n)+tx*ty);
}

int spearman_test(double t,double alpha,char opt)
{
  double t_p;

  t_p=(opt=='0') ? fabs(t_point(n-2,alpha/2)) : fabs(t_point(n-2,alpha));

  if (opt=='0') return (t>=-t_p && t<=t_p) ? 1 : 0;
  else if (opt=='1') return (t<=t_p) ? 1 : 0;
  else return (t>=-t_p) ? 1 : 0;
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

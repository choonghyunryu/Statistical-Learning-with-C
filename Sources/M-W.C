/* Program : M-W.C ver 0.1     */
/* Author  : Ryu choong hyun   */
/* Date    : 94.8.17.	       */
/* Note    : Mann-Whitney Test */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define FACTOR 30

int mw_test(double alpha,double u,char opt);
void sort(int left,int right);
void rank(void);
double nor_point(double point);

typedef struct{
  int id;
  double no;
  double x;
} data;

int n1,n2;
data x[2*FACTOR+1];

void main(int argc,char *argv[])
{
  FILE *stream;
  char opt;
  int i,j,h0;
  double num,x1[FACTOR],x2[FACTOR],u,sum1=0,alpha;

  clrscr();

  if (argc<=1) {
    puts("Usage : m-w datafile");
    exit(EXIT_FAILURE);
  }

  stream=fopen(argv[1],"rt");
  if (stream==NULL) printerrmsg("File not found !!");

  fscanf(stream,"%d\n",&n1);
  fscanf(stream,"%d\n",&n2);

  for (i=0;i<n2;i++) {
    for (j=0;j<2;j++) {
      if (i<n1) {
	fscanf(stream,"%lf\n",&num);
	if (j==0) x1[i]=num;
	else x2[i]=num;
      }
      else break;
    }
    if (i>=n1) {
      fscanf(stream,"%lf\n",&num);
      x2[i]=num;
    }
  }

  for (i=0;i<n1;i++) {
    x[i].x=x1[i];
    x[i].id=1;
  }
  for (i=0;i<n2;i++) {
    x[n1+i].x=x2[i];
    x[n1+i].id=2;
  }

  sort(0,n1+n2);
  rank();

  for (i=1;i<=n1+n2;i++) if (x[i].id==1) sum1+=x[i].no;

  u=n1*n2+n1*(n1+1)/2-sum1;

  printf("\t\t   ** Mann-Whitney Test **\n\n");
  printf("\t\t1. H0 : M1 \xf3 M2  H1 : M1 > M2\n");
  printf("\t\t2. H0 : M1 \xf2 M2  H1 : M1 < M2\n");
  printf("\t\t3. H0 : M1 = M2  H1 : M1 \xd8 M2\n\n");
  do
  switch (opt=getch()) {
    case '1':
    case '2':
    case '3': printf("\t\t   ID     Data     Rank\n");
	      for (i=1;i<=n1+n2;i++)
		printf("\t\t    %d   %7.2f   %4.1f\n",x[i].id,x[i].x,x[i].no);
	      if (n2<11) {
		printf("\n\t\t\tN1 = %3d\n",n1);
		printf("\t\t\tN2 = %3d\n",n2);
		printf("\t\t\tU0 = %3d\n",(int)u);
	      }
	      else {
		printf("\n\t\t  Significance Level = ");
		scanf("%lf",&alpha);
		h0=mw_test(alpha,u,opt);
		if (h0==1) printf("\n\t\tNull Hypothesis(H0) is Accept\n");
		else printf("\n\t\tNull Hypothesis(H0) is Reject\n");
	      }
	      break;
    default : break;
  } while (!(opt=='1' || opt=='2' || opt=='3'));
  getch();
}

void sort(int left,int right)
{
  int i=left,j=right;
  double temp,mid=x[((left+right)/2)].x;

  do {
    while (x[i].x < mid && i<right)
      i++;
    while (mid < x[j].x && j>left)
      j--;
    if (i<=j) {
      temp=x[i].x;
      x[i].x=x[j].x;
      x[j].x=temp;
      temp=x[i].id;
      x[i].id=x[j].id;
      x[j].id=temp;
      i++;
      j--;
    }
  } while (i<=j);
  if (left<j)
    sort(left,j);
  if (i<right)
    sort(i,right);
}

void rank(void)
{
  int i,j,k,same,sum;

  for (i=1;i<=n1+n2;i++) {
    sum=0;
    same=1;
    if (x[i].x!=x[i+1].x) x[i].no=i;
    else {
      for (j=i;;j++) {
	if (x[j+1].x!=x[j+2].x) {
	  same++;
	  for (k=i;k<i+same;k++) sum+=k;
	  for (k=i;k<i+same;k++) x[k].no=(double)sum/same;
	  i=i+same-1;
	  break;
	}
	same++;
      }
    }
  }
}

int mw_test(double alpha,double u,char opt)
{
  double z,mean,dev,z_p;

  mean=n1*n2/2.;
  dev=sqrt(n1*n2*(n1+n2+1)/12.);
  z=(u-mean)/dev;
  z_p=(opt=='0') ? fabs(nor_point(alpha/2)) : fabs(nor_point(alpha));

  if (opt=='1') return (z<=z_p) ? 1 : 0;
  else if (opt=='2') return (z>=-z_p) ? 1 : 0;
  else return (z>=-z_p && z<=z_p) ? 1 : 0;
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

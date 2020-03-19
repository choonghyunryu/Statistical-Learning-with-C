/* Program : MREGTEST.C ver 0.1            */
/* Author  : Ryu choong hyun               */
/* Date    : 94.8.3.   	                   */
/* Note    : Multiple Regression Analysis  */
/*           Hypothesis Tets               */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define PI M_PI
#define M 50
#define N 6

void get_xmatrix(void);
void matrix_inverse(void);
void get_bvector(void);
void get_yhat(void);
int get_rank(void);
int bi_one_test(char opt,int nb,int df,double b0,double mse,
		double alpha,double *ts,double *t);
int bi_two_test(int nb,int df,double b0,double mse,
		double alpha,double *ts,double *t);
int y_one_test(char opt,int df,double y,double mse,
	       double alpha,double *ts,double *t);
int y_two_test(int df,double y,double mse,double alpha,double *ts,double *t);

double nor_point(double p);
double t_q(int df,double t);
double t_point(int df,double point);
double t_sub(int df,double p);

int row,col;
double y[M],yhat[M],x[M][N],xmatrix[N][N],b[N],xvector[N],ymean;

void main(int argc,char *argv[])
{
  FILE *stream;
  char job,opt;
  int i,j,rank,dfr,dfe,nb,h0;
  double num,ysum=0,sse=0,mse,sy_x,b0,mu,alpha,l,u,ts,t;

  clrscr();

  if (argc<=1) {
    puts("Usage : mregtest datafile");
    exit(EXIT_FAILURE);
  }

  stream=fopen(argv[1],"rt");
  if (stream==NULL) printerrmsg("File not found !!");

  fscanf(stream,"%d\n",&row);
  fscanf(stream,"%d\n",&col);

  for (i=0;i<row;i++)
    for (j=0;j<col;j++) {
      fscanf(stream,"%lf\n",&num);
      if (j==0) {
	y[i]=num;
	ysum+=y[i];
	x[i][j]=1;
      }
      else x[i][j]=num;
    }

  ymean=ysum/row;
  get_xmatrix();
  matrix_inverse();
  get_bvector();
  get_yhat();
  rank=get_rank();

  for (i=0;i<row;i++) sse+=pow(y[i]-yhat[i],2);
  dfe=row-rank;
  mse=sse/dfe;

  printf("\t *** Multiple Regression Analysis (Hypothesis Test) ***\n\n");
  printf("\t\t\t 1.\xe1i Test  2.\xe6y.x Test\n\n");

  do
  switch (job=getch()) {
    case '1': printf("\t\t1.H0 : \xe1i \xf3 \xe1  H1 : \xe1i > \xe1\n");
	      printf("\t\t2.H0 : \xe1i \xf2 \xe1  H1 : \xe1i < \xe1\n");
	      printf("\t\t3.H0 : \xe1i = \xe1  H1 : \xe1i \xd8 \xe1\n\n");
	      do
		opt=getch();
	      while (!(opt=='1' || opt=='2' || opt=='3'));
	      printf("\t\tNumber of \xe1        = ");
	      scanf("%d",&nb);
	      printf("\t\tb                  = ");
	      scanf("%lf",&b0);
	      printf("\t\tSignificance Level = ");
	      scanf("%lf",&alpha);
	      if (opt=='3') h0=bi_two_test(nb,dfe,b0,mse,alpha,&ts,&t);
	      else  h0=bi_one_test(opt,nb,dfe,b0,mse,alpha,&ts,&t);
	      break;
    case '2': printf("\t\t1.H0 : \xe6y.x \xf3 \xe6  H1 : \xe6y.x > \xe1\n");
	      printf("\t\t2.H0 : \xe6y.x \xf2 \xe6  H1 : \xe6y.x < \xe1\n");
	      printf("\t\t3.H0 : \xe6y.x = \xe6  H1 : \xe6y.x \xd8 \xe1\n\n");
	      do
		opt=getch();
	      while (!(opt=='1' || opt=='2' || opt=='3'));
	      for (i=0;i<col;i++) {
		if (i==0) xvector[i]=1;
		else {
		  printf("\t\tX%d                  = ",i);
		  scanf("%lf",&xvector[i]);
		}
	      }
	      printf("\t\t\xe6                   = ");
	      scanf("%lf",&mu);
	      printf("\t\tSignificance Level  = ");
	      scanf("%lf",&alpha);

	      if (opt=='3') h0=y_two_test(dfe,mu,mse,alpha,&ts,&t);
	      else  h0=y_one_test(opt,dfe,mu,mse,alpha,&ts,&t);
	      break;
    default : break;
  } while (!(job=='1' || job=='2'));

  if (opt=='3') {
    printf("\n\t\tTest Statistic       = %lf\n",ts);
    printf("\t\tLeft Reject Value    = %lf\n",-t);
    printf("\t\tRight Reject Value   = %lf\n",t);
  }
  else {
    t=(opt==2) ? -t : t;
    printf("\n\t\tTest Statistic = %lf\n",ts);
    printf("\t\tReject Value   = %lf\n",t);
  }

  if (h0==1) printf("\t\tNull Hypothesis(H0) is Accept\n");
  else printf("\t\tNull Hypothesis(H0) is Reject\n");

  getch();
}

void get_xmatrix(void)
{
  int i,j,k;

  double sum;

  for (i=0;i<col;i++)
    for (j=0;j<col;j++) {
      sum=0;
      for (k=0;k<row;k++) sum+=x[k][i]*x[k][j];
      xmatrix[i][j]=sum;
    }
}

void matrix_inverse(void)
{
  int i,j,k;
  double temp,u;

  for (k=0;k<col;k++) {
    temp=xmatrix[k][k];
    for (i=0;i<col;i++) xmatrix[k][i]/=temp;
    xmatrix[k][k]=1/temp;
    for (j=0;j<col;j++)
      if (j!=k) {
	u=xmatrix[j][k];
	for (i=0;i<col;i++)
	  if (i!=k) xmatrix[j][i]-=xmatrix[k][i]*u;
	  else xmatrix[j][i]=-u/temp;
      }
  }
}

void get_bvector(void)
{
  int i,j;

  double sum,xy[N];

  for (i=0;i<col;i++) {
    sum=0;
    for (j=0;j<row;j++) sum+=x[j][i]*y[j];
    xy[i]=sum;
  }

  for (i=0;i<col;i++) {
    sum=0;
    for (j=0;j<col;j++) sum+=xmatrix[j][i]*xy[j];
    b[i]=sum;
  }
}

void get_yhat(void)
{
  int i,j;
  double sum;

  for (i=0;i<row;i++) {
    sum=0;
    for (j=0;j<col;j++) sum+=x[i][j]*b[j];
    yhat[i]=sum;
  }
}

int get_rank(void)
{
  int i,j,k,rank=row;
  double xx[M][N];

  for (i=0;i<row;i++)
    for (j=0;j<col;j++) xx[i][j]=x[i][j];

  for (i=0;i<col;i++)
    for (j=i+1;j<row;j++)
      xx[j][i]-=xx[0][i]*xx[j][i]/xx[0][i];

  for (i=0;i<row;i++) {
    k=0;
    for (j=0;j<col;j++) if (xx[i][j]==0) k++;
    if (k==col) rank--;
  }
  return rank;
}

int bi_one_test(char opt,int nb,int df,double b0,double mse,
		double alpha,double *ts,double *t)
{
  *t=fabs(t_point(df,alpha));
  *ts=(b[nb]-b0)/sqrt(mse*xmatrix[nb][nb]);

  if (opt=='1') return (*ts<=*t) ? 1 : 0;
  else return (*ts>=-*t) ? 1 : 0;
}

int bi_two_test(int nb,int df,double b0,double mse,
		double alpha,double *ts,double *t)
{
  *t=fabs(t_point(df,alpha/2));
  *ts=(b[nb]-b0)/sqrt(mse*xmatrix[nb][nb]);

  return (fabs(*ts)<*t) ? 1 : 0;
}

int y_one_test(char opt,int df,double y,double mse,
	       double alpha,double *ts,double *t)
{
  int i,j;
  double y0=0,h=0,temp[N]={0,};

  for (i=0;i<col;i++) {
    y0+=b[i]*xvector[i];
    for (j=0;j<col;j++) temp[i]+=xvector[j]*xmatrix[j][i];
  }
  for (i=0;i<col;i++) h+=temp[i]*xvector[i];

  *t=fabs(t_point(df,alpha));
  *ts=(y0-y)/sqrt(mse*h);

  if (opt=='1') return (*ts<=*t) ? 1 : 0;
  else return (*ts>=-*t) ? 1 : 0;
}

int y_two_test(int df,double y,double mse,double alpha,double *ts,double *t)
{
  int i,j;
  double y0=0,h=0,temp[N]={0,};

  for (i=0;i<col;i++) {
    y0+=b[i]*xvector[i];
    for (j=0;j<col;j++) temp[i]+=xvector[j]*xmatrix[j][i];
  }
  for (i=0;i<col;i++) h+=temp[i]*xvector[i];

  *t=fabs(t_point(df,alpha/2));
  *ts=(y0-y)/sqrt(mse*h);

  return (fabs(*ts)<*t) ? 1 : 0;
}

double nor_point(double p)

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

  q=p-.5;

  if (fabs(q)<=.425) {
    r=.425*.425-q*q;
    px=(((((((a[7]*r+a[6])*r+a[5])*r+a[4])*r+a[3])*r+a[2])*r+a[1])*r+a[0])/
       (((((((b[6]*r+b[5])*r+b[4])*r+b[3])*r+b[2])*r+b[1])*r+b[0])*r+1);
    return (q*px);
  }
  else {
    r=(q<0.) ? p : 1-p;
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

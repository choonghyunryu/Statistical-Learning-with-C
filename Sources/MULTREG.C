/* Program : MULTREG.C ver 0.1         */
/* Author  : Ryu choong hyun           */
/* Date    : 94.8.1.   	               */
/* Note    : MultiRegression Analysis  */

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
void est_y(int df,double mse,double alpha,double *l,double *u);
void est_bi(int i,int df,double mse,double alpha,double *l,double *u);

double nor_point(double p);
double f_q(int df1,int df2,double f);
double f_point(int df1,int df2,double p);
double t_q(int df,double t);
double t_point(int df,double point);
double t_sub(int df,double p);

int row,col;
double y[M],yhat[M],x[M][N],xmatrix[N][N],b[N],xvector[N],ymean;

void main(int argc,char *argv[])
{
  FILE *stream;
  int i,j,rank,dfr,dfe;
  double num,ysum=0,sst=0,sse=0,ssr,msr,mse,
	 f0,f_5,f_1,rsquare,sy_x,alpha,l,u;

  clrscr();

  if (argc<=1) {
    puts("Usage : multreg datafile");
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

  for (i=0;i<row;i++) {
    sst+=pow(y[i]-ymean,2);
    sse+=pow(y[i]-yhat[i],2);
  }
  ssr=sst-sse;
  rsquare=ssr/sst;
  dfr=rank-1;
  dfe=row-rank;
  msr=ssr/dfr;
  mse=sse/dfe;
  sy_x=sqrt(mse);
  f0=msr/mse;
  f_5=f_point(dfr,dfe,.05);
  f_1=f_point(dfr,dfe,.01);

  printf("\t  *** Multiple Regression Analysis ***\n\n");
  printf("\t << Regression Model >>\n");
  printf("\tY = %lf",b[0]);
  for (i=1;i<col;i++) printf(" + %lfX%d",b[i],i);
  printf(" + E\n\n");
  printf("\t << Precison >>\n");
  printf("\tStandard Error of Estimate       = %lf\n",sy_x);
  printf("\tCoefficient of Determination(R\xfd) = %lf\n\n",rsquare);

  printf("\t << Analysis of Variance Tabel >>\n");
  printf("Source      D.F.     S.S.       M.S.       F0        F(.05)   "
	 "F(.01)\n");
  printf("Regression   %d  %10.3f%10.3f %10.3f",dfr,ssr,msr,f0);
  if (f0>f_1) printf("**");
  if (f0>f_5 && f0<f_1) printf("*");
  printf("%10.3f%10.3f\n",f_5,f_1);
  printf("Error        %d  %10.3f%10.3f\n",dfe,sse,mse);
  printf("Total        %d  %10.3f\n\n",row-1,sst);

  printf("\t << Confidence Interval >>\n");
  printf("\tConfidence Level = ");
  scanf("%lf",&alpha);
  for (i=0;i<col;i++) {
    est_bi(i,dfe,mse,alpha,&l,&u);
    printf("\t\xe1%d Confidence Interval   : (%lf , %lf)\n",i,l,u);
  }

  for (i=0;i<col;i++)
    if (i==0) xvector[i]=1;
    else {
      printf("\tX%d               = ",i);
      scanf("%lf",&xvector[i]);
    }

  est_y(dfe,mse,alpha,&l,&u);
  printf("\tY  Confidence Interval   : (%lf , %lf)\n",l,u);

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

  for (i=0;i<col;i++)              /*deagak */
    for (j=i+1;j<row;j++)
      xx[j][i]-=xx[0][i]*xx[j][i]/xx[0][i];

  for (i=0;i<row;i++) {
    k=0;
    for (j=0;j<col;j++) if (xx[i][j]==0) k++;
    if (k==col) rank--;
  }
  return rank;
}

void est_y(int df,double mse,double alpha,double *l,double *u)
{
  int i,j;
  double t,y=0,h=0,temp[N]={0,};

  t=t_point(df,(1-alpha)/2);
  for (i=0;i<col;i++) {
    y+=b[i]*xvector[i];
    for (j=0;j<col;j++) temp[i]+=xvector[j]*xmatrix[j][i];
  }
  for (i=0;i<col;i++) h+=temp[i]*xvector[i];

  *l=y-fabs(t)*sqrt(mse*h);
  *u=y+fabs(t)*sqrt(mse*h);
}

void est_bi(int i,int df,double mse,double alpha,double *l,double *u)
{
  double t;

  t=t_point(df,(1-alpha)/2);

  *l=b[i]-fabs(t)*sqrt(mse*xmatrix[i][i]);
  *u=b[i]+fabs(t)*sqrt(mse*xmatrix[i][i]);
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

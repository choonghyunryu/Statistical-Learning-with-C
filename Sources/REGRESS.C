/* Program : REGRESS.C ver 0.1          */
/* Note    : Simple Regression Analysis */
/*           (Hypothesis Test)          */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define PI M_PI
#define FACTOR 50

double nor_point(double point);
double t_q(int df,double t);
double t_point(int df,double point);
double t_sub(int df,double p);

void get_mean_ss(int col,double *mean,double *ss);
double get_sxy(void);
double get_ssr(double b0,double b1);
double get_sse(double b0,double b1);
int b0_one_test(char opt,double b0,double b,double mse,
                double alpha,double *ts,double *t);
int b0_two_test(double b0,double b,double mse,
                double alpha,double *ts,double *t);
int b1_one_test(char opt,double b1,double b,double mse,
                double alpha,double *ts,double *t);
int b1_two_test(double b1,double b,double mse,
                double alpha,double *ts,double *t);
int y_one_test(char opt,double b0,double b1,double x,double y,double mse,
               double alpha,double *ts,double *t);
int y_two_test(double b0,double b1,double x,double y,double mse,
               double alpha,double *ts,double *t);
int rho_one_test(char opt,double r,double rho,
                 double alpha,double *ts,double *z);
int rho_two_test(double r,double rho,double alpha,double *ts,double *z);

int row;
double data[FACTOR][2],mean[2],ss[2];

void main(int argc,char *argv[])
{
  FILE *stream;
  char opt,job;
  int i,j,h0;
  double num,b0,b1,b,x,mu,s_xy,sst,ssr,sse,rsquare,r,rho,alpha,mse,ts,t;

  clrscr();

  if (argc<=1) {
    puts("Usage : stat datafile");
    exit(EXIT_FAILURE);
  }
  stream=fopen(argv[1],"rt");
  if (stream==NULL) printerrmsg("File not found !!");

  fscanf(stream,"%d\n",&row);

  for (i=0;i<row;i++)
    for (j=0;j<2;j++) {
      fscanf(stream,"%lf\n",&num);
      data[i][j]=num;
    }

  for (i=0;i<2;i++) get_mean_ss(i,&mean[i],&ss[i]);

  s_xy=get_sxy();
  b1=s_xy/ss[0];
  b0=mean[1]-b1*mean[0];
  sse=get_sse(b0,b1);
  mse=sse/(row-2);
  ssr=get_ssr(b0,b1);
  sst=ssr+sse;
  rsquare=ssr/sst;
  r=(b1>0) ? sqrt(rsquare) : -sqrt(rsquare);

  printf("\t *** Regression Analysis (Hypothesis Test) ***\n\n");
  printf("\t 1.\xe1%d Test  2.\xe1%d Test  3.\xe6y.x Test  4.rho Test\n\n",
          0,1);

  do
  switch (job=getch()) {
    case '1': printf("\t1.H0 : \xe1%d \xf3 \xe1  H1 : \xe1%d > \xe1\n",0,0);
              printf("\t2.H0 : \xe1%d \xf2 \xe1  H1 : \xe1%d < \xe1\n",0,0);
              printf("\t3.H0   :  \xe1%d   =   \xe1     H1   :   \xe1%d   \xd8 
\xe1\n\n",0,0);
              do
               opt=getch();
              while (!(opt=='1' || opt=='2' || opt=='3'));
              printf("\t\tb                  = ");
              scanf("%lf",&b);
              printf("\t\tSignificance Level = ");
              scanf("%lf",&alpha);

              if (opt=='3') h0=b0_two_test(b0,b,mse,alpha,&ts,&t);
              else  h0=b0_one_test(opt,b0,b,mse,alpha,&ts,&t);
              break;
    case '2': printf("\t1.H0 : \xe1%d \xf3 \xe1  H1 : \xe1%d > \xe1\n",1,1);
              printf("\t2.H0 : \xe1%d \xf2 \xe1  H1 : \xe1%d < \xe1\n",1,1);
              printf("\t3.H0   :  \xe1%d   =   \xe1     H1   :   \xe1%d   \xd8 
\xe1\n\n",1,1);
              do
               opt=getch();
              while (!(opt=='1' || opt=='2' || opt=='3'));
              printf("\t\tb                  = ");
              scanf("%lf",&b);
              printf("\t\tSignificance Level = ");
              scanf("%lf",&alpha);

              if (opt=='3') h0=b1_two_test(b1,b,mse,alpha,&ts,&t);
              else  h0=b1_one_test(opt,b1,b,mse,alpha,&ts,&t);
              break;
    case '3': printf("\t1.H0 : \xe6y.x \xf3 \xe6  H1 : \xe6y.x > \xe1\n");
              printf("\t2.H0 : \xe6y.x \xf2 \xe6  H1 : \xe6y.x < \xe1\n");
              printf("\t3.H0 : \xe6y.x = \xe6  H1 : \xe6y.x \xd8 \xe1\n\n");
              do
               opt=getch();
              while (!(opt=='1' || opt=='2' || opt=='3'));
              printf("\t\tX                  = ");
              scanf("%lf",&x);
              printf("\t\t\xe6                  = ");
              scanf("%lf",&mu);
              printf("\t\tSignificance Level = ");
              scanf("%lf",&alpha);

              if (opt=='3') h0=y_two_test(b0,b1,x,mu,mse,alpha,&ts,&t);
              else  h0=y_one_test(opt,b0,b1,x,mu,mse,alpha,&ts,&t);
              break;
    case '4': printf("\t1.H0 : rho \xf3 rho0  H1 : rho > rho0\n");
              printf("\t2.H0 : rho \xf2 rho0  H1 : rho < rho0\n");
              printf("\t3.H0 : rho = rho0  H1 : rho \xd8 rho0\n\n");
              do
               opt=getch();
              while (!(opt=='1' || opt=='2' || opt=='3'));
              printf("\t\tRho                = ");
              scanf("%lf",&rho);
              printf("\t\tSignificance Level = ");
              scanf("%lf",&alpha);

              if (opt=='3') h0=rho_two_test(r,rho,alpha,&ts,&t);
              else  h0=rho_one_test(opt,r,rho,alpha,&ts,&t);
              break;
    default : break;
  } while (!(job=='1' || job=='2' || job=='3' || job=='4'));

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

void get_mean_ss(int col,double *mean,double *ss)
{
  int i;
  double sum=0,sos=0;

  for (i=0;i<row;i++) sum+=data[i][col];

  *mean=sum/row;

  for (i=0;i<row;i++) sos+=pow(data[i][col]-*mean,2);

  *ss=sos;
}

double get_sxy(void)
{
  int i;
  double sum=0;

  for (i=0;i<row;i++) sum+=(data[i][0]-mean[0])*(data[i][1]-mean[1]);

  return sum;
}

double get_ssr(double b0,double b1)
{
  int i;
  double sum=0;

  for (i=0;i<row;i++) sum+=pow((b0+b1*data[i][0])-mean[1],2);

  return sum;
}

double get_sse(double b0,double b1)
{
  int i;
  double sum=0;

  for (i=0;i<row;i++) sum+=pow(data[i][1]-(b0+b1*data[i][0]),2);

  return sum;
}

int b0_one_test(char opt,double b0,double b,double mse,
                double alpha,double *ts,double *t)
{
  *t=fabs(t_point(row-2,alpha));
  *ts=(b0-b)/sqrt(mse*(1./row+pow(mean[0],2)/ss[0]));

  if (opt=='1') return (*ts<=*t) ? 1 : 0;
  else return (*ts>=-*t) ? 1 : 0;
}

int b0_two_test(double b0,double b,double mse,
                double alpha,double *ts,double *t)
{
  *t=fabs(t_point(row-2,alpha/2));
  *ts=(b0-b)/sqrt(mse*(1./row+pow(mean[0],2)/ss[0]));

  return (fabs(*ts)<*t) ? 1 : 0;
}

int b1_one_test(char opt,double b1,double b,double mse,
                double alpha,double *ts,double *t)
{
  *t=fabs(t_point(row-2,alpha));
  *ts=(b1-b)/sqrt(mse/ss[0]);

  if (opt=='1') return (*ts<=*t) ? 1 : 0;
  else return (*ts>=-*t) ? 1 : 0;
}

int b1_two_test(double b1,double b,double mse,
                double alpha,double *ts,double *t)
{
  *t=fabs(t_point(row-2,alpha));
  *ts=(b1-b)/sqrt(mse/ss[0]);

  return (fabs(*ts)<*t) ? 1 : 0;
}

int y_one_test(char opt,double b0,double b1,double x,double y,double mse,
               double alpha,double *ts,double *t)
{
  double y0;

  y0=b0+b1*x;

  *t=fabs(t_point(row-2,alpha));
  *ts=(y0-y)/sqrt(mse*(1./row+pow(x-mean[0],2)/ss[0]));

  if (opt=='1') return (*ts<=*t) ? 1 : 0;
  else return (*ts>=-*t) ? 1 : 0;
}

int y_two_test(double b0,double b1,double x,double y,double mse,
               double alpha,double *ts,double *t)
{
  double y0;

  y0=b0+b1*x;

  *t=fabs(t_point(row-2,alpha/2));
  *ts=(y0-y)/sqrt(mse*(1./row+pow(x-mean[0],2)/ss[0]));

  return (fabs(*ts)<*t) ? 1 : 0;
}

int rho_one_test(char opt,double r,double rho,
                 double alpha,double *ts,double *z)
{
  double zz,ez;

  zz=log((1+r)/(1-r))/2;
  ez=log((1+rho)/(1-rho))/2;

  *z=fabs(nor_point(alpha));
  *ts=(zz-ez)/sqrt(1./(row-3));

  if (opt=='1') return (*ts<=*z) ? 1 : 0;
  else return (*ts>=-*z) ? 1 : 0;
}

int rho_two_test(double r,double rho,double alpha,double *ts,double *z)
{
  double zz,ez;

  zz=log((1+r)/(1-r))/2;
  ez=log((1+rho)/(1-rho))/2;

  *z=fabs(nor_point(alpha/2));
  *ts=(zz-ez)/sqrt(1./(row-3));

  return (fabs(*ts)<*t) ? 1 : 0;
}


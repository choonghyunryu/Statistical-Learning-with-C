/* Program : DIAGON.C ver 0.1          */
/* Author  : Ryu choong hyun           */
/* Date    : 94.7.30.	              */
/* Note    : Regression Diagonistics    */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define PI M_PI
#define FACTOR 50

void get_mean_ss(int col,double *mean,double *ss);
double get_sxy(void);
double get_sse(double b0,double b1);

int row;
double data[FACTOR][2],mean[2],ss[2];

void main(int argc,char *argv[])
{
  FILE *stream;
  int i,j;
  double num,b0,b1,s_xy,sse,mse,e[FACTOR],h[FACTOR],r[FACTOR],rr[FACTOR],
         dff[FACTOR],d[FACTOR],m[FACTOR],ap[FACTOR],cov[FACTOR],
         fav[FACTOR],k=0,s;

  clrscr();

  if (argc<=1) {
    puts("Usage : diagon datafile");
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

  for (i=0;i<row;i++) {
    e[i]=data[i][1]-(b0+b1*data[i][0]);
    h[i]=1./row+pow(data[i][0]-mean[0],2)/ss[0];
    r[i]=e[i]/(sqrt(mse)*sqrt(1-h[i]));
    k+=h[i];
  }
  k--;

  for (i=0;i<row;i++) {
    s=sqrt(((row-k-1)*mse-pow(e[i],2)/(1-h[i]))/(row-k-2));
    rr[i]=e[i]/(s*sqrt(1-h[i]));
  }

  for (i=0;i<row;i++) {
    dff[i]=sqrt(h[i]/(1-h[i]))*rr[i];
    d[i]=h[i]/((k+1)*(1-h[i]))*pow(r[i],2);
    m[i]=(row-1)*(h[i]-1./row);
    ap[i]=1-h[i]-pow(e[i],2)/((row-k-1)*(mse));
    cov[i]=1/(pow(1+(pow(rr[i],2)-1)/(row-k-1),k+1)*(1-h[i]));
    fav[i]=pow(e[i],2)/(pow(rr[i],2)*pow(1-h[i],2)*mse);
  }
  printf("\t     *** Regression Analysis Diagonistics ***\n\n");

  printf("\t  NO    Residual      hii    Standard R. Student R.  DFFITS\n");
  for (i=0;i<row;i++) printf("\t%3d  %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                             i+1,e[i],h[i],r[i],rr[i],dff[i]);
  getch();
  printf("\n\n");
  printf("\t  NO    Cook's D    M-diff    A-P Stat   COVRATIO   FVARATIO\n");
  for (i=0;i<row;i++) printf("\t%3d  %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                             i+1,d[i],m[i],ap[i],cov[i],fav[i]);
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

double get_sse(double b0,double b1)
{
  int i;
  double sum=0;

  for (i=0;i<row;i++) sum+=pow(data[i][1]-(b0+b1*data[i][0]),2);

  return sum;
}


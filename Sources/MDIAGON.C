/* Program : MDIAGON.C ver 0.1         */
/* Author  : Ryu choong hyun           */
/* Date    : 94.8.1.   	               */
/* Note    : MultiRegression Analysis  */
/*           Diagonistics              */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define M 50
#define N 6

void get_xmatrix(void);
void matrix_inverse(void);
void get_bvector(void);
void get_yhat(void);
int get_rank(void);
double get_hii(int i);

int row,col;
double y[M],yhat[M],x[M][N],xmatrix[N][N],b[N],xvector[N],ymean;

void main(int argc,char *argv[])
{
  FILE *stream;
  int i,j,rank,dfe;
  double num,ysum=0,sse=0,mse,e[M],h[M],r[M],rr[M],dff[M],d[M],m[M],
         ap[M],cov[M],fav[M];

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

  for (i=0;i<row;i++) sse+=pow(y[i]-yhat[i],2);

  dfe=row-rank;
  mse=sse/dfe;

  for (i=0;i<row;i++) {
    e[i]=y[i]-yhat[i];
    h[i]=get_hii(i);
    r[i]=e[i]/(sqrt(mse)*sqrt(1-h[i]));
    rr[i]=e[i]/(sqrt((dfe*mse-pow(e[i],2)/(1-h[i]))/(dfe-1))*sqrt(1-h[i]));
    dff[i]=sqrt(h[i]/(1-h[i]))*rr[i];
    d[i]=h[i]/(rank*(1-h[i]))*pow(r[i],2);
    m[i]=(row-1)*(h[i]-1./row);
    ap[i]=1-h[i]-pow(e[i],2)/(dfe*(mse));
    cov[i]=1/(pow(1+(pow(rr[i],2)-1)/dfe,rank)*(1-h[i]));
    fav[i]=pow(e[i],2)/(pow(rr[i],2)*pow(1-h[i],2)*mse);
  }
  printf("\t*** Multple Regression Analysis Diagonistics ***\n\n");

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

double get_hii(int i)
{
  int j,k;

  double sum,temp[N];

  for (j=0;j<col;j++) {
    sum=0;
    for (k=0;k<col;k++) sum+=x[i][k]*xmatrix[k][j];
    temp[j]=sum;
  }

  sum=0;
  for (j=0;j<col;j++) sum+=temp[j]*x[i][j];

  return sum;
}

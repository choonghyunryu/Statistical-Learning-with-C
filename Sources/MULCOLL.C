/* Program : MULCOLL.C ver 0.1         */
/* Author  : Ryu choong hyun           */
/* Date    : 94.8.7.   	               */
/* Note    : Multicollinearity Check   */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define PI M_PI
#define M 50
#define N 5+1

void get_xstar_matrix(void);
void xstar_matrix_inverse(void);
void get_xstar(void);
double get_vif1(int n);

int row,col;
double xstar[M][N],xstar_matrix[N][N],x[M][N];

void main(int argc,char *argv[])
{
  FILE *stream;
  int i,j;
  double num,vif,rank;

  clrscr();

  if (argc<=1) {
    puts("Usage : mulcoll datafile");
    exit(EXIT_FAILURE);
  }

  stream=fopen(argv[1],"rt");
  if (stream==NULL) printerrmsg("File not found !!");

  fscanf(stream,"%d\n",&row);
  fscanf(stream,"%d\n",&col);

  for (i=0;i<row;i++) {
    for (j=0;j<col;j++) {
    fscanf(stream,"%lf\n",&num);
      if (j==0) x[i][j]=1;
      else x[i][j]=num;
    }
  }

  printf("\t\t *** Variance Inflation ***\n\n");
  printf("\t\t     Variable       VIF\n");

  for (i=0;i<col-1;i++) {
    vif=get_vif1(i);
    printf("\t\t\tX%d        %lf\n",i+1,vif);
  }

  getch();
}

void get_xstar_matrix(void)
{
  int i,j,k;

  double sum;

  for (i=0;i<col;i++)
    for (j=0;j<col;j++) {
      sum=0;
      for (k=0;k<row;k++) sum+=xstar[k][i]*xstar[k][j];
      xstar_matrix[i][j]=sum;
    }
}

void xstar_matrix_inverse(void)
{
  int i,j,k;
  double temp,u;

  for (k=0;k<col;k++) {
    temp=xstar_matrix[k][k];
    for (i=0;i<col;i++) xstar_matrix[k][i]/=temp;
    xstar_matrix[k][k]=1/temp;
    for (j=0;j<col;j++)
      if (j!=k) {
	u=xstar_matrix[j][k];
	for (i=0;i<col;i++)
	  if (i!=k) xstar_matrix[j][i]-=xstar_matrix[k][i]*u;
	  else xstar_matrix[j][i]=-u/temp;
      }
  }
}

void get_xstar(void)
{
  int i,j,k;
  double sxx[N-1]={0,},xsum[N-1]={0,},xmean[N-1];

  for (i=0;i<row;i++)
    for (j=0;j<col;j++) {
      xstar[i][j]=x[i][j];
      xsum[j]+=x[i][j+1];
    }

  for (i=1;i<col;i++) xmean[i-1]=xsum[i-1]/row;

  for (i=0;i<row;i++)
    for (j=1;j<col;j++) sxx[j-1]+=pow(xstar[i][j]-xmean[j-1],2);

  for (i=0;i<row;i++)
    for (j=1;j<col;j++)
      xstar[i][j]=(xstar[i][j]-xmean[j-1])/sqrt(sxx[j-1]);
}

double get_vif1(int n)
{
  int i,j;

  get_xstar();
  get_xstar_matrix();
  xstar_matrix_inverse();

  return xstar_matrix[n+1][n+1];
}


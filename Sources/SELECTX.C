/* Program : SELECTX.C ver 0.1      */
/* Author  : Ryu choong hyun        */
/* Date    : 94.8.7.   	            */
/* Note    : Selction of Variables  */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define PI M_PI
#define M 15
#define N 6

void get_xmatrix(int p,int n);
void matrix_inverse(int p);
double ssep(int p);
double get_sst(void);
double get_variance(void);
void criteron(int p);

int row,col,kk[30],kk1[30],m;
double y[M],x[M][N],xx[M][N],xmatrix[M][M],sst,var;

void main(int argc,char *argv[])
{
  FILE *stream;
  int i,j,p;
  double num,sse,mse,rsquare,adr,cp;

  clrscr();

  if (argc<=1) {
    puts("Usage : selectx datafile");
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
	x[i][j]=1;
      }
      else x[i][j]=num;
    }

  sst=get_sst();
  var=get_variance();

  printf("       ** Selection of Variable of Criterion **\n\n");
  printf("In   Rsquare   Adj-Rsquare     CP         MSE      "
	 "Varables in Model\n");

  for (i=1;i<col;i++) criteron(i);

  getch();
}

void get_xmatrix(int p,int n)
{
  int i,j,k;

  k=(p==1) ? p-1 : p;

  for (i=0;i<row;i++)
    for (j=k;j<p+1;j++) {
      if (j==0) xx[i][j]=1;
      else xx[i][j]=x[i][n];
    }
}

void matrix_inverse(int p)
{
  int i,j,k;
  double temp,u;

  for (k=0;k<p+1;k++) {
    temp=xmatrix[k][k];
    for (i=0;i<p+1;i++) xmatrix[k][i]/=temp;
    xmatrix[k][k]=1/temp;
    for (j=0;j<p+1;j++)
      if (j!=k) {
	u=xmatrix[j][k];
	for (i=0;i<p+1;i++)
	  if (i!=k) xmatrix[j][i]-=xmatrix[k][i]*u;
	  else xmatrix[j][i]=-u/temp;
      }
  }
}

double ssep(int p)
{
  int i,j,k;

  double sum,temp[M][M],u;

  for (i=0;i<p+1;i++)
    for (j=0;j<p+1;j++) {
      sum=0;
      for (k=0;k<row;k++) sum+=xx[k][i]*xx[k][j];
      xmatrix[i][j]=sum;
    }

  matrix_inverse(p);

  for (i=0;i<row;i++)
    for (j=0;j<p+1;j++) {
      sum=0;
      for (k=0;k<p+1;k++) sum+=xx[i][k]*xmatrix[k][j];
      temp[i][j]=sum;
    }

  for (i=0;i<row;i++)
    for (j=0;j<row;j++) {
      sum=0;
      for (k=0;k<p+1;k++) sum+=temp[i][k]*xx[j][k];
      xmatrix[i][j]=sum;
    }

  for (i=0;i<row;i++)
    for (j=0;j<row;j++)
      if (i==j) xmatrix[i][j]=1-xmatrix[i][j];
      else xmatrix[i][j]*=-1;

  for (i=0;i<row;i++) {
    sum=0;
    for (j=0;j<row;j++) {
      sum+=y[j]*xmatrix[j][i];
    }
    xmatrix[i][0]=sum;
    }

  sum=0;
  for (i=0;i<row;i++) sum+=y[i]*xmatrix[i][0];

  return sum;
}

double get_sst(void)
{
  int i;
  double sum=0,ymean,sst=0;

  for (i=0;i<row;i++) sum+=y[i];
  ymean=sum/row;

  for (i=0;i<row;i++) sst+=pow(y[i]-ymean,2);

  return sst;
}

double get_variance(void)
{
  int i,j,k;

  double sum,xtemp[M][M],xtemp1[M][M],temp,u;

  for (i=0;i<col;i++)
    for (j=0;j<col;j++) {
      sum=0;
      for (k=0;k<row;k++) sum+=x[k][i]*x[k][j];
      xtemp[i][j]=sum;
    }

  for (k=0;k<col;k++) {
    temp=xtemp[k][k];
    for (i=0;i<col;i++) xtemp[k][i]/=temp;
    xtemp[k][k]=1/temp;
    for (j=0;j<col;j++)
      if (j!=k) {
	u=xtemp[j][k];
	for (i=0;i<col;i++)
	  if (i!=k) xtemp[j][i]-=xtemp[k][i]*u;
	  else xtemp[j][i]=-u/temp;
      }
  }

  for (i=0;i<row;i++)
    for (j=0;j<col;j++) {
      sum=0;
      for (k=0;k<col;k++) sum+=x[i][k]*xtemp[k][j];
      xtemp1[i][j]=sum;
    }

  for (i=0;i<row;i++)
    for (j=0;j<row;j++) {
      sum=0;
      for (k=0;k<col;k++) sum+=xtemp1[i][k]*x[j][k];
      xtemp[i][j]=sum;
    }

  for (i=0;i<row;i++)
    for (j=0;j<row;j++)
      if (i==j) xtemp[i][j]=1-xtemp[i][j];
      else xtemp[i][j]*=-1;

  for (i=0;i<row;i++) {
    sum=0;
    for (j=0;j<row;j++) {
      sum+=y[j]*xtemp[j][i];
    }
    xtemp[i][0]=sum;
    }

  sum=0;
  for (i=0;i<row;i++) sum+=y[i]*xtemp[i][0];

  return sum/(row-col);
}

void criteron(int p)
{
  int i,j,k,l;
  double sse,mse,rsquare,adr,cp;

  if (p==1) {
    for (i=1;i<col;i++) {
      get_xmatrix(p,i);
      sse=ssep(p);
      mse=sse/(row-p-1);
      rsquare=1-((row-p-1)/sst)*mse;
      adr=1-(row-1)*(1-rsquare)/(row-p-1);
      cp=sse/var+2*(p+1)-row;
      printf("%d %10.5f  %10.5f  %10.5f  %10.5f     X%d\n",
	      p,rsquare,adr,cp,mse,i);
    }
  }
  else if (p==2) {
    for (i=1;i<col;i++) {
      get_xmatrix(1,i);
      for (j=i+1;j<col;j++) {
	get_xmatrix(p,j);
	kk[m++]=i;
	kk[m++]=j;
	sse=ssep(p);
	mse=sse/(row-p-1);
	rsquare=1-((row-p-1)/sst)*mse;
	adr=1-(row-1)*(1-rsquare)/(row-p-1);
	cp=sse/var+2*(p+1)-row;
	printf("%d %10.5f  %10.5f  %10.5f  %10.5f     X%d",
	       p,rsquare,adr,cp,mse,i);
	printf(" X%d\n",j);
      }
    }
  }
  else {
    m=0;
    for (i=1;i<col;i++) {
      for (j=0;j<30;j+=p-1) {
	l=1;
	if (kk[j]>i) {
	  get_xmatrix(l++,i);
	  kk1[m++]=i;
	  for (k=j;k<j+p-1;k++) {
	    get_xmatrix(l++,kk[k]);
	    kk1[m++]=kk[k];
	  }
	  sse=ssep(p);
	  mse=sse/(row-p-1);
	  rsquare=1-((row-p-1)/sst)*mse;
	  adr=1-(row-1)*(1-rsquare)/(row-p-1);
	  cp=sse/var+2*(p+1)-row;
	  printf("%d %10.5f  %10.5f  %10.5f  %10.5f     X%d",
		 p,rsquare,adr,cp,mse,i);
	  for (k=j;k<j+p-1;k++) printf(" X%d",kk[k]);
	  printf("\n");
	  if (p==col-1) return;
	}
      }
    }

    for (i=0;i<30;i++) kk[i]=kk1[i];
  }
}

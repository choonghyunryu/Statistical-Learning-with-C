/* Program : NORND.C ver 0.1      */
/* Author  : Ryu choong hyun      */
/* Date    : 94.8.21.		  */
/* Note    : Normal Random Number */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <time.h>

#define PI M_PI
#define N 1000

double rnd(void);
double nornd1(void);
double nornd2(void);
double nornd3(void);
void probability(void);

int n;
double x[N],pro[11],mu,sigma;

void main(void)
{
  char opt;
  int i;
  double temp,sigma2,sum=0,sum_s=0,mean,dev;

  clrscr();

  randomize();

  printf("\t\t** Normal Random Numbers **\n\n");
  printf("     1.Central Limit Theory  2.Box-Muller  3.Marsaglia\n\n");

  do
  switch (opt=getch()) {
    case '1' :
    case '2' :
    case '3' : break;
    default  : break;
  } while (!(opt=='1' || opt=='2' || opt=='3'));

  printf("\t\tNumber of Random-Number = ");
  scanf("%d",&n);
  printf("\t\tNormal(\xe6,\xe5\xfd) --->\xe6      = ");
  scanf("%lf",&mu);
  printf("\t\tNormal(\xe6,\xe5\xfd) --->\xe5\xfd     = ");
  scanf("%lf",&sigma2);

  for (i=0;i<n;i++) {
    switch (opt) {
      case '1': temp=nornd1();
		break;
      case '2': temp=nornd2();
		break;
      case '3': temp=nornd3();
		break;
    }

    sigma=sqrt(sigma2);
    x[i]=sigma*temp+mu;
    sum+=x[i];
  }

  mean=sum/n;

  for (i=0;i<n;i++) {
    sum_s+=pow(x[i]-mean,2);
  }

  dev=sqrt(sum_s/n);

  printf("\n\t\t\t << Random-Number >>\n");
  for (i=0;i<n;i++) printf("\t\t\t    %8.5f\n",x[i]);
  printf("\n\t\t    Mean  = %8.5f(%lf)\n",mean,mu);
  printf("\t\t    Var   = %8.5f(%lf)\n",pow(dev,2),sigma2);
  printf("\t\t    S.D.  = %8.5f(%lf)\n\n",dev,sigma);

  probability();
  printf("\t\t\t  Interval       P\n");
  for (i=0;i<10;i++) {
    if (i<5) printf("\t\t\t[\xe6-%d\xe5 ,\xe6-%d\xe5)   %5.3f\n",
		    5-i,4-i,pro[i]);
    else printf("\t\t\t[\xe6+%d\xe5 ,\xe6+%d\xe5)   %5.3f\n",i-5,i-4,pro[i]);
  }

  getch();
}

double rnd(void)
{
  return (1.0/(RAND_MAX+1.0))*(rand()+0.5);
}

double nornd1(void)
{
  int i;
  double sum=0;

  for (i=0;i<12;i++) sum+=rnd();

  return sum-6;
}

double nornd2(void)
{
  static int i=0;
  static double u,temp;

  if (i==0) {
    i++;
    temp=sqrt(-2*log(rnd()));
    u=2*PI*rnd();
    return temp*cos(u);
  }
  else {
    i=0;
    return temp*sin(u);
  }
}

double nornd3(void)
{
  static int i=0;
  static double v1,v2,s;

  if (i==0) {
    i++;
    do {
      v1=2*rnd()-1;
      v2=2*rnd()-1;
      s=pow(v1,2)+pow(v2,2);
    } while (s>1 || s==0);
    s=sqrt(-2*log(s)/s);
    return v1*s;
  }
  else {
    i=0;
    return v2*s;
  }
}

void probability(void)
{
  int i,j,temp[10]={0,};

  for (i=0;i<n;i++) {
    j=((x[i]-mu)/sigma)+5;
    temp[j]++;
  }

  for (i=0;i<10;i++) pro[i]=(double)temp[i]/n;
}

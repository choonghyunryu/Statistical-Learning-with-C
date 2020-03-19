/* Program : DISCRETE.C ver 0.2               */
/* Author  : Ryu choong hyun                  */
/* Date    : 94.3.30.		              */
/* Note    : Binimial Distribution &          */
/*	     Poisson Distribution &           */
/*           Geometric Distribution           */
/*	     Hypergeometric Distribtion &     */
/*           Negative Binomial Distribution & */
/*           Poly & Etc.                      */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <conio.h>
#include <math.h>
#include <graphics.h>

#define round(x) ((x>0) ? floor(x+.5) : ceil(x-.5))
#define MAX 1000

void go_binomial(void);
void go_poisson(void);
void go_geo(void);
void go_hyper(void);
void go_negbin(void);

double fac(int x);
double com(int n,int r);
double rnd(double x);
void init_array(void);

void binomial(int n,double p);
void poisson(unsigned long n,double mean);
void geometric(int n,double p);
void hyper(int n,int a,int s);
void negbinomial(int n,double p);

int gprint(int x,int y,char *fmt, ...);
void init_graph(void);
void board(void);
void ppoly(char *str);
void cpoly(char *str);

char go;
unsigned long n;
int start,last,maxy;
double pdf[MAX],cdf[MAX];
double *pp,*pc;
int x0,y0,x_range,y_range;

void main(void)
{
  pp=pdf,pc=cdf;
  while(1) {
    clrscr();
    init_array();
    printf("\t\t        ** Discrete Probability Distribution **\n\n\n");
    printf("\t\t    1.  Binomial Probability Distribution\n\n");
    printf("\t\t    2.  Poisson Probability Distribution\n\n");
    printf("\t\t    3.  Geometric Probability Distribution\n\n");
    printf("\t\t    4.  Hypergeometric Probability Distribution\n\n");
    printf("\t\t    5.  Negative Binomial Probability Distribution\n\n");
    printf("\t\t    6.  Quit");
    switch (go=getch()) {
      case '1' : go_binomial();
		 break;
      case '2' : go_poisson();
		 break;
      case '3' : go_geo();
		 break;
      case '4' : go_hyper();
		 break;
      case '5' : go_negbin();
		 break;
      case '6' : exit(0);
		 break;
      default  : break;
    }
  }
}

void go_binomial(void)
{
  int i,j=0;
  float p;
  double mean,var,dev,skewness,kurtosis;
  char str[]="Binomial Distribution";
  char str1[]="Cumul-Binomial Distribution";

  clrscr();
  printf("\n\t ** Binomial Distribution **\n");
  printf("\n\t    N = ");
  scanf("%d",&n);
  printf("\n\t    P = ");
  scanf("%f",&p);
  mean=n*p;
  var=n*p*(1-p);
  dev=sqrt(var);
  skewness=(1-p)-p/dev;
  kurtosis=3+(1-6*(p*(1-p)))/dev;
  clrscr();
  printf("\n\t\t ** Binomial Distribution **\n\n");
  printf("\t\t     N         = %8d\n",n);
  printf("\t\t     P         = %8.6f\n",p);
  printf("\t\t     Mean      = %f\n",mean);
  printf("\t\t     Variance  = %f\n",var);
  printf("\t\t     Deviation = %f\n",dev);
  printf("\t\t     Skewness  = %f\n",skewness);
  printf("\t\t     Kurtosis  = %f\n\n",kurtosis);
  getch();

  binomial(n,p);

  clrscr();
  printf("\t\t     X      PDF        CDF\n");
  if (start==-1)
    for (i=start+1;i<=last;i++,j++) {
      printf("\t\t  %4d   %f    %f\n",i,*(pp+i),*(pc+i));
      if (j % 20==19) getch();
  }
  else
    for (i=start;i<=last;i++,j++) {
      printf("\t\t  %4d   %f    %f\n",i,*(pp+i),*(pc+i));
      if (j % 20==19) getch();
    }

  getch();

  init_graph();
  ppoly(str);
  cpoly(str1);
  closegraph();
}

void go_poisson(void)
{
  int i,j=0;
  char c;
  float p,mean,var,dev,skewness,kurtosis;
  char str[]="Poisson Distribution";
  char str1[]="Cumul-Poisson Distribution";

  clrscr();
  printf("\n\t ** Poisson Distribution **\n");
  printf("\n\t What is Data ?  (1.N,P  2.MEAN)\n");

  if (getch()=='1') {
    printf("\n\t    N = ");
    scanf("%ld",&n);
    printf("\n\t    P = ");
    scanf("%f",&p);
    mean=n*p;
    var=mean;
    dev=sqrt(var);
    skewness=1/dev;
    kurtosis=3+1/var;
    clrscr();
    printf("\n\t\t ** Poisson Distribution **\n\n");
    printf("\t\t     N         = %ld\n",n);
    printf("\t\t     P         = %8.6f\n",p);
    printf("\t\t     Mean      = %f\n",mean);
    printf("\t\t     Variance  = %f\n",var);
    printf("\t\t     Deviation = %f\n",dev);
    printf("\t\t     Skewness  = %f\n",skewness);
    printf("\t\t     Kurtosis  = %f\n\n",kurtosis);
    getch();
    poisson(n,mean);
  }
  else {
    printf("\n\t  Mean =  ");
    scanf("%f",&mean);
    var=mean;
    dev=sqrt(var);
    skewness=1/dev;
    kurtosis=3+1/var;
    clrscr();
    printf("\n\t\t ** Poisson Distribution **\n\n");
    printf("\t\t     Mean      = %f\n",mean);
    printf("\t\t     Variance  = %f\n",var);
    printf("\t\t     Deviation = %f\n",dev);
    printf("\t\t     Skewness  = %f\n",skewness);
    printf("\t\t     Kurtosis  = %f\n\n",kurtosis);
    getch();
    poisson(1000,mean);
  }

  clrscr();
  printf("\t\t     X      PDF        CDF\n");
  for (i=start+1;i<=last;i++,j++) {
    printf("\t\t  %4d   %f    %f\n",i,*(pp+i),*(pc+i));
    if (j % 20==19) getch();
  }
  getch();

  init_graph();
  ppoly(str);
  cpoly(str1);
  closegraph();
}

void go_geo(void)
{
  int i;
  double p,mean,var,dev;
  char str[]="Geometric Distribution";
  char str1[]="Cumul-Geometric Distribution";

  clrscr();
  printf("\n\t ** Geometic Distribution **\n");
  printf("\n\t    First Succes N = ");
  scanf("%d",&n);
  printf("\n\t    Probability    = ");
  scanf("%lf",&p);
  mean=1/p;
  var=(1-p)/p*p;
  dev=sqrt(var);
  clrscr();
  printf("\n\t\t ** Geometric Distribution **\n\n");
  printf("\t\t     First Succes N         = %8d\n",n);
  printf("\t\t     Probability            = %8.6f\n",p);
  printf("\t\t     Mean                   = %lf\n",mean);
  printf("\t\t     Variance               = %lf\n",var);
  printf("\t\t     Standard Deviation     = %lf\n\n",dev);
  getch();

  geometric(n,p);

  clrscr();
  printf("\t\t     X      PDF        CDF\n");
  for (i=1;i<=last;i++) {
    printf("\t\t  %4d   %lf    %lf\n",i,pdf[i],cdf[i]);
    if (i % 21==20) getch();
  }

  getch();

  init_graph();
  ppoly(str);
  cpoly(str1);
  closegraph();
}

void go_hyper(void)
{
  int i,j=0,n,a,s;
  double mean,var,dev;
  char str[]="Hypergeometric Distribution";
  char str1[]="Cumul-Hypergeometic Distribution";

  clrscr();
  printf("\n\t ** Hypergeometric Distribution **\n");
  printf("\n\t    N      = ");
  scanf("%d",&n);
  printf("\n\t    A-N    = ");
  scanf("%d",&a);
  printf("\n\t    Sample = ");
  scanf("%d",&s);

  mean=(double)s*a/n;
  var=mean*(1-a/n)*(n-s)/(n-1);
  dev=sqrt(var);
  clrscr();
  printf("\n\t\t ** Hypergeometric Distribution **\n\n");
  printf("\t\t         N         = %4d\n",n);
  printf("\t\t         A-N       = %4d\n",a);
  printf("\t\t         Sample    = %4d\n",s);
  printf("\t\t         Mean      = %f\n",mean);
  printf("\t\t         Variance  = %f\n",var);
  printf("\t\t         Deviation = %f\n\n",dev);
  getch();

  hyper(n,a,s);

  clrscr();
  printf("\t\t     X      PDF        CDF\n");
  for (i=start+1;i<=last;i++,j++) {
    printf("\t\t  %4d   %f    %f\n",i,*(pp+i),*(pc+i));
    if (j % 20==19) getch();
  }
  getch();

  init_graph();
  ppoly(str);
  cpoly(str1);
  closegraph();
}

void go_negbin(void)
{
  int i;
  double p,mean,var,dev;
  char str[]="Negative Binomial Distribution";
  char str1[]="Cumul-Negative Binomial Distribution";

  clrscr();
  printf("\n\t ** Negative Binomial Distribution **\n");
  printf("\n\t    Sucess K = ");
  scanf("%d",&n);
  printf("\n\t    P        = ");
  scanf("%lf",&p);
  mean=n/p;
  var=n*(1-p)/p*p;
  dev=sqrt(var);
  clrscr();
  printf("\n\t\t ** Negative Binomial Distribution **\n\n");
  printf("\t\t     Sucess K  = %8d\n",n);
  printf("\t\t     P         = %8.6f\n",p);
  printf("\t\t     Mean      = %lf\n",mean);
  printf("\t\t     Variance  = %lf\n",var);
  printf("\t\t     Deviation = %lf\n\n",dev);
  getch();

  negbinomial(n,p);

  clrscr();
  printf("\t\t     X      PDF        CDF\n");
  for (i=n;i<=last;i++) {
    printf("\t\t  %4d   %lf    %lf\n",i,*(pp+i),*(pc+i));
    if (i % 21==20) getch();
  }
  getch();

  init_graph();
  ppoly(str);
  cpoly(str1);
  closegraph();
}

double fac(int x)
{
  return (x==0 ? 1 : x*fac(x-1));
}

double com(int n,int r)
{
  int i,k=n-r;
  long double npr=1,ir=1;

  if (n<r) return 0;
  for (i=n;i>r;i--) {
    npr*=i;
    ir*=k--;
  }

  return (npr/ir);
}


double rnd(double x)
{
  return round(x*1000000)/1000000;
}

void init_array(void)
{
  int i;
  for (i=0;i<=MAX;i++) {
    *(pp+i)=0;
    *(pc+i)=0;
  }
}

void binomial(int n,double p)
{
  int i;
  double max=0,sum=0;

  for (i=0;i<=n;i++) {
    *(pp+i)=com(n,i)*pow(p,i)*pow(1-p,n-i);
    sum+=*(pp+i);
    *(pc+i)=sum;
    if (rnd(*(pp+i))==0 && rnd(*(pc+i))==1) {
      if (rnd(*(pc+(i-1)))==1) last=i-1;
      else last=i;
      break;
    }
  }

  if (last==0) last=n;

  for (i=0;i<=last;i++)
    if (*(pp+i)>max) {
      max=*(pp+i);
      maxy=i;
    }
  for (i=0;i<maxy;i++)
    if (rnd(*(pp+i))>0) {
      start=i-1;
      break;
    }
}

void poisson(unsigned long n,double mean)
{
  int i;
  double max=0,sum=0;

  for (i=0;i<=n;i++) {
    *(pp+i)=exp(-mean)*pow(mean,i)/fac(i);
    sum+=*(pp+i);
    *(pc+i)=sum;
    if (rnd(*(pp+i))==0 && rnd(cdf[i])==1) {
      if (rnd(cdf[i-1])==1) last=i-1;
      else last=i;
      break;
    }
  }

  if (last==0) last=n;

  for (i=0;i<=last;i++)
    if (*(pp+i)>max) {
      max=*(pp+i);
      maxy=i;
    }
  for (i=0;i<maxy;i++)
    if (rnd(*(pp+i))>0) {
      start=i-1;
      break;
    }
}

void geometric(int n,double p)
{
  int i;
  double max=0,sum=0;

  for (i=1;i<=n;i++) {
    *(pp+i)=p*pow(1-p,i-1);
    sum+=*(pp+i);
    *(pc+i)=sum;
    if (rnd(*(pp+i))==0 && rnd(*(pc+i))==1) {
      if (rnd(cdf[i-1])==1) last=i-1;
      else last=i;
      break;
    }
  }

  if (last==0) last=n;

  for (i=0;i<=last;i++)
    if (*(pp+i)>max) {
      max=*(pp+i);
      maxy=i;
    }
  for (i=0;i<maxy;i++)
    if (rnd(*(pp+i))>0) {
      start=i-1;
      break;
    }
}

void hyper(int n,int a,int s)
{
  int i,end;
  double sum=0,max=0;

  end=min(a,s);
  for (i=0;i<=end;i++) {
    *(pp+i)=com(a,i)*com(n-a,s-i)/com(n,s);
    sum+=*(pp+i);
    *(pc+i)=sum;
    if (rnd(*(pp+i))==0 && rnd(*(pc+i))==1) {
      if (rnd(*(pc+(i-1)))==1) last=i-1;
      else last=i;
      break;
    }
  }

  if (last==0) last=end;

  for (i=0;i<=last;i++)
    if (*(pp+i)>max) {
      max=*(pp+i);
      maxy=i;
    }
  for (i=0;i<maxy;i++)
    if (rnd(*(pp+i))>0) {
      start=i-1;
      break;
    }
}

void negbinomial(int n,double p)
{
  int i;
  double max=0,sum=0;

  for (i=n;;i++) {
    *(pp+i)=com(i-1,n-1)*pow(p,n)*pow(1-p,i-n);
    sum+=*(pp+i);
    *(pc+i)=sum;
    if (rnd(*(pp+i))==0 && rnd(*(pc+i))==1) {
      if (rnd(*(pc+(i-1)))==1) last=i-1;
      else last=i;
      break;
    }
  }

  if (last==0) last=n;

  for (i=0;i<=last;i++)
    if (*(pp+i)>max) {
      max=*(pp+i);
      maxy=i;
    }
  for (i=0;i<maxy;i++)
    if (rnd(*(pp+i))>0) {
      start=i-1;
      break;
    }
}

void init_graph(void)
{
  int xasp,yasp,graphdrive=DETECT,graphmode;
  initgraph(&graphdrive,&graphmode,"");
}

void board(void)
{
  int xn,yn;
  char str[]="Any key pressed continue !!";
  x0=60;
  y0=getmaxy()-50;
  xn=getmaxx()-x0;
  yn=getmaxy()-y0;
  x_range=abs(xn-x0);
  y_range=abs(yn-y0);
  cleardevice();
  rectangle(0,0,getmaxx(),getmaxy());
  setcolor(12);
  outtextxy(getmaxx()/2-(textwidth(str)/2),getmaxy()-textheight(str)-2,str);
  setcolor(15);
  line(x0,y0,x0,yn);
  line(x0,y0,xn,y0);
  gprint(x0-20,y0-3,"%1d",0);
}

int gprint(int x,int y,char *fmt, ...)
{
  va_list arg;
  char str[50];
  int cnt;

  va_start(arg,fmt);
  cnt=vsprintf(str,fmt,arg);
  outtextxy(x,y,str);
  va_end(arg);
  return cnt;
}

void ppoly(char *str)
{
  int i,j,y,temp=0;

  board();

  setcolor(11);
  settextstyle(0,0,2);
  outtextxy(getmaxx()/2-(textwidth(str)/2),textheight(str)-10,str);
  settextstyle(0,0,0);
  setcolor(10);
  gprint(x0-20,y0-y_range-15,"P");
  gprint(x0+x_range+20,y0+20,"X");
  setcolor(15);

  for (i=1;i<=20;i++) {
    line(x0,y0-y_range*i/20,x0-4,y0-y_range*i/20);
    gprint(x0-55,y0-y_range*i/20-3,"%6.4f",*(pp+maxy)*i/20);
    for (j=1;j<=50;j++)
      putpixel(x0+x_range*j/50,y0-y_range*i/20,12);
  }

  if (last<=25) {
    if (go=='3') temp=1;
    for (i=temp;i<=last;i++) {
      circle(x0+x_range/last*i,y0-y_range/(*(pp+maxy))**(pp+i),3);
      setfillstyle(1,9);
      floodfill(x0+x_range/last*i-1,y0-y_range/(*(pp+maxy))**(pp+i)-1,WHITE);
    }

    if (go!='3')
      for (i=1;i<=80;i++)
	putpixel(x0+x_range/last*maxy,y0-y_range*i/80,14);
    for (i=1;i<=last;i++) {
      if (last>=20 && i % 2==1) y=20;
      else y=8;

      line(x0+x_range/last*i,y0,x0+x_range/last*i,y0+4);
      gprint(x0+x_range/last*i-3,y0+y,"%d",i);
      if (go!='3' || i!=1)
	line(x0+x_range/last*(i-1),y0-y_range/(*(pp+maxy))**(pp+(i-1)),
	     x0+x_range/last*i,y0-y_range/(*(pp+maxy))**(pp+i));
    }
  }
  else {
    if (start!=-1)
      for (i=start;i<=last;i++) {
	  circle(x0+x_range/(last-start)*(i-start),
		 y0-y_range/(*(pp+maxy))**(pp+i),3);
	  setfillstyle(1,9);
	  floodfill(x0+x_range/(last-start)*(i-start)-1,
		    y0-y_range/(*(pp+maxy))**(pp+i)-1,WHITE);
      }
    else
      for (i=0;i<=last;i++) {
	circle(x0+x_range/(last-start)*i,
	       y0-y_range/(*(pp+maxy))**(pp+i),3);
	setfillstyle(1,9);
	floodfill(x0+x_range/(last-start)*i-1,
		  y0-y_range/(*(pp+maxy))**(pp+i)-1,WHITE);
     }

    if (last-start>30) {
      if (start+1!=0) {
	line(x0+x_range/(last-start),y0,x0+x_range/(last-start),y0+6);
	line(x0+x_range/(last-start)*(maxy-start),y0,
	     x0+x_range/(last-start)*(maxy-start),y0+6);
	line(x0+x_range/(last-start)*(last-start),y0,
	   x0+x_range/(last-start)*(last-start),y0+6);
	gprint(x0+x_range/(last-start)-3,y0+10,"%d",start+1);
	gprint(x0+x_range/(last-start)*(maxy-start)-3,y0+10,"%d",maxy);
	gprint(x0+x_range/(last-start)*(last-start)-8,y0+10,"%d",last);
      }
      else {
	line(x0+x_range/(last-start)*(maxy-start-1),y0,
	     x0+x_range/(last-start)*(maxy-start-1),y0+6);
	line(x0+x_range/(last-start)*(last-start-1),y0,
	   x0+x_range/(last-start)*(last-start-1),y0+6);
	gprint(x0+x_range/(last-start)*(maxy-start-1)-3,y0+10,"%d",maxy);
	gprint(x0+x_range/(last-start)*(last-start-1)-8,y0+10,"%d",last);
      }
      if (go!='3')
	for (i=1;i<=80;i++)
	  if (start+1!=0)
	    putpixel(x0+x_range/(last-start)*(maxy-start),y0-y_range*i/80,14);
	  else
	    putpixel(x0+x_range/(last-start)*(maxy-start-1),
		     y0-y_range*i/80,14);
    }
    else if (start==-1) {
      for (i=1;i<=last-start-1;i++) {
	if (last-start>=15 && i % 2==1) y=20;
	else y=8;
	line(x0+x_range/(last-start)*i,y0,x0+x_range/(last-start)*i,y0+4);
	gprint(x0+x_range/(last-start)*i-7,y0+y,"%d",start+i+1);
      }
    }
    else
      for (i=1;i<=last-start;i++) {
	if (last-start>=15 && i % 2==1) y=20;
	else y=8;
	line(x0+x_range/(last-start)*i,y0,x0+x_range/(last-start)*i,y0+4);
	gprint(x0+x_range/(last-start)*i-7,y0+y,"%d",start+i);
      }

    if (start!=-1) {
      for (i=start+1;i<=last;i++)
	line(x0+x_range/(last-start)*(i-start-1),
	     y0-y_range/(*(pp+maxy))**(pp+(i-1)),
	     x0+x_range/(last-start)*(i-start),
	     y0-y_range/(*(pp+maxy))**(pp+i));
    }
    else {
       for(i=0;i<=last-start;i++)
	line(x0+x_range/(last-start)*(i-start-1),
	     y0-y_range/(*(pp+maxy))**(pp+i),
	     x0+x_range/(last-start)*(i-start),
	     y0-y_range/(*(pp+maxy))**(pp+(i+1)));
    }

    for (i=1;i<=80;i++)
      if (start!=-1)
	putpixel(x0+x_range/(last-start)*(maxy-start),y0-y_range*i/80,14);
      else
	putpixel(x0+x_range/(last-start)*(maxy-start-1),y0-y_range*i/80,14);
  }
  getch();
}

void cpoly(char *str)
{
  int i,j,y,temp=0;

  board();

  setcolor(11);
  settextstyle(0,0,2);
  outtextxy(getmaxx()/2-(textwidth(str)/2),textheight(str)-10,str);
  settextstyle(0,0,0);
  setcolor(10);
  gprint(x0-20,y0-y_range-15,"P");
  gprint(x0+x_range+20,y0+20,"X");
  setcolor(15);

  for (i=1;i<=20;i++) {
    line(x0,y0-y_range*i/20,x0-4,y0-y_range*i/20);
    gprint(x0-55,y0-y_range*i/20-3,"%6.4f",(double)1/20*i);
    for (j=1;j<=50;j++)
      putpixel(x0+x_range*j/50,y0-y_range*i/20,12);
  }

  if (last<=25) {
    if (go=='3') temp=1;
    for (i=temp;i<=last;i++) {
      circle(x0+x_range/last*i,y0-y_range**(pc+i),3);
      setfillstyle(1,9);
      floodfill(x0+x_range/last*i-1,y0-y_range**(pc+i)-1,WHITE);
    }

    if (go!='3')
      for (i=1;i<=80;i++)
	putpixel(x0+x_range/last*maxy,y0-y_range*i/80,14);

    for (i=1;i<=last;i++) {
      if (last>=20 && i % 2==1) y=20;
      else y=8;

      line(x0+x_range/last*i,y0,x0+x_range/last*i,y0+4);
      gprint(x0+x_range/last*i-3,y0+y,"%d",i);
      if (go!='3' || i!=1)
	line(x0+x_range/last*(i-1),y0-y_range**(pc+(i-1)),
    	     x0+x_range/last*i,y0-y_range**(pc+i));
    }
  }
  else {
    if (start!=-1)
      for (i=start;i<=last;i++) {
	circle(x0+x_range/(last-start)*(i-start),y0-y_range**(pc+i),3);
	setfillstyle(1,9);
	floodfill(x0+x_range/(last-start)*(i-start)-1,
		  y0-y_range**(pc+i)-1,WHITE);
      }

    else
      for (i=0;i<=last;i++) {
	circle(x0+x_range/(last-start)*i,y0-y_range**(pc+i),3);
	setfillstyle(1,9);
	floodfill(x0+x_range/(last-start)*i-1,y0-y_range**(pc+i)-1,WHITE);
      }

    if (last-start>30) {
      if (start+1!=0) {
	line(x0+x_range/(last-start),y0,x0+x_range/(last-start),y0+6);
	line(x0+x_range/(last-start)*(maxy-start),y0,
	     x0+x_range/(last-start)*(maxy-start),y0+6);
	line(x0+x_range/(last-start)*(last-start),y0,
	   x0+x_range/(last-start)*(last-start),y0+6);
	gprint(x0+x_range/(last-start)-3,y0+10,"%d",start+1);
	gprint(x0+x_range/(last-start)*(maxy-start)-3,y0+10,"%d",maxy);
	gprint(x0+x_range/(last-start)*(last-start)-8,y0+10,"%d",last);
      }
      else {
	line(x0+x_range/(last-start)*(maxy-start-1),y0,
	     x0+x_range/(last-start)*(maxy-start-1),y0+6);
	line(x0+x_range/(last-start)*(last-start-1),y0,
	   x0+x_range/(last-start)*(last-start-1),y0+6);
	gprint(x0+x_range/(last-start)*(maxy-start-1)-3,y0+10,"%d",maxy);
	gprint(x0+x_range/(last-start)*(last-start-1)-8,y0+10,"%d",last);
      }

      for (i=1;i<=80;i++)
	if (start+1!=0)
	  putpixel(x0+x_range/(last-start)*(maxy-start),y0-y_range*i/80,14);
	else
	  putpixel(x0+x_range/(last-start)*(maxy-start-1),y0-y_range*i/80,14);
    }
    else
      if (start==-1) {
	for (i=1;i<=last-start-1;i++) {
	  if (last-start>=15 && i % 2==1) y=20;
	  else y=8;
	  line(x0+x_range/(last-start)*i,y0,x0+x_range/(last-start)*i,y0+4);
	  gprint(x0+x_range/(last-start)*i-7,y0+y,"%d",start+i+1);
	}
      }
      else
	for (i=1;i<=last-start;i++) {
	  if (last-start>=15 && i % 2==1) y=20;
	  else y=8;
	  line(x0+x_range/(last-start)*i,y0,x0+x_range/(last-start)*i,y0+4);
	  gprint(x0+x_range/(last-start)*i-7,y0+y,"%d",start+i);
	}

      if (start!=-1) {
	for (i=start+1;i<=last;i++)
	line(x0+x_range/(last-start)*(i-start-1),y0-y_range**(pc+(i-1)),
	     x0+x_range/(last-start)*(i-start),y0-y_range**(pc+i));
      }
      else {
	for(i=0;i<=last-start-2;i++)
	  line(x0+x_range/(last-start)*(i-start-1),
	       y0-y_range**(pc+i),x0+x_range/(last-start)*(i-start),
	       y0-y_range**(pc+(i+1)));
      }

      for (i=1;i<=80;i++)
      if (start!=-1)
	putpixel(x0+x_range/(last-start)*(maxy-start),y0-y_range*i/80,14);
      else
	putpixel(x0+x_range/(last-start)*(maxy-start-1),y0-y_range*i/80,14);

  }
  getch();
}
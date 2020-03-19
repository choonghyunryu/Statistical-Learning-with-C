/* Program : BINOMIAL.C ver 0.1       */
/* Author  : Ryu choong hyun          */
/* Date    : 93.7.22		      */
/* Note    : Binomial Distribtion &   */
/*           Histogram & Poly,etc.    */

#include <stdio.h>
#include <stdarg.h>
#include <conio.h>
#include <math.h>
#include <graphics.h>

#define round(x) ((x>0) ? floor(x+.5) : ceil(x-.5))
#define MAX 1000

double fac(unsigned x);
double com(unsigned n,unsigned r);
double rnd(double x);
void get_bdf(int n,double p);

int gprint(int x,int y,char *fmt, ...);
void init_graph(void);
void board(void);
void bpoly(void);
void cpoly(void);

int n,start,last,pmax,maxy;
double bpdf[MAX],bcdf[MAX];

int x0,y0,x_range,y_range;

void main(void)
{
  int i;
  double p,mean,var,dev;
  clrscr();
  printf("\n\t ** Binomial Distribution **\n");
  printf("\n\t    N = ");
  scanf("%d",&n);
  printf("\n\t    P = ");
  scanf("%lf",&p);
  mean=n*p;
  var=n*p*(1-p);
  dev=sqrt(var);
  clrscr();
  printf("\n\t\t ** Binomial Distribution **\n\n");
  printf("\t\t     N         = %4d\n",n);
  printf("\t\t     P         = %8.6f\n",p);
  printf("\t\t     Mean      = %lf\n",mean);
  printf("\t\t     Variance  = %lf\n",var);
  printf("\t\t     Deviation = %lf\n\n",dev);
  getch();

  get_bdf(n,p);

  clrscr();
  printf("\t\t     X      PDF        CDF\n");
  if (start==-1)
    for (i=start+1;i<=last;i++) {
      printf("\t\t  %4d   %lf    %lf\n",i,bpdf[i],bcdf[i]);
      if (i % 21==20) getch();
    }
  else
    for (i=start;i<=last;i++) {
      printf("\t\t  %4d   %lf    %lf\n",i,bpdf[i],bcdf[i]);
      if (i % 21==20) getch();
    }

  getch();

  init_graph();
  bpoly();
  cpoly();
  closegraph();
}

double fac(unsigned x)
{
  return ((x==0) ? 1 : x*fac(x-1));
}

double com(unsigned n,unsigned r)
{
  return ((n==r || r==0) ? 1 : (n<r) ? 0 : fac(n)/(fac(r)*fac(n-r)));
}

double rnd(double x)
{
  return round(x*1000000)/1000000;
}

void get_bdf(int n,double p)
{
  int i;
  double max=0,sum=0;

  for (i=0;i<=n;i++) {
    bpdf[i]=com(n,i)*pow(p,i)*pow(1-p,n-i);
    sum+=bpdf[i];
    bcdf[i]=sum;
    if (rnd(bpdf[i])==0 && rnd(bcdf[i])==1) {
      if (rnd(bcdf[i-1])==1) last=i-1;
      else last=i;
      break;
    }
  }

  if (last==0) last=n;

  for (i=0;i<=last;i++)
    if (bpdf[i]>max) {
      max=bpdf[i];
      maxy=i;
    }
  for (i=0;i<maxy;i++)
    if (rnd(bpdf[i])>0) {
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

void bpoly(void)
{
  int i,j,y;
  char str[]="Binomial Distribution";

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
    gprint(x0-55,y0-y_range*i/20-3,"%6.4f",bpdf[maxy]*i/20);
    for (j=1;j<=50;j++)
      putpixel(x0+x_range*j/50,y0-y_range*i/20,12);
  }

  if (last<=50) {
    for (i=0;i<=last;i++) {
      circle(x0+x_range/last*i,y0-y_range/bpdf[maxy]*bpdf[i],3);
      setfillstyle(1,9);
      floodfill(x0+x_range/last*i-1,y0-y_range/bpdf[maxy]*bpdf[i]-1,WHITE);
    }

    for (i=1;i<=last;i++) {
      if (last>=20 && i % 2==1) y=20;
      else y=8;

      line(x0+x_range/last*i,y0,x0+x_range/last*i,y0+4);
      gprint(x0+x_range/last*i-3,y0+y,"%d",i);
      line(x0+x_range/last*(i-1),y0-y_range/bpdf[maxy]*bpdf[i-1],
	   x0+x_range/last*i,y0-y_range/bpdf[maxy]*bpdf[i]);
    }
  }
  else {
    for (i=start;i<=last;i++) {
      circle(x0+x_range/(last-start)*(i-start),
	     y0-y_range/bpdf[maxy]*bpdf[i],3);
      setfillstyle(1,9);
      floodfill(x0+x_range/(last-start)*(i-start)-1,
		y0-y_range/bpdf[maxy]*bpdf[i]-1,WHITE);
    }

    line(x0+x_range/(last-start),y0,x0+x_range/(last-start),y0+6);
    line(x0+x_range/(last-start)*(maxy-start),y0,
	 x0+x_range/(last-start)*(maxy-start),y0+6);

    for (i=1;i<=80;i++)
      putpixel(x0+x_range/(last-start)*(maxy-start),y0-y_range*i/80,14);

    line(x0+x_range/(last-start)*(last-start),y0,
	 x0+x_range/(last-start)*(last-start),y0+6);
    gprint(x0+x_range/(last-start)-10,y0+10,"%d",start);
    gprint(x0+x_range/(last-start)*(maxy-start)-10,y0+10,"%d",maxy);
    gprint(x0+x_range/(last-start)*(last-start)-10,y0+10,"%d",last);

    for (i=start+1;i<=last;i++)
      line(x0+x_range/(last-start)*(i-start-1),
	   y0-y_range/bpdf[maxy]*bpdf[i-1],x0+x_range/(last-start)*(i-start),
	   y0-y_range/bpdf[maxy]*bpdf[i]);
  }
  getch();
}

void cpoly(void)
{
  int i,j,y;
  char str[]="Cumul-Binomial Distribution";

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

  if (last<=50) {
    for (i=0;i<=last;i++) {
      circle(x0+x_range/last*i,y0-y_range*bcdf[i],3);
      setfillstyle(1,9);
      floodfill(x0+x_range/last*i-1,y0-y_range*bcdf[i]-1,WHITE);
    }

    for (i=1;i<=last;i++) {
      if (last>=20 && i % 2==1) y=20;
      else y=8;

      line(x0+x_range/last*i,y0,x0+x_range/last*i,y0+4);
      gprint(x0+x_range/last*i-3,y0+y,"%d",i);
      line(x0+x_range/last*(i-1),y0-y_range*bcdf[i-1],
	   x0+x_range/last*i,y0-y_range*bcdf[i]);
    }
  }
  else {
    for (i=start;i<=last;i++) {
      circle(x0+x_range/(last-start)*(i-start),y0-y_range*bcdf[i],3);
      setfillstyle(1,9);
      floodfill(x0+x_range/(last-start)*(i-start)-1,
		y0-y_range*bcdf[i]-1,WHITE);
    }

    line(x0+x_range/(last-start),y0,x0+x_range/(last-start),y0+6);
    line(x0+x_range/(last-start)*(maxy-start),y0,
	 x0+x_range/(last-start)*(maxy-start),y0+6);

    for (i=1;i<=80;i++)
      putpixel(x0+x_range/(last-start)*(maxy-start),y0-y_range*i/80,14);

    line(x0+x_range/(last-start)*(last-start),y0,
	 x0+x_range/(last-start)*(last-start),y0+6);
    gprint(x0+x_range/(last-start)-10,y0+10,"%d",start);
    gprint(x0+x_range/(last-start)*(maxy-start)-10,y0+10,"%d",maxy);
    gprint(x0+x_range/(last-start)*(last-start)-10,y0+10,"%d",last);

    for (i=start+1;i<=last;i++)
      line(x0+x_range/(last-start)*(i-start-1),y0-y_range*bcdf[i-1],
	   x0+x_range/(last-start)*(i-start),y0-y_range*bcdf[i]);
  }
  getch();
}
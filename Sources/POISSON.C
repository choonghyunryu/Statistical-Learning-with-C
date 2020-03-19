/* Program : POISSON.C ver 0.1        */
/* Author  : Ryu choong hyun          */
/* Date    : 93.7.22		      */
/* Note    : Poisson Distribtion &    */
/*           Histogram & Poly,etc.    */

#include <stdio.h>
#include <stdarg.h>
#include <conio.h>
#include <math.h>
#include <graphics.h>

#define round(x) ((x>0) ? floor(x+.5) : ceil(x-.5))
#define MAX 1000

double fac(int x);
double rnd(double x);
void poisson(unsigned long n,double mean);

int gprint(int x,int y,char *fmt, ...);
void init_graph(void);
void end_graph(void);
void board(void);
void box(int x1,int y1,int x2,int y2);
void ppoly(void);
void cpoly(void);

unsigned long n;
int start,last,pmax,maxy;
double ppdf[MAX],pcdf[MAX];

int x0,y0,x_range,y_range;

void main(void)
{
  int i;
  char c;
  double p,mean;

  clrscr();
  printf("\n\t ** Poisson Distribution **\n");
  printf("\n\t What is Data ?  (1.N,P  2.MEAN)\n");

  if (getch()=='1') {
    printf("\n\t    N = ");
    scanf("%ld",&n);
    printf("\n\t    P = ");
    scanf("%lf",&p);
    mean=n*p;
    clrscr();
    printf("\n\t\t ** Poisson Distribution **\n\n");
    printf("\t\t     N         = %ld\n",n);
    printf("\t\t     P         = %8.6f\n",p);
    printf("\t\t     Mean      = %lf\n",mean);
    printf("\t\t     Variance  = %lf\n",mean);
    printf("\t\t     Deviation = %lf\n\n",sqrt(mean));
    getch();
    poisson(n,mean);
  }
  else {
    printf("\n\t  Mean =  ");
    scanf("%lf",&mean);
    clrscr();
    printf("\n\t\t ** Poisson Distribution **\n\n");
    printf("\t\t     Mean      = %lf\n",mean);
    printf("\t\t     Variance  = %lf\n",mean);
    printf("\t\t     Deviation = %lf\n\n",sqrt(mean));
    getch();
    poisson(1000,mean);
  }

  clrscr();
  printf("\t\t     X      PDF        CDF\n");
  for (i=0;i<=last;i++) {
    printf("\t\t  %4d   %lf    %lf\n",i,ppdf[i],pcdf[i]);
    if (i % 21==20) getch();
  }
  getch();

  init_graph();
  ppoly();
  cpoly();
  end_graph();
}

double fac(int x)
{
  int i;
  double fac=1;

  if (x==0) return 1;

  for (i=1;i<=x;i++)
    fac*=i;
  return fac;
}

double rnd(double x)
{
  return round(x*1000000)/1000000;
}

void poisson(unsigned long n,double mean)
{
  int i;
  double max=0,sum=0;

  for (i=0;i<=n;i++) {
    ppdf[i]=exp(-mean)*pow(mean,i)/fac(i);
    sum+=ppdf[i];
    pcdf[i]=sum;
    if (rnd(ppdf[i])==0 && rnd(pcdf[i])==1) {
      if (rnd(pcdf[i-1])==1) last=i-1;
      else last=i;
      break;
    }
  }

  if (last==0) last=n;

  for (i=0;i<=last;i++)
    if (ppdf[i]>max) {
      max=ppdf[i];
      maxy=i;
    }
  for (i=0;i<maxy;i++)
    if (rnd(ppdf[i])>0) {
      start=i-1;
      break;
    }
}

void init_graph(void)
{
  int xasp,yasp,graphdrive=DETECT,graphmode;
  initgraph(&graphdrive,&graphmode,"");
}

void end_graph(void)
{
  closegraph();
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
  box(0,0,getmaxx(),getmaxy());
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

void box(int x1,int y1,int x2,int y2)
{
  line(x1,y1,x1,y2);
  line(x1,y1,x2,y1);
  line(x1,y2,x2,y2);
  line(x2,y1,x2,y2);
}

void ppoly(void)
{
  int i,j,y;
  char str[]="Poisson Distribution";

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
    gprint(x0-55,y0-y_range*i/20-3,"%6.4f",ppdf[maxy]*i/20);
    for (j=1;j<=50;j++)
      putpixel(x0+x_range*j/50,y0-y_range*i/20,12);
  }

  if (last<=50) {
    for (i=0;i<=last;i++) {
      circle(x0+x_range/last*i,y0-y_range/ppdf[maxy]*ppdf[i],3);
      setfillstyle(1,9);
      floodfill(x0+x_range/last*i-1,y0-y_range/ppdf[maxy]*ppdf[i]-1,WHITE);
    }

    for (i=1;i<=last;i++) {
      if (last>=20 && i % 2==1) y=20;
      else y=8;

      line(x0+x_range/last*i,y0,x0+x_range/last*i,y0+4);
      gprint(x0+x_range/last*i-3,y0+y,"%d",i);
      line(x0+x_range/last*(i-1),y0-y_range/ppdf[maxy]*ppdf[i-1],
	   x0+x_range/last*i,y0-y_range/ppdf[maxy]*ppdf[i]);
    }
  }
  else {
    for (i=start;i<=last;i++) {
      circle(x0+x_range/(last-start)*(i-start),
	     y0-y_range/ppdf[maxy]*ppdf[i],3);
      setfillstyle(1,9);
      floodfill(x0+x_range/(last-start)*(i-start)-1,
		y0-y_range/ppdf[maxy]*ppdf[i]-1,WHITE);
    }

    line(x0+x_range/(last-start),y0,x0+x_range/(last-start),y0+6);
    line(x0+x_range/(last-start)*(maxy-start),y0,
	 x0+x_range/(last-start)*(maxy-start),y0+6);
    line(x0+x_range/(last-start)*(last-start),y0,
	 x0+x_range/(last-start)*(last-start),y0+6);
    gprint(x0+x_range/(last-start)-10,y0+10,"%d",start);
    gprint(x0+x_range/(last-start)*(maxy-start)-10,y0+10,"%d",maxy);
    gprint(x0+x_range/(last-start)*(last-start)-10,y0+10,"%d",last);

    for (i=start+1;i<=last;i++)
      line(x0+x_range/(last-start)*(i-start-1),
	   y0-y_range/ppdf[maxy]*ppdf[i-1],x0+x_range/(last-start)*(i-start),
	   y0-y_range/ppdf[maxy]*ppdf[i]);
  }
  getch();
}

void cpoly(void)
{
  int i,j,y;
  char str[]="Cumul-Poisson Distribution";

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
      circle(x0+x_range/last*i,y0-y_range*pcdf[i],3);
      setfillstyle(1,9);
      floodfill(x0+x_range/last*i-1,y0-y_range*pcdf[i]-1,WHITE);
    }

    for (i=1;i<=last;i++) {
      if (last>=20 && i % 2==1) y=20;
      else y=8;

      line(x0+x_range/last*i,y0,x0+x_range/last*i,y0+4);
      gprint(x0+x_range/last*i-3,y0+y,"%d",i);
      line(x0+x_range/last*(i-1),y0-y_range*pcdf[i-1],
	   x0+x_range/last*i,y0-y_range*pcdf[i]);
    }
  }
  else {
    for (i=start;i<=last;i++) {
      circle(x0+x_range/(last-start)*(i-start),y0-y_range*pcdf[i],3);
      setfillstyle(1,9);
      floodfill(x0+x_range/(last-start)*(i-start)-1,
		y0-y_range*pcdf[i]-1,WHITE);
    }

    line(x0+x_range/(last-start),y0,x0+x_range/(last-start),y0+6);
    line(x0+x_range/(last-start)*(maxy-start),y0,
	 x0+x_range/(last-start)*(maxy-start),y0+6);
    line(x0+x_range/(last-start)*(last-start),y0,
	 x0+x_range/(last-start)*(last-start),y0+6);
    gprint(x0+x_range/(last-start)-10,y0+10,"%d",start);
    gprint(x0+x_range/(last-start)*(maxy-start)-10,y0+10,"%d",maxy);
    gprint(x0+x_range/(last-start)*(last-start)-10,y0+10,"%d",last);

    for (i=start+1;i<=last;i++)
      line(x0+x_range/(last-start)*(i-start-1),y0-y_range*pcdf[i-1],
	   x0+x_range/(last-start)*(i-start),y0-y_range*pcdf[i]);
  }
  getch();
}
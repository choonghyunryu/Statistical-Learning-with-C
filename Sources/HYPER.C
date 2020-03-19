/* Program : HYPER.C ver 0.2               */
/* Author  : Ryu choong hyun               */
/* Date    : 93.7.23		           */
/* Note    : Hypergeometric Distribtion &  */
/*           Histogram & Poly.             */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <conio.h>
#include <math.h>
#include <graphics.h>

#define round(x) ((x>0) ? floor(x+.5) : ceil(x-.5))
#define MAX 1000

double com(int n,int r);
double rnd(double x);
void get_hyper(int n,int a,int s);

int gprint(int x,int y,char *fmt, ...);
void init_graph(void);
void end_graph(void);
void board(void);
void box(int x1,int y1,int x2,int y2);
void hpoly(void);
void cpoly(void);

int start,last,end,maxy;
double hpdf[MAX],hcdf[MAX];

int x0,y0,x_range,y_range;

void main(void)
{
  int i,n,a,s;
  double mean,var,dev;
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
  printf("\n\t\t **Hypergeometric Distribution **\n\n");
  printf("\t\t         N         = %4d\n",n);
  printf("\t\t         A-N       = %4d\n",a);
  printf("\t\t         Sample    = %4d\n",s);
  printf("\t\t         Mean      = %lf\n",mean);
  printf("\t\t         Variance  = %lf\n",var);
  printf("\t\t         Deviation = %lf\n\n",dev);
  getch();

  get_hyper(n,a,s);

  clrscr();
  printf("\t\t     X      PDF        CDF\n");
  for (i=start+1;i<=last;i++) {
    printf("\t\t  %4d   %lf    %lf\n",i,hpdf[i],hcdf[i]);
    if ((i-start) % 20==0) getch();
  }
  getch();

  init_graph();
  hpoly();
  cpoly();
  end_graph();
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

void get_hyper(int n,int a,int s)
{
  int i;
  double sum=0,max=0;

  end=min(a,s);
  for (i=0;i<=end;i++) {
    hpdf[i]=com(a,i)*com(n-a,s-i)/com(n,s);
    sum+=hpdf[i];
    hcdf[i]=sum;
    if (rnd(hpdf[i])==0 && rnd(hcdf[i])==1) {
      if (rnd(hcdf[i-1])==1) last=i-1;
      else last=i;
      break;
    }
  }

  if (last==0) last=end;

  for (i=0;i<=last;i++)
    if (hpdf[i]>max) {
      max=hpdf[i];
      maxy=i;
    }
  for (i=0;i<maxy;i++)
    if (rnd(hpdf[i])>0) {
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

void hpoly(void)
{
  int i,j,y;
  char str[]="Hypergeometric Distribution";

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
    gprint(x0-55,y0-y_range*i/20-3,"%6.4f",hpdf[maxy]*i/20);
    for (j=1;j<=50;j++)
      putpixel(x0+x_range*j/50,y0-y_range*i/20,12);
  }

  if (last<=25) {
    for (i=0;i<=last;i++) {
      circle(x0+x_range/last*i,y0-y_range/hpdf[maxy]*hpdf[i],3);
      setfillstyle(1,9);
      floodfill(x0+x_range/last*i-1,y0-y_range/hpdf[maxy]*hpdf[i]-1,WHITE);
    }

    for (i=1;i<=last;i++) {
      if (last>=20 && i % 2==1) y=20;
      else y=8;

      line(x0+x_range/last*i,y0,x0+x_range/last*i,y0+4);
      gprint(x0+x_range/last*i-3,y0+y,"%d",i);
      line(x0+x_range/last*(i-1),y0-y_range/hpdf[maxy]*hpdf[i-1],
	   x0+x_range/last*i,y0-y_range/hpdf[maxy]*hpdf[i]);
    }
  }
  else {
    for (i=start;i<=last;i++) {
      circle(x0+x_range/(last-start)*(i-start),
	     y0-y_range/hpdf[maxy]*hpdf[i],3);
      setfillstyle(1,9);
      floodfill(x0+x_range/(last-start)*(i-start)-1,
		y0-y_range/hpdf[maxy]*hpdf[i]-1,WHITE);
    }

    if (last-start>30) {
      line(x0+x_range/(last-start),y0,x0+x_range/(last-start),y0+6);
      line(x0+x_range/(last-start)*(maxy-start),y0,
	   x0+x_range/(last-start)*(maxy-start),y0+6);

      for (i=1;i<=80;i++)
	putpixel(x0+x_range/(last-start)*(maxy-start),y0-y_range*i/80,14);

      line(x0+x_range/(last-start)*(last-start),y0,
	   x0+x_range/(last-start)*(last-start),y0+6);
      gprint(x0+x_range/(last-start)-10,y0+10,"%d",start+1);
      gprint(x0+x_range/(last-start)*(maxy-start)-10,y0+10,"%d",maxy);
      gprint(x0+x_range/(last-start)*(last-start)-10,y0+10,"%d",last);
    }
    else
      for (i=1;i<=last-start;i++) {
	if (last-start>=15 && i % 2==1) y=20;
	else y=8;
	line(x0+x_range/(last-start)*i,y0,x0+x_range/(last-start)*i,y0+4);
	gprint(x0+x_range/(last-start)*i-10,y0+y,"%d",last+i);
      }

    for (i=start+1;i<=last;i++)
      line(x0+x_range/(last-start)*(i-start-1),
	   y0-y_range/hpdf[maxy]*hpdf[i-1],x0+x_range/(last-start)*(i-start),
	   y0-y_range/hpdf[maxy]*hpdf[i]);
  }
  getch();
}

void cpoly(void)
{
  int i,j,y;
  char str[]="Cumul-Hypergeometric Distribution";

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
    for (i=0;i<=last;i++) {
      circle(x0+x_range/last*i,y0-y_range*hcdf[i],3);
      setfillstyle(1,9);
      floodfill(x0+x_range/last*i-1,y0-y_range*hcdf[i]-1,WHITE);
    }

    for (i=1;i<=last;i++) {
      if (last>=20 && i % 2==1) y=20;
      else y=8;

      line(x0+x_range/last*i,y0,x0+x_range/last*i,y0+4);
      gprint(x0+x_range/last*i-3,y0+y,"%d",i);
      line(x0+x_range/last*(i-1),y0-y_range*hcdf[i-1],
	   x0+x_range/last*i,y0-y_range*hcdf[i]);
    }
  }
  else {
    for (i=start;i<=last;i++) {
      circle(x0+x_range/(last-start)*(i-start),y0-y_range*hcdf[i],3);
      setfillstyle(1,9);
      floodfill(x0+x_range/(last-start)*(i-start)-1,
		y0-y_range*hcdf[i]-1,WHITE);
    }

    if (last-start>30) {
      line(x0+x_range/(last-start),y0,x0+x_range/(last-start),y0+6);
      line(x0+x_range/(last-start)*(maxy-start),y0,
	 x0+x_range/(last-start)*(maxy-start),y0+6);

      for (i=1;i<=80;i++)
	putpixel(x0+x_range/(last-start)*(maxy-start),y0-y_range*i/80,14);

      line(x0+x_range/(last-start)*(last-start),y0,
	   x0+x_range/(last-start)*(last-start),y0+6);
      gprint(x0+x_range/(last-start)-10,y0+10,"%d",start+1);
      gprint(x0+x_range/(last-start)*(maxy-start)-10,y0+10,"%d",maxy);
      gprint(x0+x_range/(last-start)*(last-start)-10,y0+10,"%d",last);
    }
    else
      for (i=1;i<=last-start;i++) {
	if (last-start>=15 && i % 2==1) y=20;
	else y=8;
	line(x0+x_range/(last-start)*i,y0,x0+x_range/(last-start)*i,y0+4);
	gprint(x0+x_range/(last-start)*i-10,y0+y,"%d",last+i);
      }

    for (i=start+1;i<=last;i++)
    line(x0+x_range/(last-start)*(i-start-1),y0-y_range*hcdf[i-1],
	 x0+x_range/(last-start)*(i-start),y0-y_range*hcdf[i]);
  }
  getch();
}
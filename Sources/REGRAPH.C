/* Program : REGRAPH.C ver 0.1                */
/* Author  : Ryu choong hyun                  */
/* Date    : 94.7.27.		      	      */
/* Note    : Regression Analysis Graph &      */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <stdarg.h> /* for gprint() */
#include <graphics.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define PI M_PI
#define FACTOR 50

void get_minmax(int col,double *min,double *max);
void get_mean_ss(int col,double *mean,double *ss);
double get_sxy(void);
double get_maxe(void);

int gprint(int x,int y,char *fmt, ...);
void init_graph(void);
void board(void);
void regress(void);
void residual(void);

int row,x0,y0,x_range,y_range;
double data[FACTOR][2],min[2],max[2],mean[2],ss[2],
       b0,b1,x_temp,y_temp;

void main(int argc,char *argv[])
{
  FILE *stream;
  int i,j;
  double num,s_xy;

  clrscr();

  if (argc<=1) {
    puts("Usage : regraph datafile");
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

  for (i=0;i<2;i++) {
    get_minmax(i,&min[i],&max[i]);
    get_mean_ss(i,&mean[i],&ss[i]);
  }

  s_xy=get_sxy();
  b1=s_xy/ss[0];
  b0=mean[1]-b1*mean[0];

  init_graph();
  regress();
  closegraph();
}

void get_minmax(int col,double *min,double *max)
{
  int i;

  *min=data[0][col];
  *max=data[0][col];

  for (i=0;i<row;i++) {
    if (*min>data[i][col]) *min=data[i][col];
    if (*max<data[i][col]) *max=data[i][col];
  }
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

void init_graph(void)
{
  int graphdrive=DETECT,graphmode;
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

  x_temp=(x_range-10)/(max[0]-min[0]);
  y_temp=(y_range-10)/(max[1]-min[1]);

  rectangle(0,0,getmaxx(),getmaxy());
  setcolor(LIGHTRED);
  outtextxy(getmaxx()/2-(textwidth(str)/2),getmaxy()-textheight(str)-2,str);
  setcolor(WHITE);
  line(x0,y0,x0,yn);
  line(x0,y0,xn,y0);
  gprint(x0-20,y0+3,"%1d",0);
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

void regress(void)
{
  int i,j;
  char str[]="REGRESSION ANALYSIS";

  board();

  setcolor(LIGHTCYAN);
  settextstyle(DEFAULT_FONT,HORIZ_DIR,2);
  outtextxy(getmaxx()/2-(textwidth(str)/2),textheight(str)-10,str);
  settextstyle(DEFAULT_FONT,HORIZ_DIR,0);

  setcolor(LIGHTBLUE);
  outtextxy(x0-35,y0-y_range-15,"Y");
  outtextxy(x0+x_range-30,y0+30,"X");

  setcolor(YELLOW);

  line(x0,y0-10,x0-4,y0-10);
  gprint(x0-50,y0-15,"%5.1f",min[1]);
  line(x0+10,y0,x0+10,y0+4);
  gprint(x0-5,y0+10,"%5.1f",min[0]);

  for (j=1;j<=3;j++) {
    line(x0,y0-(y_range-10)/3*j,x0-4,y0-(y_range-10)/3*j);
    gprint(x0-50,y0-(y_range-10)/3*j-2,"%5.1f",min[1]+(max[1]-min[1])/3*j);
  }

  for (i=1;i<=5;i++) {
    line(x0+(x_range-10)/5*i,y0,x0+(x_range-10)/5*i,y0+4);
    gprint(x0+(x_range-10)/5*i-23,y0+10,"%5.1f",min[0]+(max[0]-min[0])/5*i);
  }

  for (i=0;i<row;i++) {
    circle((x0+10)+x_temp*(data[i][0]-min[0]),
           (y0-10)-y_temp*(data[i][1]-min[1]),3);
    setfillstyle(SOLID_FILL,LIGHTBLUE);
    floodfill((x0+10)+x_temp*(data[i][0]-min[0]),
              (y0-10)-y_temp*(data[i][1]-min[1]),YELLOW);
  }

  setcolor(WHITE);
  if (b1>0) line(x0+10+x_temp*((min[1]-b0)/b1-min[0]),y0-10,x0+x_range-10,
                 y0-10-y_temp*((b0+b1*max[0])-min[1]));
  else if (b1!=0) line(x0+10,y0-10-y_temp*((b0+b1*min[0])-min[1]),
                  x0+x_range,y0-10-y_temp*((b0+b1*max[0])-min[1]));
  else line(x0,y0-10-y_temp*(b0-min[1]),x0+x_range,y0-10-y_temp*(b0-min[1]));

  getch();
}


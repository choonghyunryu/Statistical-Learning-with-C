/* Program : TSRIES.C ver 0.1   */
/* Author  : Ryu choong hyun    */
/* Date    : 94.7.27.		*/
/* Note    : Time Series Graph  */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <stdarg.h> /* for gprint() */
#include <graphics.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define FACTOR 50

void get_minmax(double *min,double *max);
void time_graph(int interval);

int gprint(int x,int y,char *fmt, ...);
void init_graph(void);
void board(void);

int row,x0,y0,x_range,y_range;
double data[FACTOR],time[FACTOR],min,max,x_temp,y_temp;

void main(int argc,char *argv[])
{
  FILE *stream;
  int i,j,interval;
  double num;

  clrscr();

  if (argc<=1) {
    puts("Usage : tsries datafile");
    exit(EXIT_FAILURE);
  }

  stream=fopen(argv[1],"rt");
  if (stream==NULL) printerrmsg("File not found !!");

  fscanf(stream,"%d\n",&row);

  for (i=0;i<row;i++)
    for (j=0;j<2;j++) {
      fscanf(stream,"%lf\n",&num);
      if (j==0) time[i]=num;
      else data[i]=num;
    }

  get_minmax(&min,&max);

  printf("\t\t *** Time Series Graph ***\n\n");
  printf("\t\t       Line Interval = ");
  scanf("%d",&interval);

  init_graph();
  time_graph(interval);
  closegraph();
}

void get_minmax(double *min,double *max)
{
  int i;

  *min=data[0];
  *max=data[0];

  for (i=0;i<row;i++) {
    if (*min>data[i]) *min=data[i];
    if (*max<data[i]) *max=data[i];
  }
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

  x_temp=(x_range-10)/(double)row;
  y_temp=(y_range-10)/(max-min);

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

void time_graph(int interval)
{
  int i,j;
  char str[]="TIME SERIES GRAPH";

  board();

  setcolor(LIGHTCYAN);
  settextstyle(DEFAULT_FONT,HORIZ_DIR,2);
  outtextxy(getmaxx()/2-(textwidth(str)/2),textheight(str)-10,str);
  settextstyle(DEFAULT_FONT,HORIZ_DIR,0);

  setcolor(LIGHTBLUE);
  outtextxy(x0-35,y0-y_range-15,"DATA");
  outtextxy(x0+x_range-30,y0+30,"TIME");

  setcolor(RED);
  for (i=interval;i<=row;i+=interval)
    line(x0+10+x_temp*(i-1),y0,x0+10+x_temp*(i-1),y0-y_range);

  setcolor(WHITE);
  line(x0,y0-10,x0-4,y0-10);
  gprint(x0-60,y0-15,"%5.1f",min);

  for (j=1;j<=3;j++) {
    line(x0,y0-(y_range-10)/3*j,x0-4,y0-(y_range-10)/3*j);
    gprint(x0-60,y0-(y_range-10)/3*j-2,"%5.1f",min+(max-min)/3*j);
  }

  j=0;
  for (i=interval;i<=row;i+=interval) {
    line(x0+10+x_temp*(i-1),y0,x0+10+x_temp*(i-1),y0+4);
    if (row/interval>10 && j%2==1) {
      gprint(x0+10+x_temp*(i-1)-10,y0+20,"%d",(int)time[i-1]);
    }
    else gprint(x0+10+x_temp*(i-1)-10,y0+10,"%d",(int)time[i-1]);
    j++;
  }

  for (i=0;i<row;i++) {
    if (i<row-1) line((x0+10)+x_temp*i,(y0-10)-y_temp*(data[i]-min),
		      (x0+10)+x_temp*(i+1),(y0-10)-y_temp*(data[i+1]-min));
    circle((x0+10)+x_temp*i,(y0-10)-y_temp*(data[i]-min),3);
    setfillstyle(SOLID_FILL,LIGHTBLUE);
    floodfill((x0+10)+x_temp*i-1,(y0-10)-y_temp*(data[i]-min)-1,WHITE);

  }
  getch();
}


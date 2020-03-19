/* Program : TANAL.C ver 0.1     */
/* Note    : Time Siries Analysis */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h> /* for gprint() */
#include <graphics.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define FACTOR 50

void get_bi(double *b0,double *b1);                          /* Â­A¬å· ‰® 
*/
void get_moving_mean(int n);                                 /* ·¡•·Íw‹E */
void get_s(int n);                                           /* ‰¸é»¡® */
void get_ci(int n);
void get_minmax(int opt,double x[],double *min,double *max);
void time_graph(char *c,double data[]);

int gprint(int x,int y,char *fmt, ...);
void init_graph(void);
void board(void);

int row,n,x0,y0,x_range,y_range;
double data[FACTOR],time[FACTOR],move[FACTOR],tc[FACTOR],si[FACTOR],s[12],
       t[FACTOR],c[FACTOR],ii[FACTOR],min,max,x_temp,y_temp;

void main(int argc,char *argv[])
{
  FILE *stream;
  int i,j;
  double num,b0,b1,mafe,mape,rmse,ytemp,des_data[FACTOR],ec,et,es,ey;

  clrscr();

  if (argc<=1) {
    puts("Usage : tanal datafile");
    exit(EXIT_FAILURE);
  }

  stream=fopen(argv[1],"rt");
  if (stream==NULL) printerrmsg("File not found !!");

  fscanf(stream,"%d\n",&row);                                 /* ‰ÅÀiÃ¡· ˆ•® 
*/

  for (i=0;i<row;i++)                                          /* ¸aža· ·³b 
*/
    for (j=0;j<2;j++) {
      fscanf(stream,"%lf\n",&num);
      if (j==0) time[i]=num;
      else data[i]=num;
    }

  printf("\t\t Moving Average Number = ");                    /* ·¡•·Íw‹E ¹·ŸA 
¬å¸÷ */
  scanf("%d",&n);

  get_bi(&b0,&b1);                                             /* Â­A¬å  ‰® 
*/
  for (i=0;i<row;i++) t[i]=b0+b1*(i+1);

  get_moving_mean(n);                                         /* ·¡•·Íw‹E ‰¬e 
*/
  get_s(n);                                                  /* ‰¸÷»¡® •¡Â‰ 
*/
  get_ci(n);                                                 /* C*I ‰¬e */

  clrscr();

  printf("\t ** Time Siries Analysis **\n\n");
  printf("\t  << Simple Trend Line >>\n");                   /* Â­A¬å */ 
  printf("\tT=%lf+%lfX  (X=1 to %d)\n\n",b0,b1,row);

  ytemp=0;
  for (i=0;i<row;i++) ytemp+=fabs(data[i]-t[i]);
  mafe=ytemp/row;

  ytemp=0;
  for (i=0;i<row;i++) ytemp+=fabs(data[i]-t[i])/data[i];
  mape=ytemp/row*100;

  ytemp=0;
  for (i=0;i<row;i++) ytemp+=pow(data[i]-t[i],2);
  rmse=sqrt(ytemp/row);

  printf("\n\t << Trend Line's Measures >>\n");
  printf("\t  MAFE     MAPE     RMSE\n");
  printf("\t %4.3f    %4.1f    %4.3f\n",mafe,mape,rmse);

  printf("\n\n\t << Moving Average  & Seasonal Index >>\n"); /* ·¡•·Íw‹E‰Á  ‰
¸é»¡® */
  printf("\t  Time       Y       M-M      T*C      S*I");

  for (i=0;i<row;i++) {
    if (i<n/2 || i>row-n/2) {
      printf("\n");
      printf("\t%6.f   %8.3f\n",time[i],data[i]);
    }
    else if (i!=row-n/2) {
      printf("\t\t\t   %7.1f\n",move[i]);
      printf("\t%6.f   %8.3f\t    %7.1f  %7.1f\n",
             time[i],data[i],tc[i],si[i]);
    }
    else {
      printf("\t\t\t   %7.1f\n",move[i]);
      printf("\t%6.f   %8.3f\n",time[i],data[i]);
    }
  }

  printf("\n\n\t << Adjusted Seasonal Index >>\n");            /* ¹¡¸÷–E  ‰¸é
»¡® */  
  for (i=0;i<n;i++) printf("\t\t  S%d = %5.1f\n",i+1,s[i]);

  printf("\n\n\t\t       << Time Siries's Element >>\n");
  printf("\t  Time       Y       D-Y      T       S       C       I\n");

  for (i=0;i<row;i++) {
    des_data[i]=data[i]/s[i%n]*100;
    if (i==0 || i==row-1)
      printf("\t%6.f  %8.3f  %8.3f  %5.1f   %5.1f     -       -\n",
             time[i],data[i],des_data[i],t[i],s[i%n]);
    else printf("\t%6.f  %8.3f  %8.3f  %5.1f   %5.1f   %5.1f   %5.1f\n",
                time[i],data[i],des_data[i],t[i],s[i%n],c[i],ii[i]);
  }

  getch();

  init_graph();
  get_minmax(1,t,&min,&max);         
  time_graph("T",t);                 /* Â­A¥e•· ‹aœÏa */
  get_minmax(0,s,&min,&max);
  time_graph("S",s);                 /* ‰¸é¥e•· ‹aœÏa */ 
  get_minmax(1,c,&min,&max);
  time_graph("C",c);                 /* º‹¡¥e•· ‹aœÏa */
  get_minmax(1,ii,&min,&max); 
  time_graph("I",ii);                /* ¦‰‹AÃ¢¥e•· ‹aœÏa */
  closegraph();

  printf("\n\n\t\t << Estimate Time Siries >>\n");
  printf("\t\t   Estimate C      =  ");                /* µÃb º‹¡¥e•· ˆt */   
    
  scanf("%lf",&ec);

  printf("\n\t    Time    T      C      S       Y\n");
  for (i=row;i<row+n;i++) {
    et=b0+b1*(i+1);
    es=s[i%n];
    ey=et*ec*es/10000;
    if (n==4) printf("\t  %6.f  %5.1f  %5.1f  %5.1f %8.3f\n",
                     time[row-n+(i-row)]+10,et,ec,es,ey);
    else printf("\t  %6.f  %5.1f  %5.1f  %5.1f %8.3f\n",
                time[row-n+(i-row)]+100,et,ec,es,ey);
  }
  getch();
}

void get_bi(double *b0,double *b1)
{
  int i;
  double sum_sxy=0,sum_sxx=0,sum_x=0,sum_y=0,xmean,ymean;

  for (i=0;i<row;i++) {
    sum_x+=i+1;
    sum_y+=data[i];
  }

  xmean=sum_x/row;
  ymean=sum_y/row;

  for (i=0;i<row;i++) {
    sum_sxy+=(i+1-xmean)*(data[i]-ymean);
    sum_sxx+=pow(i+1-xmean,2);
  }

  *b1=sum_sxy/sum_sxx;
  *b0=ymean-*b1*xmean;
}

void get_moving_mean(int n)
{
  int i,j;
  double sum;

  for (i=n/2;i<=row-n/2;i++) {
    sum=0;
    for (j=i-n/2;j<i+n/2;j++) sum+=data[j];
    move[i]=sum/n;
  }

  for (i=n/2;i<row-n/2;i++) {
    tc[i]=(move[i]+move[i+1])/2;
    si[i]=data[i]/tc[i]*100;
  }
}

void get_s(int n)
{
  int i;
  double sum=0,temp[12]={0,};

  for (i=0;i<row;i++) temp[i%n]+=si[i];

  for (i=0;i<n;i++) {
    s[i]=temp[i]/(row/n-1);
    sum+=s[i];
  }

  for (i=0;i<n;i++) s[i]*=n*100/sum;

}

void get_ci(int n)
{

  int i,j;
  double sum,ci[FACTOR];

  for (i=0;i<row;i++) ci[i]=data[i]/(t[i]*s[i%n])*10000;

  for (i=1;i<row-1;i++) {
    sum=0;
    for (j=i-1;j<=i+1;j++) sum+=ci[j];
    c[i]=sum/3;
    ii[i]=ci[i]/c[i]*100;
  }
}

void get_minmax(int opt,double x[],double *min,double *max)
{
  int i;

  if (x[0]!=0) *min=x[0];
  else *min=x[1];
  *max=x[0];

  if (opt==1) {
    if (x[0]==0)
      for (i=1;i<row-1;i++) {
        if (*min>x[i]) *min=x[i];
        if (*max<x[i]) *max=x[i];
      }
    else
      for (i=0;i<row;i++) {
        if (*min>x[i]) *min=x[i];
        if (*max<x[i]) *max=x[i];
      }
  }
  else {
    for (i=0;i<n;i++) {
      if (*min>x[i]) *min=x[i];
      if (*max<x[i]) *max=x[i];
    }
  }
}

void time_graph(char *c,double data[])
{
  int i,j;
  char str[]="TIME SERIES GRAPH";

  board();

  setcolor(LIGHTCYAN);
  settextstyle(DEFAULT_FONT,HORIZ_DIR,2);
  outtextxy(getmaxx()/2-(textwidth(str)/2),textheight(str)-10,str);
  settextstyle(DEFAULT_FONT,HORIZ_DIR,0);

  setcolor(LIGHTBLUE);
  outtextxy(x0-35,y0-y_range-15,c);
  outtextxy(x0+x_range-30,y0+30,"TIME");

  setcolor(RED);
  for (i=n;i<=row;i+=n)
    line(x0+10+x_temp*(i-1),y0,x0+10+x_temp*(i-1),y0-y_range);

  setcolor(WHITE);
  line(x0,y0-10,x0-4,y0-10);
  gprint(x0-60,y0-15,"%5.1f",min);

  for (j=1;j<=3;j++) {
    line(x0,y0-(y_range-10)/3*j,x0-4,y0-(y_range-10)/3*j);
    gprint(x0-60,y0-(y_range-10)/3*j-2,"%5.1f",min+(max-min)/3*j);
  }

  j=0;
  for (i=n;i<=row;i+=n) {
    line(x0+10+x_temp*(i-1),y0,x0+10+x_temp*(i-1),y0+4);
    if (row/n>10 && j%2==1) {
      gprint(x0+10+x_temp*(i-1)-10,y0+20,"%d",(int)time[i-1]);
    }
    else gprint(x0+10+x_temp*(i-1)-10,y0+10,"%d",(int)time[i-1]);
    j++;
  }

  if (!strcmp(c,"C") || !strcmp(c,"I")) {
    for (i=1;i<row-2;i++) {
      circle((x0+10)+x_temp*i,(y0-10)-y_temp*(data[i]-min),3);
      setfillstyle(SOLID_FILL,LIGHTBLUE);
      floodfill((x0+10)+x_temp*i-1,(y0-10)-y_temp*(data[i]-min)-1,WHITE);
      line((x0+10)+x_temp*i,(y0-10)-y_temp*(data[i]-min),
           (x0+10)+x_temp*(i+1),(y0-10)-y_temp*(data[i+1]-min));
    }
    circle((x0+10)+x_temp*i,(y0-10)-y_temp*(data[i]-min),3);
    setfillstyle(SOLID_FILL,LIGHTBLUE);
    floodfill((x0+10)+x_temp*i-1,(y0-10)-y_temp*(data[i]-min)-1,WHITE);
  }
  else if (strcmp(c,"S")) {
    for (i=0;i<row;i++) {
      circle((x0+10)+x_temp*i,(y0-10)-y_temp*(data[i]-min),3);
      setfillstyle(SOLID_FILL,LIGHTBLUE);
      floodfill((x0+10)+x_temp*i-1,(y0-10)-y_temp*(data[i]-min)-1,WHITE);
      if (i<row-1) line((x0+10)+x_temp*i,(y0-10)-y_temp*(data[i]-min),
                        (x0+10)+x_temp*(i+1),(y0-10)-y_temp*(data[i+1]-min));
    }
  }
  else {
    for (i=0;i<row;i++) {
      circle((x0+10)+x_temp*i,(y0-10)-y_temp*(data[i%n]-min),3);
      setfillstyle(SOLID_FILL,LIGHTBLUE);
      floodfill((x0+10)+x_temp*i-1,(y0-10)-y_temp*(data[i%n]-min)-1,WHITE);
      if (i<row-1)
        line((x0+10)+x_temp*i,(y0-10)-y_temp*(data[i%n]-min),
             (x0+10)+x_temp*(i+1),(y0-10)-y_temp*(data[(i+1)%n]-min));
    }
  }

  getch();
}


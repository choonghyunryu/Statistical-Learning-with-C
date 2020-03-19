/* Program : RESANAL.C ver 0.1   */
/* Author  : Ryu choong hyun     */
/* Date    : 94.7.30.	         */
/* Note    : Residual Analysis   */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <graphics.h>
#include <stdarg.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define PI M_PI
#define FACTOR 50

void get_mean_ss(int col,double *mean,double *ss);
double get_sxy(void);
double get_sse(double b0,double b1);
void get_minmax(int col,double *min,double *max);
double get_maxr(void);
void sort(int left,int right,double *pdata);
void residual(void);
double nor_point(double point);

void nor_plot(void);
void init_graph(void);
void board(void);
int gprint(int x,int y,char *fmt, ...);

int row,x0,y0,x_range,y_range;
double data[FACTOR][2],mean[2],ss[2],e[FACTOR],r[FACTOR],g[FACTOR],
       min[2],max[2],x_temp,y_temp,r_max;

void main(int argc,char *argv[])
{
  FILE *stream;
  int i,j;
  double num,b0,b1,s_xy,h[FACTOR],e_sum1=0,e_sum2=0,
         dw,*rdata,sse,mse;

  clrscr();

  rdata=r;

  if (argc<=1) {
    puts("Usage : resanal datafile");
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
  sse=get_sse(b0,b1);
  mse=sse/(row-2);

  for (i=0;i<row;i++) {
    e[i]=data[i][1]-(b0+b1*data[i][0]);
    if (i!=0) e_sum1+=pow(e[i]-e[i-1],2);
    e_sum2+=pow(e[i],2);
    g[i]=nor_point(((i+1)-1/4.)/(row+3/8.));
    h[i]=1./row+pow(data[i][0]-mean[0],2)/ss[0];
    r[i]=e[i]/(sqrt(mse)*sqrt(1-h[i]));
  }

  dw=e_sum1/e_sum2;

  sort(0,row-1,rdata);

  r_max=get_maxr();

  printf("\t    ** Independent Test **\n\n");
  printf("\tDurbin-Watson Statistic = %lf\n",dw);
  printf("\tFor Number of Object    = %d\n",row);
  if (0<dw && dw<2) printf("\tPositive Antocorrelation\n");
  else if (dw==2) printf("\tNot Antocorrelation\n");
  else printf("\tNegative Antocorrelation\n");

  getch();
  init_graph();
  residual();
  nor_plot();
  closegraph();
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

double get_sse(double b0,double b1)
{
  int i;
  double sum=0;

  for (i=0;i<row;i++) sum+=pow(data[i][1]-(b0+b1*data[i][0]),2);

  return sum;
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

double get_maxr(void)
{
  int i;
  double max=0;

  for (i=0;i<row;i++) if (max<fabs(r[i])) max=fabs(r[i]);

  return max;
}

void sort(int left,int right,double *pdata) /* Quick sort */
{
  int i=left,j=right;
  double temp,mid=*(pdata+((left+right)/2));

  do {
    while (*(pdata+i) < mid && i<right)
      i++;
    while (mid < *(pdata+j) && j>left)
      j--;
    if (i<=j) {
      temp=*(pdata+i);
      *(pdata+i)=*(pdata+j);
      *(pdata+j)=temp;
      i++;
      j--;
    }
  } while (i<=j);
  if (left<j)
    sort(left,j,pdata);
  if (i<right)
    sort(i,right,pdata);
}

double nor_point(double point)
{
  double q,r,px;
  double a[8]={3.3871328727963666080e0,
	       1.3314166789178437745e2,
	       1.9715909503065514427e3,
	       1.3731693765509461125e4,
	       4.5921953931549871457e4,
	       6.7265770927008700853e4,
	       3.3430575583588128105e4,
	       2.5090809287301226727e3};
  double b[7]={4.2313330701600911252e1,
	       6.8718700749205790830e2,
	       5.3941960214247511077e3,
	       2.1213794301586595867e4,
	       3.9307895800092710610e4,
	       2.8729085735721942674e4,
	       5.2264952788528545610e3};
  double c[8]={1.42343711074968357734e0,
	       4.63033784615654529590e0,
	       5.76949722146069140550e0,
	       3.64784832476320460504e0,
	       1.27045825245236838258e0,
	       2.41780725177450611770e-1,
	       2.27238449892691845833e-2,
	       7.74545014278341407640e-4};
  double d[7]={2.05319162663775882187e0,
	       1.67638483018380384940e0,
	       6.89767334985100004550e-1,
	       1.48103976427480074590e-1,
	       1.51986665636164571966e-2,
	       5.47593808499534494600e-4,
	       1.05075007164441684324e-9};
  double e[8]={6.65790464350110377720e0,
	       5.46378491116411436990e0,
	       1.78482653991729133580e0,
	       2.96560571828504891230e-1,
	       2.65321895265761230930e-2,
	       1.24266094738807843860e-3,
	       2.71155556874348757815e-5,
	       2.01033439929228813265e-7};
  double f[7]={5.99832206555887937690e-1,
	       1.36929880922735805310e-1,
	       1.48753612908506148525e-2,
	       7.86869131145613259100e-4,
	       1.84631831751005468180e-5,
	       1.42151175831644588870e-7,
	       2.04426310338993978564e-15};

  q=point-.5;

  if (fabs(q)<=.425) {
    r=.425*.425-q*q;
    px=(((((((a[7]*r+a[6])*r+a[5])*r+a[4])*r+a[3])*r+a[2])*r+a[1])*r+a[0])/
       (((((((b[6]*r+b[5])*r+b[4])*r+b[3])*r+b[2])*r+b[1])*r+b[0])*r+1);
    return (q*px);
  }
  else {
    r=(q<0.) ? point : 1-point;
    if (r<=0) return 0;
    r=sqrt(-log(r));

    if (r<=5.0) {
      r=r-1.6;
      px=(((((((c[7]*r+c[6])*r+c[5])*r+c[4])*r+c[3])*r+c[2])*r+c[1])*r+c[0])/
	 (((((((d[6]*r+d[5])*r+d[4])*r+d[3])*r+d[2])*r+d[1])*r+d[0])*r+1);
    }
    else {
      r=r-5.0;
      px=(((((((e[7]*r+e[6])*r+e[5])*r+e[4])*r+e[3])*r+e[2])*r+e[1])*r+e[0])/
	 (((((((f[6]*r+f[5])*r+f[4])*r+f[3])*r+f[2])*r+f[1])*r+f[0])*r+1);
    }
  }
  return (q<0) ? -px : px;
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

  x_temp=(x_range-10)/(g[row-1]-g[0]);
  y_temp=(y_range-10)/(r[row-1]-r[0]);

  rectangle(0,0,getmaxx(),getmaxy());
  setcolor(LIGHTRED);
  outtextxy(getmaxx()/2-(textwidth(str)/2),getmaxy()-textheight(str)-2,str);
  setcolor(WHITE);
  line(x0,y0,x0,yn);
  line(x0,y0,xn,y0);
  gprint(x0-20,y0+3,"%1d",0);
}

void nor_plot(void)
{
  int i,j;
  char str[]="Normal Probability Plot";

  board();

  setcolor(LIGHTCYAN);
  settextstyle(DEFAULT_FONT,HORIZ_DIR,2);
  outtextxy(getmaxx()/2-(textwidth(str)/2),textheight(str)-10,str);
  settextstyle(DEFAULT_FONT,HORIZ_DIR,0);

  setcolor(LIGHTBLUE);
  outtextxy(x0-55,y0-y_range-15,"Standard Residual");
  outtextxy(x0+x_range-45,y0+30,"Normal Score");

  setcolor(YELLOW);

  line(x0,y0-10,x0-4,y0-10);
  gprint(x0-50,y0-15,"%5.1f",r[0]);
  line(x0+10,y0,x0+10,y0+4);
  gprint(x0-5,y0+10,"%5.1f",g[0]);

  for (j=1;j<=3;j++) {
    line(x0,y0-(y_range-10)/3*j,x0-4,y0-(y_range-10)/3*j);
    gprint(x0-50,y0-(y_range-10)/3*j-2,"%5.1f",r[0]+(r[row-1]-r[0])/3*j);
  }

  for (i=1;i<=5;i++) {
    line(x0+(x_range-10)/5*i,y0,x0+(x_range-10)/5*i,y0+4);
    gprint(x0+(x_range-10)/5*i-23,y0+10,"%5.1f",g[0]+(g[row-1]-g[0])/5*i);
  }

  for (i=0;i<row;i++) {
    circle((x0+10)+x_temp*(g[i]-g[0]),(y0-10)-y_temp*(r[i]-r[0]),3);
    setfillstyle(SOLID_FILL,LIGHTBLUE);
    floodfill((x0+10)+x_temp*(g[i]-g[0]),(y0-10)-y_temp*(r[i]-r[0]),YELLOW);
  }
  getch();
}

void residual(void)
{
  int i,j;
  char str[]="RESIDUAL GRAPH";

  board();

  x_temp=(x_range-10)/(max[0]-min[0]);
  y_temp=(y_range-10)/(max[1]-min[1]);

  setcolor(LIGHTCYAN);
  settextstyle(DEFAULT_FONT,HORIZ_DIR,2);
  outtextxy(getmaxx()/2-(textwidth(str)/2),textheight(str)-10,str);
  settextstyle(DEFAULT_FONT,HORIZ_DIR,0);

  setcolor(LIGHTBLUE);
  outtextxy(x0-45,y0-y_range-15,"Standard Residual");
  outtextxy(x0+x_range-30,y0+30,"X");

  setcolor(YELLOW);

  line(x0+10,y0,x0+10,y0+4);
  gprint(x0-5,y0+10,"%5.1f",min[0]);
  setcolor(WHITE);
  line(x0+10,y0-y_range/2.,x0+x_range,y0-y_range/2.);
  setcolor(YELLOW);
  gprint(x0-10,y0-y_range/2.,"%c",0);

  line(x0,y0-y_range/2.-100,x0-5,y0-y_range/2.-100);
  gprint(x0-55,y0-y_range/2.-103,"%5.1f",r_max);
  line(x0,y0-y_range/2.+100,x0-5,y0-y_range/2.+100);
  gprint(x0-55,y0-y_range/2.+97,"%5.1f",-r_max);

  for (i=1;i<=5;i++) {
    line(x0+(x_range-10)/5*i,y0,x0+(x_range-10)/5*i,y0+4);
    gprint(x0+(x_range-10)/5*i-23,y0+10,"%5.1f",min[0]+(max[0]-min[0])/5*i);
  }

  for (i=0;i<row;i++) {
    circle((x0+10)+x_temp*(data[i][0]-min[0]),
           y0-y_range/2.-(r[i]/r_max*100),3);
    setfillstyle(SOLID_FILL,LIGHTBLUE);
    floodfill((x0+10)+x_temp*(data[i][0]-min[0]),
              y0-y_range/2.-(r[i]/r_max*100),YELLOW);
  }
  getch();
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


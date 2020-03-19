/* Program : NORMAL.C ver 0.2         */
/* Author  : Ryu choong hyun          */
/* Date    : 94.6.20.		      */
/* Note    : Normal Distribtion &     */
/*           Plot.                    */

#include <stdio.h>
#include <stdarg.h>
#include <conio.h>
#include <math.h>
#include <graphics.h>

#define round(x) ((x>0) ? floor(x+.5) : ceil(x-.5))
#define PI M_PI
#define MAX 61

double normal(double z);
double nor_p(double z);
double nor_point(double point);

int gprint(int x,int y,char *fmt, ...);
void init_graph(void);
void end_graph(void);
void board(void);
void curve(double a,double b,double aval,double bval,char c);

char opt;
double pdf[MAX],cdf[13],mean,dev,point;

int x0,y0,x_range,y_range;

void main(void)
{
  char c,job;
  double i,a,b,sa,sb=1000,aval,bval,value,z;

  clrscr();
  printf("\n\t ** Normal Distribution **\n");
  printf("\n\t    Mean = ");
  scanf("%lf",&mean);
  printf("\n\t    S.D. = ");
  scanf("%lf",&dev);
  printf("\n\t1. Probability  2. Percentage Point\n");

  do
  switch (job=getch()) {
    case '1' :
      printf("\n\t1. (a<X<b)  2. (X>a)  3. (X<a)\n");
      do
      switch (c=getch()) {
      case '1' : printf("\n\t    From = ");
		 scanf("%lf",&a);
		 sa=(a-mean)/dev;
		 printf("\n\t    To   = ");
		 scanf("%lf",&b);
		 sb=(b-mean)/dev;
		 value=nor_p(sb)-nor_p(sa);
		 value=(b>a) ? value : -value;
		 clrscr();
		 printf("\t      ** Normal Probability Distribution **\n\n");
		 printf("\t\tMean = %-10.4lf   S.D. = %-10.4lf\n",mean,dev);
		 if (a>b)
		 printf("\t\tP(%lf \x3c X \x3c %lf) = %lf\n\n",b,a,value);
		 else
		 printf("\t\tP(%lf \x3c X \x3c %lf) = %lf\n\n",a,b,value);
		 break;
      case '2' : printf("\n\t    a    = ");
		 scanf("%lf",&a);
		 sa=(a-mean)/dev;
		 value=1-nor_p(sa);
		 b=mean-10*dev;
		 clrscr();
		 printf("\t      ** Normal Probability Distribution **\n\n");
		 printf("\t\tMean = %-10.4lf   S.D. = %-10.4lf\n",mean,dev);
		 printf("\t\tP(%lf \x3c X) = %lf\n\n",a,value);
		 break;
      case '3' : printf("\n\t    a    = ");
		 scanf("%lf",&a);
		 sa=(a-mean)/dev;
		 value=nor_p(sa);
		 b=mean-10*dev;
		 clrscr();
		 printf("\t      ** Normal Probability Distribution **\n\n");
		 printf("\t\tMean = %-10.4lf   S.D. = %-10.4lf\n",mean,dev);
		 printf("\t\tP(X \x3c %lf) = %lf\n\n",a,value);
		 break;
      default :  break;
    } while (!(c=='1' || c=='2' || c=='3'));
    break;
    case '2' : printf("\n\t1. Lower Percentage Point "
		      "2. Upper Percentage Point\n\n");
	       do
	       switch (opt=getch()) {
		 case '1' :
		 case '2' : printf("\n\tPercentage Point : ");
			    scanf("%lf",&point);
			    break;
		 default  : break;
	       } while (!(opt=='1' || opt=='2'));

	       clrscr();
	       printf("\t      ** Normal Probability Distribution **\n\n");
	       printf("\t\tMean = %-10.4lf  S.D. = %-10.4lf\n",mean,dev);
	       if (opt=='1')
		 printf("\t%5.2f Lower Percentage Point = %lf(Z=%lf)\n\n",
			point,nor_point(point)*dev+mean,
			nor_point(point));
	       else
		 printf("\t%5.2f Upper Percentage Point = %lf(Z=%lf)\n\n",
			point,nor_point(1-point)*dev+mean,
			nor_point(1-point));
	       break;
    default  : break;
  } while (!(job=='1' || job=='2'));

  printf("\t      X \t     Z  \t  P( Z<X )\n\n");
  for (i=-3;i<=3;i+=.5)
    printf("\t%10lf\t%10lf\t%10lf\n",i*dev+mean,i,nor_p(i));
  getch();

  for (i=0;i<=60;i++)
    pdf[i]=normal(-3+i/10);

  if (job=='1') {
     aval=normal(sa);
     bval=normal(sb);
  }

  init_graph();

  if (job=='1')  curve(a,b,aval,bval,c);
  else {
    if (opt=='1') {
      aval=normal(nor_point(point));
      curve(nor_point(point)*dev+mean,99999,aval,99999,'3');
    }
    else {
      aval=normal(nor_point(1-point));
      curve(nor_point(1-point)*dev+mean,99999,aval,99999,'2');
    }
  }

  closegraph();
}

double normal(double z)
{
  return 1/(sqrt(2*PI))*exp(-z*z/2);
}

double nor_p(double z)
{
  double x,t,temp;
  double b[5]={ 0.319381530, -.356563782, 1.781477937,
	       -1.821255978, 1.330274429};

  x=fabs(z);
  t=1/(1+.2316419*x);
  temp=b[0]*t+b[1]*pow(t,2)+b[2]*pow(t,3)+b[3]*pow(t,4)+b[4]*pow(t,5);
  return (z<0) ? normal(x)*temp : 1-normal(x)*temp;
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
  int xasp,yasp,graphdrive=DETECT,graphmode;
  initgraph(&graphdrive,&graphmode,"");
}

void board(void)
{
  int xn,yn;
  char str[]="Any key pressed continue !!";
  x0=70;
  y0=getmaxy()-90;
  xn=getmaxx()-40;
  yn=getmaxy()-y0;
  x_range=abs(xn-x0);
  y_range=abs(yn-y0);
  cleardevice();
  rectangle(0,50,getmaxx(),getmaxy()-40);
  setcolor(12);
  outtextxy(getmaxx()/2-(textwidth(str)/2),getmaxy()-textheight(str)-2,str);
  setcolor(15);
  line(x0,y0,x0,yn);
  line(x0,y0,xn,y0);
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

void curve(double a,double b,double aval,double bval,char c)
{
  int i,j;
  double atemp,btemp;
  char str[]="Normal Distribution";

  atemp=(a-(mean-3*dev))*10/dev;
  btemp=(b-(mean-3*dev))*10/dev;

  board();

  setcolor(11);
  settextstyle(0,0,2);
  outtextxy(getmaxx()/2-(textwidth(str)/2),textheight(str)-10,str);
  settextstyle(0,0,0);
  setcolor(10);
  gprint(x0-30,y0-y_range-15,"f(X)");
  gprint(x0+x_range-7,y0+33,"X");
  setcolor(15);

  for (i=1;i<=10;i++) {
    line(x0,y0-y_range*i/10,x0-4,y0-y_range*i/10);
    if (dev>30) gprint(x0-67,y0-y_range*i/10-3,"%.1e",pdf[30]/10*i);
    else gprint(x0-67,y0-y_range*i/10-3,"%7.3f",pdf[30]/10*i);
    for (j=1;j<=50;j++)
      putpixel(x0+x_range*j/50,y0-y_range*i/10,12);
  }

  gprint(x0-60,y0+20,"-\xec");
  gprint(x0+x_range+12,y0+20,"\xec");

   for (i=1;i<=59;i++)
     line(x0+x_range/60*(i),y0-y_range/pdf[30]*pdf[i],
	  x0+x_range/60*(i+1),y0-y_range/pdf[30]*pdf[i+1]);

   line(x0+x_range/60,y0-y_range/pdf[30]*pdf[1],x0+x_range/60,y0);
   line(x0+x_range/60*60,y0-y_range/pdf[30]*pdf[60],x0+x_range/60*60,y0);

   for (i=0;i<=6;i++) {
     if (i==0) {
       line(x0+x_range/60,y0,x0+x_range/60,y0+4);
       gprint(x0+x_range/60-40,y0+8,"%-10.3lf",mean-3*dev);
     }
     else {
       line(x0+x_range/60*(i*10),y0,x0+x_range/60*(i*10),y0+4);
       gprint(x0+x_range/60*(i*10)-30,y0+8,"%-10.3lf",mean+(-3+i)*dev);
     }

     if (i==3) gprint(x0+x_range/60*(i*10)-2,y0+20,"\xe6");
     else if (i<3) gprint(x0+x_range/60*(i*10)-12,y0+20,"\xe6%d\xe5",i-3);
	  else gprint(x0+x_range/60*(i*10)-12,y0+20,"\xe6+%d\xe5",i-3);
  }

  if (mean-3*dev<=a<=mean+3*dev) {
    if (a==mean-3*dev) atemp=1;
    line(x0+x_range/60*atemp,y0,x0+x_range/60*atemp,y0-y_range/pdf[30]*aval);
    gprint(x0+x_range/60*atemp,y0-y_range/pdf[30]*aval-30,"a");
  }

  if (mean-3*dev<=b<=mean+3*dev) {
    if (b==mean-3*dev) btemp=1;
    line(x0+x_range/60*btemp,y0,x0+x_range/60*btemp,y0-y_range/pdf[30]*bval);
    gprint(x0+x_range/60*btemp,y0-y_range/pdf[30]*bval-30,"b");
  }

  setfillstyle(6,BLUE);
  switch (c) {
    case '1' : if ((a>=mean+3*dev && b>=mean+3*dev) ||
		   (a<=mean-3*dev && b<=mean-3*dev)) break;
	       floodfill(x0+x_range/60*atemp+1,y0-1,WHITE);
	       if (a<=mean-3*dev) {
		 floodfill(x0+x_range/60*btemp-1,y0-1,WHITE);
		 if (b>=mean+3*dev) floodfill(x0+x_range/60*30,y0-1,WHITE);
	       }
	       break;
    case '2' : if (a>=mean+3*dev) break;
	       if (a<=mean-3*dev) atemp=30;
	       floodfill(x0+x_range/60*atemp+1,y0-1,WHITE);
	       break;
    case '3' : if (a<=mean-3*dev) break;
	       if (a>=mean+3*dev) atemp=30;
	       floodfill(x0+x_range/60*atemp-1,y0-1,WHITE);
	       break;
    default  : break;
  }

  getch();
}


/* Program : TDIS.C ver 0.2                */
/* Author  : Ryu choong hyun               */
/* Date    : 94.7.5		      	   */
/* Note    : T-Probability Distribution &  */
/*           Poly,etc.                	   */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <conio.h>
#include <math.h>
#include <graphics.h>

#define PI M_PI
#define MAX 61

double fac(int x);
double gamma(int x);
double nor_point(double point);
double tdis(double t);
double t_q(int df,double t);
double t_point(int df,double point);
double t_sub(int df,double p);

int gprint(int x,int y,char *fmt, ...);
void init_graph(void);
void board(void);
void curve(double a,double b,double aval,double bval,char c);

int df;
double mean,dev,pdf[MAX];

int x0,y0,x_range,y_range;

void main(void)
{
  int n,i;
  char c,job,opt;
  double value,ax,bx,a,b=1000,aval,bval,cdf[13],point;

  clrscr();
  printf("\n\t ** T-Probability Distribution **\n");
  printf("\n\t1. Sample Probability  2. DF & T Value  3.Percentage Point\n");

  do
  switch (job=getch()) {
    case '1' : printf("\n\t             MEAN      =  ");
	       scanf("%lf",&mean);
	       printf("\n\t             S.Dev     =  ");
	       scanf("%lf",&dev);
	       printf("\n\t             Sample    =  ");
	       scanf("%d",&n);
	       df=n-1;
	       printf("\n\t               << Interval >>\n");
	       printf("\t        1. a<x<b   2. x>a   3.x<a\n");
	       do
	       switch (c=getch()) {
		 case '1' : printf("\n\t             FROM   = ");
			    scanf("%lf",&ax);
			    printf("\n\t             TO     = ");
			    scanf("%lf",&bx);
			    a=(ax-mean)/(dev/sqrt(n));
			    b=(bx-mean)/(dev/sqrt(n));
			    value=(1-t_q(df,b))-(1-t_q(df,a));
			    clrscr();
			    printf("\n\n\t\t** T-Probability Distribution **"
				   "\n\n");
			    printf("\t\t\tMean   = %-10.4lf\n\t\t\tSample = "
				   "%-3d\n\t\t\tS.Dev  = %-10.4lf",mean,n,dev);
			    printf("\n\t\t\tA      = %-10.4lf\n\t\t\tB      ="
				   " %-10.4lf",ax,bx);
			    printf("\n\t      P(%lf < T < %lf) = %lf\n",
				   a,b,value);
			    break;
		 case '2' : printf("\n\t             a    = ");
			    scanf("%lf",&ax);
			    a=(ax-mean)/(dev/sqrt(n));
			    b=mean-10*dev;
			    if (ax==mean) value=.5;
			    else value=t_q(df,a);
			    clrscr();
			    printf("\n\n\t\t** T-Probability Distribution **"
				   "\n\n");
			    printf("\t\t\tMean   = %-10.4lf\n\t\t\tSample = "
				   "%-3d\n\t\t\tS.Dev  = %-10.4lf",mean,n,dev);
			    printf("\n\t\t\tA      = %-10.4lf",ax);
			    printf("\n\t\t   P(T > %lf) = %lf\n",a,value);
			    break;
		 case '3' : printf("\n\t             a    = ");
			    scanf("%lf",&ax);
			    a=(ax-mean)/(dev/sqrt(n));
			    b=mean-10*dev;
			    if (ax==mean) value=.5;
			    else value=1-t_q(df,a);
			    clrscr();
			    printf("\n\n\t\t** T-Probability Distribution **"
				   "\n\n");
			    printf("\t\t\tMean   = %-10.4lf\n\t\t\tSample = "
				   "%-3d\n\t\t\tS.Dev  = %-10.4lf",mean,n,dev);
			    printf("\n\t\t\tA      = %-10.4lf",ax);
			    printf("\n\t\t   P(T < %lf) = %lf\n",a,value);
			    break;
		 default  : break;
	       } while(!(c=='1' || c=='2' || c=='3'));
	       break;
    case '2' : mean=0;
	       dev=1;
	       printf("\n\t       Degree of Freedom = ");
	       scanf("%d",&df);
	       printf("\n\t               << Interval >>\n");
	       printf("\t        1. a<t<b   2. t>a   3.t<a\n");
	       do
	       switch (c=getch()) {
		 case '1' : printf("\n\t             FROM   = ");
			    scanf("%lf",&a);
			    printf("\n\t             TO     = ");
			    scanf("%lf",&b);
			    value=fabs(t_q(df,b)-t_q(df,a));
			    clrscr();
			    printf("\n\n\t\t** T-Probability Distribution **"
				   "\n\n");
			    printf("\n\t\t\tA      = %-10.4lf\n\t\t\tB      ="
				   " %-10.4lf",a,b);
			    printf("\n\t      P(%lf < T < %lf) = %lf\n",
				   a,b,value);
			    break;
		 case '2' : printf("\n\t             T    = ");
			    scanf("%lf",&a);
			    b=99999;
			    value=t_q(df,a);
			    clrscr();
			    printf("\n\n\t\t** T-Probability Distribution **"
				   "\n\n");
			    printf("\n\t\t\tA      = %-10.4lf",a);
			    printf("\n\t\t   P(T > %lf) = %lf\n",a,value);
			    break;
		 case '3' : printf("\n\t             T    = ");
			    scanf("%lf",&a);
			    b=99999;
			    value=1-t_q(df,a);
			    clrscr();
			    printf("\n\n\t\t** T-Probability Distribution **"
				   "\n\n");
			    printf("\n\t\t\tA      = %-10.4lf",a);
			    printf("\n\t\t   P(T < %lf) = %lf\n",a,value);
			    break;
		 default  : break;
	       } while(!(c=='1' || c=='2' || c=='3'));
	       break;
    case '3' : mean=0;
	       dev=1;
	       printf("\n\t      Degree of Freedom = ");
	       scanf("%d",&df);
	       printf("\n\t  1. Lower Percentage Point"
		      "  2. Upper Percentage Point\n");
	       do
	       switch (opt=getch()) {
		 case '1' :
		 default  : printf("\n\n\tPercentage Point = ");
			    scanf("%lf",&point);
			    break;
	       } while (!(opt=='1' || opt=='2'));

	       clrscr();
	       printf("\t     ** T Probability Distribution **\n\n");
	       b=99999;
	       if (opt!='1') {
		 if (point>=.5) a=-t_point(df,1-point);
		 else a=t_point(df,point);
		 printf("\t\t%5.4f Upper Percentage Point = %lf\n\n",point,a);
	       }
	       else {
		 if (point<=.5) a=-t_point(df,point);
		 else a=t_point(df,1-point);
		 printf("\t\t%5.4f Lower Percentage Point = %lf\n\n",point,a);
	       }
	       if (opt=='1') opt+=2;
	       break;
    default  : break;
  } while(!(job=='1' || job=='2' || job=='3'));

  for (i=0;i<13;i++) {
    if (i<=6) cdf[i]=1-t_q(df,.5*(i-6));
    else cdf[i]=.5+(.5-cdf[12-i]);
  }
  printf("\n\t            X           t         P( T<t )\n");
  for (i=0;i<13;i++)
    if (job=='1') printf("\t      %10.4f  %10.4f    %10.7f\n",
			 mean+(.5*(i-6)*dev)/sqrt(n),.5*(i-6),cdf[i]);
    else printf("\t      %10.4f  %10.4f    %10.7f\n",
		.5*(i-6)/sqrt(df+1),.5*(i-6),cdf[i]);
  getch();

    for (i=0;i<=60;i++)
      pdf[i]=tdis(-3+i/10.);

  aval=tdis(a);
  bval=tdis(b);

  init_graph();

  switch (job) {
    case '1' :
    case '2' : curve(a,b,aval,bval,c);
	       break;
    case '3' : curve(a,b,aval,bval,opt);
	       break;
  }
}

double fac(int x)
{
  return (x==0 ? 1 : x*fac(x-1));
}

double gamma(int x)
{
  int i;
  double sum=1;

  if (x<=0) {
    fprintf(stderr,"Gamma-Function's Value <=0 !!!");
    exit(0);
  }

  x*=2;

  if (x % 2==0) return fac(x/2-1);
  else {
    if (x==1) return sqrt(PI);
    else
      for (i=1;i<=x/2;i++) sum*=(i*2-1);
      return sum*sqrt(PI)/pow(2,(int)x/2);
  }
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

double tdis(double t)
{
  return gamma((df+1)/2.)/(sqrt((double)df*PI)*gamma(df/2.))
	 *pow(1+t*t/df,-((double)df+1)/2);
}

double t_q(int df,double t)
{
  int i;
  double theta,s,c,sum=0,temp=1,p,tx;

  tx=t;
  t=(t<0) ? -t : t;

  theta=atan(t/sqrt(df));
  c=sqrt(df/(df+t*t));
  s=sqrt(1-df/(df+t*t));

  for (i=df%2+4;i<=df;i+=2) {
    temp*=(i-3.)/(i-2.);
    sum+=temp*s*pow(c,i-2);
  }

  if (df & 1) {
    if (df==1) p=1-(.5+theta/PI);
    else p=1-(.5+(theta+s*c+sum)/PI);
  }
  else p=1-(.5+(s+sum)/2);

  return (tx<0) ? 1-p : p;
}

double t_point(int df,double point)
{
  double u,y[5];

  if(df>5 && (point>5.e-3 || df>13)) {
  u=nor_point(1-point);
  y[0]=(pow(u,3)+u)/pow(2,2);
  y[1]=(5*pow(u,5)+16*pow(u,3)+3*u)/(pow(2,5)*3);
  y[2]=(3*pow(u,7)+19*pow(u,5)+17*pow(u,3)-15*u)/(pow(2,7)*3);
  y[3]=(79*pow(u,9)+776*pow(u,7)+1482*pow(u,5)-1920*pow(u,3)-945*u)/
       (pow(2,11)*pow(3,2)*5);
  y[4]=(27*pow(u,11)+339*pow(u,9)+930*pow(u,7)-1782*pow(u,5)-765*pow(u,3)+
	17955*u)/(pow(2,13)*pow(3,2)*5);

  return u+y[0]/df+y[1]/pow(df,2)+y[2]/pow(df,3)+y[3]/pow(df,4)+
	 y[4]/pow(df,5);
  }
  else return t_sub(df,point);
}

double t_sub(int df,double p)
{
  double  eps,q,point=0,s=10000;

  if (df==1 && p<0.01 && p>=0.001) eps = 1.0e-4;
  else if(df==2 && p< 0.0001)      eps = 1.0e-4;
  else if(df==1 && p< 0.001)       eps = 1.0e-2;
  else                             eps = 1.0e-5;

  while (1) {
    point+=s;
    if (s<=eps) return point;
    q=t_q(df,point)-p;
    if (q==0) return point;
    if (q<0) {
      point-=s;
      s/=10;
    }
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
  char str[]="T Distribution";

  atemp=(a-(mean-3*dev))*10/dev;
  btemp=(b-(mean-3*dev))*10/dev;

  board();

  setcolor(11);
  settextstyle(0,0,2);
  outtextxy(getmaxx()/2-(textwidth(str)/2),textheight(str)-10,str);
  settextstyle(0,0,0);
  setcolor(10);
  gprint(x0-30,y0-y_range-15,"f(t)");
  gprint(x0+x_range-7,y0+33,"t");
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
  }

  if (mean-3*dev<=a && a<=mean+3*dev) {
    if (a==mean-3*dev) atemp=1;
    line(x0+x_range/60*atemp,y0,x0+x_range/60*atemp,y0-y_range/pdf[30]*aval);
    gprint(x0+x_range/60*atemp,y0-y_range/pdf[30]*aval-30,"a");
  }

  if (mean-3*dev<=b && b<=mean+3*dev) {
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
  closegraph();
}


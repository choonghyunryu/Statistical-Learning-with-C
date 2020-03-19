/* Program : FDIS.C ver 0.1                */
/* Author  : Ryu choong hyun               */
/* Date    : 94.7.6		      	   */
/* Note    : F-robability Distribution &   */
/*           Poly,etc.                	   */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <conio.h>
#include <stdarg.h>
#include <graphics.h>

#define PI M_PI
#define MAX 101

double fac(int x);
double gamma(float x);
double fdis(double f);
double f_q(int df1,int df2,double f);
double f_point(int df1,int df2,double p);
double nor_point(double p);
double normal(double z);
double nor_cdf(double z);

int gprint(int x,int y,char *fmt, ...);
void init_graph(void);
void board(void);
void curve(double a,double b,double aval,double bval,char c);

int df1,df2;
double pdf[MAX],max;

int x0,y0,x_range,y_range;

void main(void)
{
  int i;
  char job,c,opt;
  double var,avar,bvar,value,a,b,avalue,bvalue,point;

  clrscr();
  printf("\n\t ** F Probability Distribution **\n\n");
  printf("\t  1. DF & f  2. Percentage Point\n");

  job=getch();

  printf("\n\tDegree of Freedom 1 = ");
  scanf("%d",&df1);
  printf("\n\tDegree of Freedom 2 = ");
  scanf("%d",&df2);

  do
  switch (job) {
    case '1' : printf("\n\t               << Interval >>\n");
	       printf("\t        1. a< f <b   2. f>a   3.f<a\n");
	       do
	       switch (c=getch()) {
		 case '1' : printf("\n\t             f-a   = ");
			    scanf("%lf",&a);
			    printf("\n\t             f-b   = ");
			    scanf("%lf",&b);
			    value=f_q(df1,df2,a)-f_q(df1,df2,b);
			    clrscr();
			    printf("\n\n\t\t ** F Probability "
				   "Distribution **\n");
			    printf("\n\n\t\t\t P(%lf < f < %lf) = %lf\n",
				   a,b,value);
			    break;
		  case '2': printf("\n\t             f-a   = ");
			    scanf("%lf",&a);
			    b=99999;
			    value=f_q(df1,df2,a);
			    clrscr();
			    printf("\n\n\t\t  ** F Probability "
				   "Distribution **");
			    printf("\n\n\t\t\t P(f > %lf) = %lf\n",a,value);
			    break;
		  case '3': printf("\n\t             f-a   = ");
			    scanf("%lf",&a);
			    b=99999;
			    value=1-f_q(df1,df2,a);
			    clrscr();
			    printf("\n\n\t\t  ** F Probability "
				   "Distribution **");
			    printf("\n\n\t\t\t P(f< %lf) = %lf\n",a,value);
			    break;
		 default  : break;
	       } while(!(c=='1' || c=='2' || c=='3'));
	       break;
    case '2' : printf("\n\t1. Lower Percentage Point "
		      "2. Upper Percentage Point\n");
	       do
	       switch (opt=getch()) {
		 case '1' :
		 default  : printf("\n\n\tPercentage Point = ");
			    scanf("%lf",&point);
			    break;
	       } while (!(opt=='1' || opt=='2'));

	       clrscr();
	       printf("\t     ** F Probability Distribution **\n\n");
	       b=-100;
	       if (opt!='1') {
		 a=f_point(df1,df2,point);
		 printf("\t\t%5.4f Upper Percentage Point= %lf\n\n",point,a);
	       }
	       else {
		 a=f_point(df1,df2,1-point);
		 printf("\t\t%5.4f Lower Percentage Point= %lf\n\n",point,a);
	       }
	       if (opt=='1') opt+=2;
	       break;
    default  : break;
  } while (!(job=='1' || job=='2'));

  getch();

  avalue=fdis(a);
  if (job=='2') ;
  else  bvalue=fdis(b);

  for (i=0;i<101;i++) {
    pdf[i]=fdis(.03*i);
    if (pdf[i]>max)
      max=pdf[i];
  }

  init_graph();
  if (job!='2') curve(a,b,avalue,bvalue,c);
  else curve(a,b,avalue,bvalue,opt);
}

double f_q(int df1,int df2,double f)
{
  int i;
  double theta,s,c,sum=0,temp=1;

  theta=atan(sqrt(df1*f/df2));
  c=sqrt(1/(1+df1*f/df2));
  s=sqrt(1-c*c);

  if ((df1*df2)%2==0 && df1 & 1) return 1-f_q(df2,df1,1/f);

  if ((df1*df2)%2==0) {
    for (i=4;i<=df1;i+=2) {
      temp*=(df2+i-4.)/(i-2.);
      sum+=temp*pow(s,i-2);
    }
    return pow(c,df2)*(1+sum);
  }
  else {
    for (i=3;i<=df2;i+=2) {
      if (i != 3) temp*=((i-3.)/(i-2.));
      sum+=temp*s*pow(c,i-2);
    }
    for (i=3;i<=df1;i+=2) {
      if (df2==1 && i==3) temp=1;
      else temp*=(i+df2-4.)/(i-2.);
      sum-=temp*pow(c,df2)*pow(s,i-2);
    }
    return 1-2*(theta+sum)/PI;
  }
}

double fac(int x)
{
  return (x==0 ? 1 : x*fac(x-1));
}

double gamma(float x)
{
  int i;
  double sum=1;

  if (x<=0) {
    fprintf(stderr,"Gamma-Function's Value <=0 !!!");
    exit(0);
  }

  x*=2;

  if ((int)x % 2==0) return fac(x/2-1);
  else {
    if (x==1) return sqrt(PI);
    else for (i=1;i<=x/2;i++) sum*=(i*2-1);
    return sum*sqrt(PI)/pow(2,(int)x/2);
  }
}

double fdis(double f)
{
  double a,b,c,temp;

  if (f<=0) return 0;
  a=df1/2.;
  b=df2/2.;
  c=a/b;

  temp=gamma(a+b)*pow(c,a)/(gamma(a)*gamma(b));
  return temp*pow(f,a-1)/pow(1+c*f,a+b);
}

double f_point(int df1,int df2,double p)
{
  double gamma,h,delta,u,z,t1,t2,t3,t4,eps,s,fw,qe;

  if (df1<=10 || df2<=10) {
    eps=1.0e-5;
    if (df2==1) eps=1.0e-4;
    s=1000;
    fw=0;

    while (1) {
      fw+=s;
      if (s<=eps) return fw;
      qe=f_q(df1,df2,fw)-p;
      if (qe==0) return fw;
      if (qe<0) {
	fw-=s;
	s/=10;
      }
    }
  }
  else {
    gamma=1./df1+1./df2;
    delta=1./df1-1./df2;
    h=2/gamma;
    u=nor_point(1-p);

    t1=(u*u+2)*delta/6.;
    t2=((pow(u,3)+3*u)/(12*h)+(pow(u,3)+11*u)/144.*h*delta*delta)/sqrt(h);
    t3=(pow(u,4)+9*u*u+8)*delta/(60*h)-(3*pow(u,4)+7*u*u-16)*h*pow(delta,3)
       /6480.;
    t4=((pow(u,5)+20*pow(u,3)+15*u)/(480*h*h)+(pow(u,5)+44*pow(u,3)+183*u)*
	 pow(delta,4)/2880.+(9*pow(u,5)-284*pow(u,3)-1513*u)*h*h*pow(delta,4)
	 /622080.)/sqrt(h);
    z=u/sqrt(h)-t1+t2-t3+t4;

    return exp(2*z);
  }
}

double nor_point(double p)
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

  q=p-.5;

  if (fabs(q)<=.425) {
    r=.425*.425-q*q;
    px=(((((((a[7]*r+a[6])*r+a[5])*r+a[4])*r+a[3])*r+a[2])*r+a[1])*r+a[0])/
       (((((((b[6]*r+b[5])*r+b[4])*r+b[3])*r+b[2])*r+b[1])*r+b[0])*r+1);
    return (q*px);
  }
  else {
    r=(q<0.) ? p : 1-p;
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

double normal(double z)
{
  return 1/(sqrt(2*PI))*exp(-z*z/2);
}

double nor_cdf(double z)
{
  double x,t,temp;
  double b[5]={ 0.319381530, -.356563782, 1.781477937,
	       -1.821255978, 1.330274429};

  x=fabs(z);
  t=1/(1+.2316419*x);
  temp=b[0]*t+b[1]*pow(t,2)+b[2]*pow(t,3)+b[3]*pow(t,4)+b[4]*pow(t,5);
  return (z<0) ? normal(x)*temp : 1-normal(x)*temp;
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

  char str[]="F Distribution";

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
    gprint(x0-67,y0-y_range*i/10-3,"%lf",max/10*i);
    for (j=1;j<=50;j++)
      putpixel(x0+x_range*j/50,y0-y_range*i/10,12);
  }

  for (i=1;i<=5;i++) {
    line(x0+x_range/100*20*i,y0,x0+x_range/100*20*i,y0+4);
    gprint(x0+x_range/100*20*i-25,y0+10,"%-6.3lf",3/5.*i);
  }
  gprint(x0-5,y0+10,"0");
  gprint(x0+x_range+15,y0+10,"\xec");

  if (max<1)
    for (i=0;i<=99;i++)
      line(x0+x_range/100*(i),y0-y_range/max*pdf[i],
	   x0+x_range/100*(i+1),y0-y_range/max*pdf[i+1]);
  else {
    gprint(x0-50,y0-y_range-15,"\xec");
    for (i=1;i<=99;i++)
    line(x0+x_range/100*(i),y0-y_range/max*pdf[i],
	 x0+x_range/100*(i+1),y0-y_range/max*pdf[i+1]);
  }
  if (0<a && a<=3) {
    line(x0+x_range/100*(a*100/3.),y0,x0+x_range/100*(a*100/3.),
	 y0-y_range/max*aval);
    gprint(x0+x_range/100*(a*100/3.)-3,y0-y_range/max*aval-30,"a");
  }

  if (0<b && b<=3) {
    line(x0+x_range/100*(b*100/3.),y0,
	 x0+x_range/100*(b*100/3.),y0-y_range/max*bval);
    gprint(x0+x_range/100*(b*100/3.)-3,y0-y_range/max*bval-30,"b");
  }

  setfillstyle(6,BLUE);
  switch (c) {
    case '1' : if (b>3) {
		 if (a>3) ;
		 else {
		   line(x0+x_range/100*100,y0-y_range/max*pdf[99],
			x0+x_range/100*100,y0);
		   floodfill(x0+x_range-90,y0-1,WHITE);
		 }
	       }
	       else floodfill(x0+x_range/100*(b*100/3.)-1,y0-1,WHITE);
	       break;
    case '2' : line(x0+x_range/100*100,y0-y_range/max*pdf[100],
		    x0+x_range/100*100,y0);
	       if (!(a>=3))
		 floodfill(x0+x_range/100*(a*100/3.)+1,y0-1,WHITE);
	       break;
    case '3' : if (max>1) line(x0,y0-y_range,x0+4,y0-y_range);
	       if (a>=3) {
		 line(x0+x_range/100*100,y0-y_range/max*pdf[100],
		      x0+x_range/100*100,y0);
		 floodfill(x0+x_range-90,y0-1,WHITE);
	       }
	       else floodfill(x0+x_range/100*(a*100/3.)-1,y0-1,WHITE);
	       break;
    default  : break;
  }

  getch();
  closegraph();
}


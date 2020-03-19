/* Program : CSQUARE.C ver 0.2                      */
/* Author  : Ryu choong hyun                        */
/* Date    : 94.6.28.		      	            */
/* Note    : Chi-Square Probability Distribution &  */
/*           Curve,etc.                		    */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <conio.h>
#include <math.h>
#include <graphics.h>

#define PI M_PI
#define MAX 101

double normal(double z);
double chi(double k);
double nor_cdf(double z);
double nor_point(double point);
double chi_q(int df,double chi2);
double chi_point(int df,double point);
double fac(int x);
double gamma(int x);

int gprint(int x,int y,char *fmt, ...);
void init_graph(void);
void board(void);
void curve(double a,double b,double aval,double bval,char c);

int df;
double pdf[MAX],max;

int x0,y0,x_range,y_range;

void main(void)
{
  int n,i;
  char job,c,opt;
  double var,avar,bvar,value,a,b,avalue,bvalue,point;

  clrscr();
  printf("\n\t ** Chi Square Probability Distribution **\n\n");
  printf("\t1. Sample Probability  2. DF & X\xfd  "
	 "3. Percentage Point\n");
  do
  switch (job=getch()) {
    case '2' : printf("\n\t Degree of Freedom = ");
	       scanf("%d",&df);
	       n=df+1;
    case '1' : if (job!='2') {
		 printf("\n\t             VAR       = ");
		 scanf("%lf",&var);
		 printf("\n\t             Sample    = ");
		 scanf("%d",&n);
		 df=n-1;
	       }
	       printf("\n\t               << Interval >>\n");
	       printf("\t        1. a<X\xfd<b   2. X\xfd>a   3.X\xfd<a\n");
	       do
	       switch (c=getch()) {
		 case '1' : if (job!='2') {
			      printf("\n\t             S.VAR a    = ");
			      scanf("%lf",&avar);
			      printf("\n\t             S.VAR b    = ");
			      scanf("%lf",&bvar);
			      a=df*avar/var;
			      b=df*bvar/var;
			    }
			    else {
			      printf("\n\t             X\xfd a   = ");
			      scanf("%lf",&a);
			      printf("\n\t             X\xfd b   = ");
			      scanf("%lf",&b);
			    }
			    value=chi_q(df,a)-chi_q(df,b);
			    clrscr();
			    printf("\n\n\t\t ** Chi Square Probability "
				   "Distribution **\n");
			    printf("\n\n\t\t\t P(%lf < X\xfd < %lf) = %lf\n",
				   a,b,value);
			    break;
		  case '2': if (job!='2') {
			      printf("\n\t             S.VAR a    = ");
			      scanf("%lf",&avar);
			      a=df*avar/var;
			    }
			    else {
			      printf("\n\t             X\xfd a   = ");
			      scanf("%lf",&a);
			    }
			    b=3*n;
			    value=chi_q(df,a);
			    clrscr();
			    printf("\n\n\t\t  ** Chi Square Probability "
				   "Distribution **");
			    printf("\n\n\t\t\t P(X\xfd > %lf) = %lf\n",
				   a,value);
			    break;
		  case '3': if (job!='2') {
			      printf("\n\t             S.VAR a    = ");
			      scanf("%lf",&avar);
			      a=df*avar/var;
			    }
			    else {
			      printf("\n\t             X\xfd a   = ");
			      scanf("%lf",&a);
			    }
			    b=3*n;
			    value=1-chi_q(df,a);
			    clrscr();
			    printf("\n\n\t\t  ** Chi Square Probability "
				   "Distribution **");
			    printf("\n\n\t\t\t P(X\xfd < %lf) = %lf\n",
				   a,value);
			    break;
		 default  : break;
	       } while(!(c=='1' || c=='2' || c=='3'));
	       break;

    case '3' : printf("\n\tDegree of Freedom = ");
	       scanf("%d",&df);
	       n=df+1;
	       printf("\n\t1. Lower Percentage Point "
		      "2. Upper Percentage Point\n");
	       do
	       switch (opt=getch()) {
		 case '1' :
		 default  : printf("\n\n\tPercentage Point = ");
			    scanf("%lf",&point);
			    break;
	       } while (!(opt=='1' || opt=='2'));

	       clrscr();
	       printf("\t     ** Chi Square Probability Distribution **\n\n");
	       b=-10;
	       if (opt!='1') {
		 a=chi_point(df,point);
		 printf("\t\t%5.4f Lower Percentage Point= %lf\n\n",point,a);
	       }
	       else {
		 a=chi_point(df,1-point);
		 printf("\t\t%5.4f Upper Percentage Point= %lf\n\n",point,a);
	       }
	       if (opt=='1') opt+=2;
	       break;
    default  : break;
  } while (!(job=='1' || job=='2' || job=='3'));

  getch();

  avalue=chi(a);
  if (job=='3') ;
  else  bvalue=chi(b);

  for (i=0;i<101;i++) {
    pdf[i]=chi(2.*n/100*i);
    if (pdf[i]>max)
      max=pdf[i];
  }

  init_graph();
  if (job!='3') curve(a,b,avalue,bvalue,c);
  else curve(a,b,avalue,bvalue,opt);
}

double normal(double z)
{
  return 1/(sqrt(2*PI))*exp(-z*z/2);
}

double chi(double x)
{
  if (df==1)
    return 1/(gamma(1/2.)*sqrt(2))*exp(-x/2)*sqrt(x)/10;
  else
    if (df==2)
      return exp(-x/2)/2;
    else
      return 1/(gamma(df/2.)*pow(2,(double)df/2))*exp(-x/2)
	     *pow(x,(double)df/2-1);
}

double nor_p(double z)
{
  double d[6]={.0498673470,.0211410061,.0032776263,
	       .0000380036,.0000488906,.0000053830};
  double x,px;

  x=(z<0) ? -z : z;
  px=1-.5*pow(1+d[0]*x+d[1]*pow(x,2)+d[2]*pow(x,3)+d[3]*pow(x,4)+
	      d[4]*pow(x,5)+d[5]*pow(x,6),-16);

  return (z<0) ? 1-px : px;
}

double chi_q(int df,double chi2)
{
  int i;
  double chi,sum=0,dev=1;

  chi=sqrt(chi2);

  if (df%2==1) {
    for (i=1;i<=(df-1)/2;i++) {
      dev*=i*2-1;
      sum+=pow(chi,2*i-1)/dev;
    }
    return 2*(1-nor_p(chi))+2*normal(chi)*sum;
  }
  else {
    for (i=1;i<=(df-1)/2;i++) {
      dev*=2*i;
      sum+=pow(chi,2*i)/dev;
    }
    return sqrt(2*PI)*normal(chi)*(1+sum);
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

double chi_point(int df, double point)
{
  int i;
  double u,c[11],temp=0;

  u=nor_point(1-point);

  switch (df) {
    case 1  : return pow(nor_point(point*.5),2);
    case 2  : return -2*log(point);
    default :
  c[0]=1;
  c[1]=u;
  c[2]=(pow(u,2)-1)/3;
  c[3]=(pow(u,3)-7*u)/(pow(2,2)*pow(3,2));
  c[4]=(-3*pow(u,4)-7*pow(u,2)+16)/(2*pow(3,4)*5);
  c[5]=(9*pow(u,5)+256*pow(u,3)-433*u)/(pow(2,5)*pow(3,5)*5);
  c[6]=(12*pow(u,6)-243*pow(u,4)-923*pow(u,2)+1472)/(pow(2,3)*pow(3,6)*5*7);
  c[7]=(-3753*pow(u,7)-4353*pow(u,5)+289517*pow(u,3)+289717*u)/
       (pow(2,7)*pow(3,8)*pow(5,2)*7);
  c[8]=(270*pow(u,8)+4614*pow(u,6)-9513*pow(u,4)-104989*pow(u,2)+35968)/
       (pow(2,4)*pow(3,9)*pow(5,2)*7);
  c[9]=(-5139*pow(u,9)-547848*pow(u,7)-2742210*pow(u,5)+7016224*pow(u,3)+
	37501325*u)/(pow(2,11)*pow(3,10)*pow(5,2)*7);
  c[10]=(-364176*pow(u,10)+6208146*pow(u,8)+125735778*pow(u,6)+303753831*
	 pow(u,4)-672186949*pow(u,2)-2432820224)/
	(pow(2,7)*pow(3,13)*pow(5,3)*7*11);
  c[11]=(199112985*pow(u,11)+1885396761*pow(u,9)-31857434154*pow(u,7)-
	 287542736226*pow(u,5)-556030221167*pow(u,3)+487855454729*u)/
	(pow(2,13)*pow(3,14)*pow(5,3)*pow(7,2)*11);

  for (i=0;i<=11;i++)
    temp+=c[i]*pow(2./df,i/2.);
  return df*temp;
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

  char str[]="X\xfd Distribution";

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
    gprint(x0+x_range/100*20*i-25,y0+10,"%-6.3lf",(double)2*(df+1)/5*i);
  }
  gprint(x0-5,y0+10,"0");
  gprint(x0+x_range+15,y0+10,"\xec");

  for (i=0;i<=99;i++)
    line(x0+x_range/100*(i),y0-y_range/max*pdf[i],
	 x0+x_range/100*(i+1),y0-y_range/max*pdf[i+1]);

  if (0<a && a<=2*(df+1)) {
    line(x0+x_range/100*(a*100/(2*(df+1))),y0,
	 x0+x_range/100*(a*100/(2*(df+1))),y0-y_range/max*aval);
    gprint(x0+x_range/100*(a*100/(2*(df+1)))-3,y0-y_range/max*aval-30,"a");
  }

  if (0<b && b<=2*(df+1)) {
    line(x0+x_range/100*(b*100/(2*(df+1))),y0,
	 x0+x_range/100*(b*100/(2*(df+1))),y0-y_range/max*bval);
    gprint(x0+x_range/100*(b*100/(2*(df+1)))-3,y0-y_range/max*bval-30,"b");
  }

  setfillstyle(6,BLUE);
  switch (c) {
    case '1' : if (b>2*(df+1)) {
		 if (a>2*(df+1)) ;
		 else {
		   line(x0+x_range/100*100,y0-y_range/max*pdf[99],
			x0+x_range/100*100,y0);
		   floodfill(x0+x_range-90,y0-1,WHITE);
		 }
	       }
	       else floodfill(x0+x_range/100*(b*100/(2*(df+1)))-1,y0-1,WHITE);
	       break;
    case '2' : line(x0+x_range/100*100,y0-y_range/max*pdf[100],
		    x0+x_range/100*100,y0);
	       if (!a>=2*(df+1))
	         floodfill(x0+x_range/100*(a*100/(2*(df+1)))+1,y0-1,WHITE);
	       break;
    case '3' : if (a>=2*(df+1)) {
		 line(x0+x_range/100*100,y0-y_range/max*pdf[100],
		      x0+x_range/100*100,y0);
		 floodfill(x0+x_range-90,y0-1,WHITE);
	       }
	       else floodfill(x0+x_range/100*(a*100/(2*(df+1)))-1,y0-1,WHITE);
	       break;
    default  : break;
  }

  getch();
  closegraph();
}


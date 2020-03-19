/* Program : STAT.C  ver 0.3          */
/* Author  : Ryu choong hyun	      */
/* Date    : 94.2.25.                 */
/* Note    : Frequency Table &        */
/*           Histogram & Ploy &       */
/*	     Basic Statistic          */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <stdarg.h> /* for gprint() */
#include <math.h>
#include <graphics.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define round(x) ((x>0) ? floor(x+.5) : ceil(x-.5)) /* round function */
#define frac(x) (x-floor(x))
#define PI M_PI
#define LOG2 0.301
#define MAX 5000 /* Maximum data */
#define CLASS 35
#define LOOF_MAX 3 /* for limt unit loof */

double unit(double x);
double g_mean(void);
double h_mean(void);
double q_mean(void);
double get_median(void);
double t_mean(double per);
double win_mean(void);
double sum_square(double mean);
double mean_deviation(double mean);
double per_range(void);
double get_mode(void);
double mx_temp(double mean,int k);
void prn_data(void);
void sort(int left,int right);
void make_table(int class);
void get_gcmdata(void);
void basic_stat(void);
void get_quart(void);

int gprint(int x,int y,char *fmt, ...);
void init_graph(void);
void board(void);
void title(char *title_name);
void fhisto(void);
void chisto(void);
void fpoly(void);
void cpoly(void);

int class_n[CLASS];
double data[MAX],max,min,range,total,class_interval,class_mlimit,
       q1,q3,*pdata,mode_temp,mode_class,max_class;

static int n,class,gcm_data;

int x0,y0,x_temp,x_range,y_range;

void main(int argc,char *argv[])
{
  FILE *stream;
  double num,small_unit,u[LOOF_MAX],range;
  int i;

  clrscr();
  pdata=data;

  if (argc<=1) {
    puts("Usage : stat datafile");
    exit(EXIT_FAILURE);
  }
  stream=fopen(argv[1],"rt");
  if (stream==NULL) printerrmsg("File not found !!");
  while (fscanf(stream,"%lf\n",&num)!=EOF) {
    ++n;
    *(pdata+n)=num;
    total+=*(pdata+n);
  }

  prn_data();
  sort(1,n);
  get_gcmdata();
  class=round(log10(n)/LOG2+1); /* 계급의 갯수를 Sturegs 공식으로 계산 */

  for (i=0;i<LOOF_MAX;i++) u[i]=unit(*(pdata+i+1));
  /* data[1] 부터 data[LOOF_MAX]까지 3개의 자료의 기본단위를 구한다. */

  small_unit=(u[0]<=u[1] && u[0]<=u[2]) ? u[0] : min(u[1],u[2]);
  /* 3개의 기본단위중 가장 작은 기본단위를 최소단위로 선정. 특별한 자료가
     아닌 이상 3개의 자료를 가지고 최소단위를 구할 수 있다. */

  max=data[n];
  min=data[1];
  range=fabs(max-min); /* 자료의 범위 */
  class_interval=(small_unit==1) ? ceil(range/class) :
		 ceil((range/class/small_unit))*small_unit; /* 계급 간격 */

  class_mlimit=(small_unit==min) ? small_unit/2 :
	       min-(((class*class_interval)-range)/2);

  if (unit(class_mlimit)>=small_unit)
    class_mlimit=(small_unit>=1) ? class_mlimit+(small_unit/2) :
    class_mlimit+(small_unit*5);
  make_table(class);

  init_graph();
  fpoly();
  cpoly();
  fhisto();
  chisto();
  closegraph();

  basic_stat();
}

void prn_data(void)
{
  int i;
  printf("\t\t\t      Data List \n");
  printf("\t\t\t     -----------\n");
  for (i=1;i<=n;i++) {
    printf("%12.6lf   ",*(pdata+i));
    if (i % 5==0) printf("\n");
  }
  getch();
  clrscr();
}

void sort(int left,int right) /* Quick sort */
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
    sort(left,j);
  if (i<right)
    sort(i,right);
}

void get_gcmdata(void)
{
  int i,j=0;
  int temp,same_data=0,sum_same_data=0;

  for (i=1;i<=n+1;i++) {
    temp=*(pdata+i);
    if (temp==*(pdata+i-1)) same_data+=1;
    else {
      class_n[j]=same_data+1;
      if (class_n[j]>=class_n[j-1])
	if (max_class<class_n[j]) {
	  max_class=class_n[j];
	  mode_class=j;
	}
      j++;
      sum_same_data+=same_data;
      same_data=0;
    }
  }
  gcm_data=n-sum_same_data;
}

double unit(double x) /* 자료의 기본단위를 구하는 함수 */
{
  int i;
  double temp;

  temp=x;
  if (temp==floor(x)) /* 주어진 자료가 정수일 경우 */
    for (i=0;;i++) {
      temp/=10;
      if ((int)temp*10!=x) return pow10(i);
    }
  else /* 주어진 자료가 소수일 경우 */
    for (i=1;;i++) {
      temp*=10;
      if (frac(temp)==0) return pow10(-i);
    }
}

void make_table(int class)
{
  int i,j=1;
  int frequency,sum_frequency=0;

  printf("\t\t\t  Sorting Data List \n");
  printf("\t\t\t -------------------\n");
  for (i=1; i<=n; i++) {
    printf("%12.6lf   ",*(pdata+i));
    if (i % 5==0) printf("\n");
  }
  getch();
  clrscr();

  if (gcm_data>class) { /* 자료의 종류가 계급의 수보다 많을 경우 */
    printf("----------------------------------------------------------------"
	   "--------------\n");
    printf("     Class Interval\t     F    Sum of F    F/N      Sum of F/N  "
	   "Class Value\n");
    printf("----------------------------------------------------------------"
	   "--------------\n");
    for (i=0;i<class;i++) { /* 계급의 수만큼 반복 */
      frequency=0;
      do {
	frequency+=1;
	j++;
	if (j>n) break;
      }	while (class_mlimit+(class_interval*(i+1))>data[j]);

      class_n[i+1]=frequency;
      if (frequency>class_n[i])
	if (max_class<frequency) max_class=frequency;

      sum_frequency+=frequency;
      printf("%-11.5f-%11.5f",class_mlimit+(class_interval*(i)),
	     class_mlimit+(class_interval*(i+1)));
      printf("    %3d    %4d       %5.4f      %5.4f     ",frequency,
	     sum_frequency,(double)frequency/n,(double)sum_frequency/n);
      printf("%-11.6f\n",class_mlimit+(class_interval/2)+(class_interval*i));
    }
    printf("----------------------------------------------------------------"
	   "--------------\n");
  }
  else { /* 자료의 종류가 계급의 수보다 적을 경우 */
    printf("----------------------------------------------------------");
    printf("\n  Data         F      Sum of F       F/N      Sum of F/N\n");
    printf("----------------------------------------------------------\n");
    for (i=1;i<=gcm_data;i++) {
      frequency=class_n[i];
      sum_frequency+=frequency;
      printf("%-11.5f  %3d      %4d\t   ",data[sum_frequency],frequency,
	     sum_frequency);
      printf("%5.4f \t%5.4f\n",(double)frequency/n,(double)sum_frequency/n);
    }
    printf("----------------------------------------------------------\n");
  }
  getch();
  clrscr();
}

void basic_stat(void)
{
  double amean,gmean,hmean,qmean,median,mid_point,tirm05_mean,tirm10_mean,
	 wmean,tirm25_mean,sumsquare,mean_dev,variance,stan_dev,
	 stan_error,coe_var,coe_mean,inter_quart,quart_dev,coe_quart,
	 quart_range,prange,uv,uv_dev,mode,pearson,kurtosis,skewness;

  amean=total/n;
  gmean=g_mean();
  hmean=h_mean();
  qmean=q_mean();
  median=get_median();
  mid_point=(max+min)/2;
  mode=get_mode();
  tirm05_mean=t_mean(5);
  tirm10_mean=t_mean(10);
  tirm25_mean=t_mean(25);
  range=max-min;
  sumsquare=sum_square(amean);
  mean_dev=mean_deviation(amean);
  variance=sumsquare/n;
  stan_dev=sqrt(variance);
  stan_error=stan_dev/sqrt(n);
  coe_var=fabs(stan_dev/amean*100);
  coe_mean=fabs(mean_dev/median*100);
  uv=sumsquare/(n-1);
  uv_dev=sqrt(uv);
  get_quart();
  inter_quart=(q1+q3)/2;
  quart_range=q3-q1;
  quart_dev=quart_range/2;
  prange=per_range();
  coe_quart=fabs(quart_dev/median*100);
  pearson=(amean-mode)/stan_dev;
  kurtosis=mx_temp(amean,4)/(n*pow(variance,2));
  wmean=win_mean();
  skewness=mx_temp(amean,3)/(n*pow(stan_dev,3));
  printf("\t\t\t Basic Statistic \n");
  printf("\t\t\t-----------------\n\n");
  printf("\t** DATA SITUATION **\n\n");
  printf("\tSum of Frequency                  : %15d\n",n);
  printf("\tMinimum Data                      : %15.6lf\n",*(pdata+1));
  printf("\tMaximum Data                      : %15.6lf\n",*(pdata+n));
  printf("\tQ1 Data                           : %15.6lf\n",q1);
  printf("\tQ3 Data                           : %15.6lf\n",q3);
  printf("\tSum of Data(Total Value)          : %15.6lf\n",total);
  getch();
  clrscr();
  printf("\t\t ** REPRESENTATIVE VALUE **\n\n");
  printf("\tArithmetic Mean                   : %15.6lf\n",amean);
  if(gmean==PI) printf("\tGeometric  Mean\t                  :   "
			 "No some data is ninus !!!\n");
  else
    printf("\tGeometric  Mean                   : %15.6lf\n",gmean);
  if(hmean==PI) printf("\tHarmonic  Mean \t                  :   "
			 "No some data is 0 !!!\n");
  else
  printf("\tHarmonic   Mean                   : %15.6lf\n",hmean);
  printf("\tQuadratic  Mean                   : %15.6lf\n",qmean);
  printf("\tMedian                            : %15.6lf\n",median);
  printf("\tMid-Point                         : %15.6lf\n",mid_point);
  printf("\tMode                              : %15.6lf\n",mode);
  printf("\tTrimmed Mean(5%c)                  : %15.6lf\n",37,tirm05_mean);
  printf("\tTirmmed Mean(10%c)                 : %15.6lf\n",37,tirm10_mean);
  printf("\tTrimmed Mean(25%c)                 : %15.6lf\n",37,tirm25_mean);
  printf("\tWinsorized Mean                   : %15.6lf\n",wmean);
  printf("\tInterquartile Value               : %15.6lf\n",inter_quart);
  getch();
  clrscr();
  printf("\n\t\t        ** MEASURE OF DISPERSION **\n\n");
  printf("\tRange                             : %15.6lf\n",range);
  printf("\tInterquartile Mid-Range           : %15.6lf\n",quart_range);
  printf("\t90-10 Percentil Mid-Range         : %15.6lf\n",prange);
  printf("\tSum of Square                     : %15.6lf\n",sumsquare);
  printf("\tMean Deviation                    : %15.6lf\n",mean_dev);
  printf("\tVariance                          : %15.6lf\n",variance);
  printf("\tStandard Deviation                : %15.6lf\n",stan_dev);
  printf("\tStandard Error of Mean            : %15.6lf\n",stan_error);
  printf("\tUnbiased Variance                 : %15.6lf\n",uv);
  printf("\tUnbiased Deviation                : %15.6lf\n",uv_dev);
  printf("\tCoefficient of Variance           : %15.6lf \x25\n",coe_var);
  printf("\tCoefficient of Mean Deviation     : %15.6lf \x25\n",coe_mean);
  printf("\tCoefficient of Quartile Deviation : %15.6lf \x25\n",coe_quart);
  printf("\tQuartely Deviation                : %15.6lf\n",quart_dev);
  printf("\tSkewness by the Moment-3          : %15.6lf\n",skewness);
  printf("\tPearson's Skewness                : %15.6lf\n",pearson);
  printf("\tKurtosis                          : %15.6lf\n",kurtosis);
  getch();
}


double g_mean(void)
{
  int i;
  double g_temp=0;

  for (i=1;i<=n;i++) {
    if (*(pdata+i)<0) return PI;
    else {
      if (*(pdata+i)==0) ;
      else
	g_temp+=log(*(pdata+i));
    }
  }
  return (exp(g_temp/n));
}

double h_mean(void)
{
  int i;
  double h_temp=0;
  for (i=1;i<=n;i++) {
    if (*(pdata+i)==0) return PI;
    else
      h_temp+=1./ *(pdata+i);
  }
  return n/h_temp;
}

double q_mean(void)
{
  int i;
  double sum_square=0;

  for (i=1;i<=n;i++)
    sum_square+=pow(*(pdata+i),2);
  return sqrt(sum_square/n);
}

double get_median(void)
{
  if (n & 1) return *(pdata+((n+1)/2));
  else
    return (*(pdata+(n/2))+*(pdata+((n+1)/2)))/2;
}

double t_mean(double per)
{
  int i,tirm_temp;
  double tirm_total=0;

  per/=100;
  tirm_temp=round(n*per);

  for (i=1;i<=tirm_temp;i++)
    tirm_total+=(*(pdata+i)+*(pdata+(n+1-i)));
  return (total-tirm_total)/(n-2*tirm_temp);
}

double win_mean(void)
{
  int i,tirm_temp;
  double tirm_total=0;

  tirm_temp=round(n*.25);

  for (i=1;i<=tirm_temp;i++) {
    tirm_total+=((*(pdata+i)+*(pdata+(n+1-i)))-(q3+q1));
  }
  return (total-tirm_total)/n;
}

double sum_square(double mean)
{
  int i;
  double sumsquare=0;
  for (i=1; i<=n ; i++)
    sumsquare+=pow((*(pdata+i)-mean),2);
  return sumsquare;
}

double mean_deviation(double mean)
{
  int i;
  double deviation=0;
  for (i=1 ; i<=n ;i++)
    deviation+=fabs(*(pdata+i)-mean);
  return deviation/n;
}

void get_quart(void)
{
  q1=data[round(n*.25)];
  q3=data[round(n*.75)];
}

double per_range(void)
{
  int i=round(n/10);
  return data[n-i]-data[i];
}

double get_mode(void)
{
  int i;
  if (gcm_data<=class) return mode_class;
  else
    return ((class_mlimit+(class_interval*(mode_temp-1)))+
	    (class_n[mode_temp+1]/(double)(class_n[mode_temp-1]+
	    class_n[mode_temp+1]))*class_interval);
}

double mx_temp(double mean,int k)
{
  int i;
  double mx_sum=0;
  for (i=1;i<=n;i++)
    mx_sum+=pow(*(pdata+i)-mean,k);
  return mx_sum;
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
  if (gcm_data>class) x_temp=x_range/class;
  else
    x_temp=(int)x_range/gcm_data;
  rectangle(0,0,getmaxx(),getmaxy());
  setcolor(LIGHTRED);
  outtextxy(getmaxx()/2-(textwidth(str)/2),getmaxy()-textheight(str)-2,str);
  setcolor(WHITE);
  line(x0,y0,x0,yn);
  line(x0,y0,xn,y0);
  gprint(x0-20,y0-3,"%1d",0);
}

void title(char *title_name)
{
  char str1[]="Frequency";
  char str2[]="Class Value";
  char str3[]="Rel-Frequency";

  board();
  setcolor(LIGHTCYAN);
  settextstyle(DEFAULT_FONT,HORIZ_DIR,2);
  outtextxy(getmaxx()/2-(textwidth(title_name)/2),
	    textheight(title_name)-10,title_name);
  settextstyle(DEFAULT_FONT,HORIZ_DIR,0);
  setcolor(LIGHTGREEN);
  outtextxy(x0-35,y0-y_range-15,str1);
  outtextxy(x0+x_range-30,y0+30,str2);
  if (*title_name=='C') outtextxy(getmaxx()-110,y0-y_range-15,str3);
  setcolor(YELLOW);
  if (*title_name=='C') {
    gprint(x0-40,y0+10,"N=%-4d",n);
    gprint(x0-40,y0+20,"F/Dot=%-6.1f",(double)n/10);
    setcolor(WHITE);
    line(x0+x_range,y0-y_range,x0+x_range,y0);
  }
  else {
    gprint(x0-40,y0+10,"Max=%d",(int)max_class);
    if (max_class>10)
      gprint(x0-40,y0+20,"F/Dot=%-4.1f",max_class/10);
    else
      gprint(x0-40,y0+20,"F/Dot=1");
    setcolor(WHITE);
  }
}

void fpoly(void)
{
  int i,j,y1,sum=0;

  title("Frequency Polygon");

  if (max_class<10)
    for (i=1;i<=max_class;i++) {
      line(x0,y0-y_range/max_class*i,x0-4,y0-y_range/max_class*i);
      gprint(x0-25,y0-y_range/max_class*i-3,"%2d",i);
    }
  else
    for (i=1;i<=3;i++) {
      line(x0,y0-y_range/3*i,x0-4,y0-y_range/3*i);
      gprint(x0-50,y0-y_range/3*i-2,"%5.1f",max_class/3*i);
    }

  if (gcm_data>class) {
    for (i=1;i<=class;i++) {
      if (i % 2==0) y1=18;
      else y1=8;

      if (class<=12) {
	line(x0+x_temp*i,y0,x0+x_temp*i,y0+4);
	gprint(x0+x_temp*i-60,y0+y1,"%11.4lf",
	       class_mlimit+class_interval/2*(i*2-1));
      }
      else {
	for (j=1;j<=3;j++) {
	  line(x0+x_range*j/3,y0,x0+x_range*j/3,y0+4);
	  gprint(x0+x_range*j/3-60,y0+8,"%11.4lf",
		 class_mlimit+(class_interval*class)/3*j);
	}

	for (j=1;j<=50;j++)
	  putpixel(x0+x_temp*i,y0-y_range*j/50,LIGHTMAGENTA);

	setcolor(YELLOW);
	gprint(x0-40,y0+30,"C1=%-11.4f",class_mlimit+(class_interval/2));
	gprint(x0-40,y0+40,"C/Dot=%-11.4f",class_interval);
	setcolor(WHITE);
	gprint(x0+x_temp-5,y0+8,"C1");
      }

      circle(x0+x_temp*i,y0-(double)(class_n[i]/max_class*y_range),3);
      setfillstyle(SOLID_FILL,LIGHTBLUE);
      floodfill(x0+x_temp*i,
		y0-(double)(class_n[i]/max_class*y_range),WHITE);
    }

    for (i=1;i<class;i++)
      line(x0+x_temp*i,y0-(double)(class_n[i]/max_class*y_range),
	   x0+x_temp*(i+1),
	   y0-(double)(class_n[i+1]/max_class*y_range));

    if (max_class<10)
      for (i=1;i<=max_class;i++)
	for (j=1;j<=50;j++)
	  putpixel(x0+x_range*j/50,y0-y_range*i/max_class,LIGHTRED);
    else
      for (i=1;i<=10;i++)
	for (j=1;j<=50;j++)
	  putpixel(x0+x_range*j/50,y0-y_range*i/10,LIGHTRED);

  }
  else {
    for (i=1;i<=gcm_data;i++) {
      sum+=class_n[i];
      if (i % 2==0) y1=18;
      else y1=8;

      if (gcm_data<=12) {
	line(x0+x_temp*i,y0,x0+x_temp*i,y0+4);
	gprint(x0+x_temp*i-60,y0+y1,"%11.4lf",data[sum]);
      }
      else {
	for (j=1;j<=3;j++) {
	  line(x0+x_range*j/3,y0,x0+x_range*j/3,y0+4);
	  gprint(x0+x_range*j/3-60,y0+8,"%11.4lf",
		 data[1]+(data[n]+data[1])/3*j);
	}

	  for (j=1;j<=50;j++)
	    putpixel(x0+x_temp*i,y0-y_range*j/50,LIGHTMAGENTA);

	setcolor(YELLOW);
	gprint(x0-40,y0+30,"C1=%-11.4f",data[1]);
	gprint(x0-40,y0+40,"C/Dot=%-11.4f",data[n]-data[1]/gcm_data);
	setcolor(WHITE);
	gprint(x0+x_temp-5,y0+8,"C1");
      }

      circle(x0+x_temp*i,y0-(class_n[i]/max_class*y_range),3);
      setfillstyle(SOLID_FILL,LIGHTBLUE);
      floodfill(x0+x_temp*i,y0-(class_n[i]/max_class*y_range),
		WHITE);
    }

    for (i=1;i<gcm_data;i++)
      line(x0+x_temp*i,y0-(double)(class_n[i]/max_class*y_range),
	   x0+x_temp*(i+1),y0-(double)(class_n[i+1]/max_class*y_range));

    if (max_class<10)
      for (i=1;i<=max_class;i++)
	for (j=1;j<=50;j++)
	  putpixel(x0+x_range*j/50,y0-y_range*i/max_class,LIGHTRED);
    else
      for (i=1;i<=10;i++)
	for (j=1;j<=50;j++)
	  putpixel(x0+x_range*j/50,y0-y_range*i/10,LIGHTRED);

  }
  getch();
}

void cpoly(void)
{
  int i,j,y1,sum=0;

  title("Cumul-Frequency polygon");

  for (i=1;i<=10;i++) {
    line(x0,y0-y_range*i/10,x0-5,y0-y_range*i/10);
    line(x0+x_range,y0-y_range*i/10,x0+x_range+4,y0-y_range*i/10);
    gprint(x0-55,y0-(y_range*i/10+3),"%6.1f",n/10.*i);
    gprint(x0+x_range+7,y0-(y_range*i/10+3),"%3.1f",0.1*i);
    for (j=1;j<=50;j++)
      putpixel(x0+x_range*j/50,y0-y_range*i/10,LIGHTRED);
  }

  if (gcm_data>class) {
    for (i=1;i<=class;i++) {
      if (i % 2==0) y1=18;
      else y1=8;

      if (class<=12) {
	line(x0+x_temp*i,y0,x0+x_temp*i,y0+4);
	gprint(x0+x_temp*i-60,y0+y1,"%11.4lf",
	       class_mlimit+(class_interval/2*(i*2-1)));
      }
      else {
	for (j=1;j<=3;j++) {
	  line(x0+x_range*j/3,y0,x0+x_range*j/3,y0+4);
	  gprint(x0+x_range*j/3-60,y0+8,"%11.4lf",
		 class_mlimit+(class_interval*class/3*j));
	}

	for (j=1;j<=50;j++)
	  putpixel(x0+x_temp*i,y0-y_range*j/50,LIGHTMAGENTA);

	setcolor(YELLOW);
	gprint(x0-40,y0+30,"C1=%-11.4f",class_mlimit+(class_interval/2));
	gprint(x0-40,y0+40,"C/Dot=%-11.4f",class_interval);
	setcolor(WHITE);
	gprint(x0+x_temp-5,y0+8,"C1");
      }

      sum+=class_n[i];
      circle(x0+x_temp*i,y0-(double)sum/n*y_range,3);
      setfillstyle(SOLID_FILL,LIGHTBLUE);
      floodfill(x0+x_temp*i,y0-(double)sum/n*y_range,WHITE);
    }
    sum=0;

    for (i=1;i<class;i++) {
      sum+=class_n[i];
      line(x0+x_temp*i,y0-(double)sum/n*y_range,
	   x0+x_temp*(i+1),y0-(double)(sum+class_n[i+1])/n*y_range);
    }
    line(x0,y0,x0+x_temp,y0-(double)class_n[1]/n*y_range);
  }
  else {
    for (i=1;i<=gcm_data;i++) {
      sum+=class_n[i];
      if (i % 2==0) y1=18;
      else y1=8;

      if (gcm_data<=12) {
	line(x0+x_temp*i,y0,x0+x_temp*i,y0+4);
	gprint(x0+x_temp*i-60,y0+y1,"%11.4lf",data[sum]);
      }
      else {
	for (j=1;j<=3;j++) {
	  line(x0+x_range*j/3,y0,x0+x_range*j/3,y0+4);
	  gprint(x0+x_range*j/3-60,y0+8,"%11.4lf",
		 data[1]+(data[n]+data[1])/3*j);
	}

	for (j=1;j<=50;j++)
	  putpixel(x0+x_temp*i,y0-y_range*j/50,LIGHTMAGENTA);

	setcolor(YELLOW);
	gprint(x0-40,y0+30,"C1=%-11.4f",class_mlimit+(class_interval/2));
	gprint(x0-40,y0+40,"C/Dot=%-11.4f",class_interval);
	setcolor(WHITE);
	gprint(x0+x_temp-5,y0+8,"C1");
      }

    }
    sum=0;

    for (i=1;i<=gcm_data;i++) {
      sum+=class_n[i];
      circle(x0+x_temp*i,y0-(double)sum/n*y_range,3);
      setfillstyle(SOLID_FILL,LIGHTBLUE);
      floodfill(x0+x_temp*i,y0-(double)sum/n*y_range,WHITE);
      }
    sum=0;

    for (i=1;i<gcm_data;i++) {
      sum+=class_n[i];
      line(x0+x_temp*i,y0-(double)sum/n*y_range,x0+x_temp*(i+1),
	   y0-(double)(sum+class_n[i+1])/n*y_range);
    }
    line(x0,y0,x0+x_temp,y0-(double)class_n[1]/n*y_range);
  }
  getch();
}

void fhisto(void)
{
  int i,j,y1,temp,sum=0;

  title("Frequency Histogram");

  if (max_class<10)
    for (i=1;i<=max_class;i++) {
      line(x0,y0-y_range/max_class*i,x0-4,y0-y_range/max_class*i);
      gprint(x0-25,y0-y_range/max_class*i-3,"%2d",i);
    }
  else
    for (i=1;i<=3;i++) {
      line(x0,y0-y_range/3*i,x0-4,y0-y_range/3*i);
      gprint(x0-50,y0-y_range/3*i-2,"%5.1f",max_class/3*i);
    }

  if (gcm_data>class) {
    temp=x_range/(class*2+1);

    for (i=1;i<=class;i++) {
      if (i % 2==0) y1=18;
      else y1=8;

      if (class<=12) {
	line(x0+temp*2*i,y0,x0+temp*2*i,y0+4);
	gprint(x0+temp*2*i-60,y0+y1,"%11.4lf",
	       class_mlimit+class_interval/2*(i*2-1));
      }
      else {
	for (j=1;j<=3;j++) {
	  line(x0+x_range*j/3,y0,x0+x_range*j/3,y0+4);
	  gprint(x0+x_range*j/3-60,y0+8,"%11.4lf",
		 class_mlimit+(class_interval*class)/3*j);
	}

	setcolor(YELLOW);
	gprint(x0-40,y0+30,"Cl=%-11.4f",class_mlimit+(class_interval/2));
	gprint(x0-40,y0+40,"C/Pole=%-11.4f",class_interval);
	setcolor(WHITE);
	line(x0+x_temp,y0,x0+x_temp,y0+4);
	gprint(x0+x_temp-5,y0+8,"C1");
      }

      rectangle(x0+temp*(i*2-1),y0-(y_range/max_class*class_n[i]),
	  x0+temp*(i*2+1),y0);
    }

    if (max_class<10)
      for (i=1;i<=max_class;i++)
	for (j=1;j<=50;j++)
	  putpixel(x0+x_range*j/50,y0-y_range*i/max_class,LIGHTRED);
    else
      for (i=1;i<=10;i++)
	for (j=1;j<=50;j++)
	  putpixel(x0+x_range*j/50,y0-y_range*i/10,LIGHTRED);

  }
  else {
    temp=x_range/(gcm_data*2+1);

    for (i=1;i<=gcm_data;i++) {
      sum+=class_n[i];
      if (i % 2==0) y1=18;
      else y1=8;

      if (gcm_data<=12) {
	line(x0+temp*2*i,y0,x0+temp*2*i,y0+4);
	gprint(x0+temp*2*i-60,y0+y1,"%11.4lf",data[sum]);
      }
      else {
	for (j=1;j<=3;j++) {
	 line(x0+x_range/3*j,y0,x0+x_range/3*j,y0+4);
	 gprint(x0+x_range/3*j-30,y0+8,"%15.6lf",
		data[1]+(data[n]+data[1])/3*j);
	}

	setcolor(YELLOW);
	gprint(x0-40,y0+30,"C1=%-11.4f",data[1]);
	gprint(x0-40,y0+40,"C/Pole=%-11.ff",data[n]-data[1]/gcm_data);
	setcolor(WHITE);
	line(x0+x_temp,y0,x0+x_temp,y0+4);
	gprint(x0+x_temp-5,y0+8,"C1");
      }

      rectangle(x0+temp*(i*2-1),y0-(y_range/max_class*class_n[i]),
	  x0+temp*(i*2+1),y0);
    }

    if (max_class<10)
      for (i=1;i<=max_class;i++)
	for (j=1;j<=50;j++)
	  putpixel(x0+x_range*j/50,y0-y_range*i/max_class,LIGHTRED);
    else
      for (i=1;i<=10;i++)
	for (j=1;j<=50;j++)
	  putpixel(x0+x_range*j/50,y0-y_range*i/10,LIGHTRED);

  }
  getch();
}

void chisto(void)
{
  int i,j,y1,temp,sum=0;

  title("Cumul-Frequency Histogram");

  for (i=1;i<=10;i++) {
    line(x0,y0-y_range*i/10,x0-5,y0-y_range*i/10);
    line(x0+x_range,y0-y_range*i/10,x0+x_range+4,y0-y_range*i/10);
    gprint(x0-55,y0-(y_range*i/10+3),"%6.1f",n/10.*i);
    gprint(x0+x_range+7,y0-(y_range*i/10+3),"%3.1f",.1*i);
    for (j=1;j<=50;j++)
     putpixel(x0+x_range*j/50,y0-y_range*i/10 ,LIGHTRED);
  }

  if (gcm_data>class) {
    temp=x_range/(class*2+1);
    for (i=1;i<=class;i++) {
      sum+=class_n[i];
      if (i % 2==0) y1=18;
      else y1=8;

      if (class<=12) {
	line(x0+temp*2*i,y0,x0+temp*2*i,y0+4);
	gprint(x0+temp*2*i-60,y0+y1,"%11.4lf",
	       class_mlimit+class_interval/2*(i*2-1));
      }
      else {
	for (j=1;j<=3;j++) {
	  line(x0+x_range/3*j,y0,x0+x_range/3*j,y0+4);
	  gprint(x0+x_range*j/3-60,y0+8,"%11.4f",
		 data[1]+((data[n]+data[1])/3*j));
	}

	setcolor(YELLOW);
	gprint(x0-40,y0+30,"C1=%-11.4f",class_mlimit+(class_interval/2));
	gprint(x0-40,y0+40,"C/Pole=%-11.f",class_interval);
	setcolor(WHITE);
	line(x0+x_temp,y0,x0+x_temp,y0+4);
	gprint(x0+x_temp-5,y0+8,"C1");
      }

      rectangle(x0+temp*(i*2-1),y0-(double)sum/n*y_range,x0+temp*(i*2+1),y0);
    }
  }
  else {
    temp=x_range/(gcm_data*2+1);
    for (i=1;i<=gcm_data;i++) {
      sum+=class_n[i+1];
      if (i % 2==0) y1=18;
      else y1=8;

      if (gcm_data<=12) {
	line(x0+temp*2*i,y0,x0+temp*2*i,y0+4);
	gprint(x0+temp*2*i-60,y0+y1,"%11.4lf",data[sum]);
      }
      else {
	for (j=1;j<=3;j++) {
	  line(x0+x_range/3*j,y0,x0+x_range/3*j,y0+4);
	  gprint(x0+x_range*j/3-60,y0+8,"%11.4f",
		 data[1]+((data[n]+data[1])/3*j));
	}

	setcolor(YELLOW);
	gprint(x0-40,y0+30,"C1=%-11.4f",class_mlimit+(class_interval/2));
	gprint(x0-40,y0+40,"C/Pole=%-11.f",class_interval);
	setcolor(WHITE);
	line(x0+x_temp,y0,x0+x_temp,y0+4);
	gprint(x0+x_temp-5,y0+8,"C1");
      }

      rectangle(x0+temp*(i*2-1),y0-(double)sum/n*y_range,x0+temp*(i*2+1),y0);
    }
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


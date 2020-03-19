/* Program : FREQUEN.C                */
/* Author  : Ryu choong hyun	      */
/* Date    : 94.2.20.                 */
/* Note    : Frequency Table ver 0.2  */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>

#define printerrmsg(s) (fputs("\n" s "\n",stderr),exit(EXIT_FAILURE))
#define round(x) ((x>0) ? floor(x+.5) : ceil(x-.5)) /* 반올림 함수 정의 */
#define frac(x) (x-floor(x)) /* 소수부를 구하는 매크로 */
#define LOG2 0.301
#define MAX 1000 /* 자료입력의 최대 갯수 */
#define CLASS 35 /* 계급의 최대 갯수 */
#define LOOF_MAX 3 /* 최소단위를 얻기 위한 루프의 반복수 */

double unit(double x);
void sort(int left,int right); /* Quick Sort */
void make_table(int class); /* 돗수분포표 작성 */
void get_gcmdata(void);

int n,gcm_data,class,class_n[CLASS];
double data[MAX],*pdata,class_interval,class_mlimit;

void main(int agrc,char *agrv[])
{
  FILE *stream;
  double num,min,max,small_unit,u[LOOF_MAX],range;
  int i;

  pdata=data; /* 메모리상에서 연속인 자료의 연산은 포인터 연산이 효율적 */
  clrscr();
  if (agrc<=1) {
    puts("Usage : frequen datafile ");
    exit(EXIT_FAILURE);
  }
  stream=fopen(agrv[1],"rt");
  while (fscanf(stream,"%lf\n",&num)!=EOF) {
    ++n; /* 자료의 갯수 */
    *(pdata+n)=num;
    /* 이 프로그램에서는 자료처리의 특성상
    자료의 첫 원소가 data[0]이 아닌 data[1]에 저장됨을 주의해야 한다. */
  }

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
}

void sort(int left,int right) /* Quick sort */
{
  int i=left,j=right;
  double temp,mid=data[(left+right)/2];

  for (;;) {
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
    else break;
  }

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
}

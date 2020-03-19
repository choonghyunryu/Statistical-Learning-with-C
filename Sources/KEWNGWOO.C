/* KEWNGWOO.C */  
#include <stdio.h>
#include <math.h>

double fac(unsigned x);
double permu(unsigned n,unsigned r);
double com(unsigned n,unsigned r);
double repeat_permu(unsigned n,unsigned r);
double repeat_com(unsigned n,unsigned r);

void main(void)
{
  printf("10 factorial               = %10.2f\n",fac(10));
  printf("10 permutation 4           = %10.2f\n",permu(10,4));
  printf("10 repeated permutation 4  = %10.2f\n",repeat_permu(10,4));
  printf("10 combination 4           = %10.2f\n",com(10,4));
  printf("10 repeated combination 4  = %10.2f\n",repeat_com(10,4));
}

double fac(unsigned x)
{
  return ((x==0) ? 1 : x*fac(x-1));
}

double permu(unsigned n,unsigned r)
{
  return ((r==0) ? 1 : (n==r) ? fac(n) : fac(n)/fac(n-r));
}

double com(unsigned n,unsigned r)
{
  return ((n==r || r==0) ? 1 : (n<r) ? 0 : fac(n)/(fac(r)*fac(n-r)));
}

double repeat_permu(unsigned n,unsigned r)
{
  return pow(n,r);
}

double repeat_com(unsigned n,unsigned r)
{
  return com(n+r-1,r);
}


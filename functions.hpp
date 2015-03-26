#ifndef __FUNCTIONS_HPP
#define __FUNCTIONS_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;


//point in the phase space
struct Point {
  double m_energy;
  vector<double> m_p;
};


//Compute the Lagrangian L[x(t),x'(t);t]
double L_func(const vector<double> x, const vector<double> xd){
  double L = 0.;
  int size = static_cast<int>(x.size());
  for(int i = 1; i < size-1; ++i){
    double d = x[0] - x[i];
    double l = 1.e12;
    if(d > 0)
      l = sqrt( (1.+xd[i]*xd[i]) / d );
    L += l;
  }
  return L;
}

//Compute the centered derivative of x
void deriv(const vector<double>& t, const vector<double>& x, vector<double>& xd){
  xd[0] = 0.;
  int size = static_cast<int>(t.size());
  for(int i = 1; i < size-1; ++i)
    xd[i] = (x[i+1]-x[i-1]) / (t[i+1]-t[i-1]);
  xd[size-1] = 0.;
}

//Inner product (integration with respect to the meassure)
void mult(const vector<vector<double> >& x, const double f[2], const double mu[2], Point& p){
  int size = static_cast<int>(x.size());
  int sz = static_cast<int>(x[0].size());
  p.m_energy = 0.;
  for(int i = 0; i < size; ++i){
    for(int j = 0; j < sz; ++j)
      p.m_p[j] = 0.;
    for(int j = 0; j < sz; ++j)
      p.m_p[j] += x[i][j] * mu[i];
    p.m_energy += f[i] * mu[i];
  }
}

//Interpolant (Lagrange 1st order). n is the number of points to plot.
void xi(const double p0[2], const double p1[2], int n, vector<double>& inter){
  double xmax = max(p0[0],p1[0]);
  double xmin = min(p0[0],p1[0]);
  double h = xmax - xmin;
  double d = h / n;
  for(int i = 0; i < n; ++i){
    double ii = xmin + i * d;
    inter[i] = (p0[1] * (xmax - ii) + p1[1] * (ii - xmin)) / h;
  }
}


#endif

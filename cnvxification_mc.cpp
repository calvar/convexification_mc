
#include <fstream>
#include <iomanip> 
#include <gsl/gsl_rng.h>
#include <cstring>
#include <sys/time.h>

#include "functions.hpp"


int main(){
  //No. of steps
  int N_steps = 500000;
  //Equilibration steps
  int eq_steps = 499000;

  //print frequency
  int Pfr = 10000;

  //read previous?
  bool read = true;

  //define "temerature" to allow better sampling
  double T = 1.e-12;

  //boundary conditions
  double b_left = 1.;
  double b_right = 0.;

  //Define the mesh. The points on the mesh are the ones that are randomly moved during the MC process
  double t_min = 0.;
  double t_max = 10.;
  int N_mesh = 11;
  double Dm = (t_max-t_min) / (N_mesh-1);
  vector<double> mesh(N_mesh, 0.);
  for(int i = 0; i < N_mesh; ++i)
    mesh[i] = i * Dm;
  // //
  // for(int i = 0; i < N_mesh; ++i)
  //   cout << mesh[i] << " ";
  // cout << "\n";
  // //

  //Number of points in which the interpolant is evaluated
  int n = 2;

  //Define range of x
  double l_lim = -10.;
  double u_lim = max(b_left, b_right);

  //initialize random num. gen.
  gsl_rng* ran = gsl_rng_alloc(gsl_rng_mt19937); //mersenne twister rng
  long seed = 1;//time(NULL);
  gsl_rng_set(ran, seed); //gsl_rng_uniform(ran);

  //Initialize evaluation functions. Two of them, to scan through the convex envelope of the functional using a linear convex combination.
  //Initialize the prob. meassure for the two points (functions).
  vector< vector<double> > xm(2, vector<double>(N_mesh, 0.));
  double mu[2] = {1., 0.};

  if(read){
    string line;
    ifstream In("final_f.dat");
    if(! In){
      cout << "Couldn't open " << "final_f.dat" << endl;
      In.close();
      return 1;
    }
    for(int i = 0; i < 2; ++i){
      for(int j = 0; j < N_mesh; ++j)
	In >> xm[i][j]; 
    }
    In >> mu[0] >> mu[1];
    In.close();
  } else {
    for(int i = 0; i < 2; ++i){
      for(int j = 0; j < N_mesh; ++j)
	xm[i][j] = (u_lim-l_lim) * gsl_rng_uniform(ran) + l_lim; 
    }
    //Boundary conditions for x: x[0] = 1., x[Nt-1] = 0.
    for(int i = 0; i < 2; ++i){
      xm[i][0] = b_left;
      xm[i][N_mesh-1] = b_right;
    }
    mu[0] = gsl_rng_uniform(ran);
    mu[1] = 1. - mu[0];
  }
  // //
  // for(int i = 0; i < 2; ++i){
  //   for(int j = 0; j < N_mesh; ++j)
  //     cout << xm[i][j] << " ";
  //   cout << endl;
  // }
  // cout << mu[0] << " " << mu[1] << endl;
  // //

  //Define the complete grid
  int Nt = n * (N_mesh-1) + 1;
  vector<double> t(Nt,0.);
  vector<vector<double> > x(2,vector<double>(Nt,0.));
  for(int i = 0; i < Nt; ++i)
    t[i] = i*Dm / n;
  // //
  // for(int i = 0; i < Nt; ++i)
  //   cout << t[i] << " ";
  // cout << endl;
  // //

  //Evaluate the interpolants
  for(int i = 0; i < 2; ++i){
    for(int j = 0; j < N_mesh-1; ++j){
      double p0[2] = {mesh[j], xm[i][j]};
      double p1[2] = {mesh[j+1], xm[i][j+1]};
      vector<double> inter(n, 0.);
      xi(p0, p1, n, inter);
      x[i][j*n] = xm[i][j];
      for(int k = 1; k < n; ++k)
	x[i][j*n+k] = inter[k];
    } 
    x[i][Nt-1] = xm[i][N_mesh-1]; 
  }
  // //
  // for(int i = 0; i < 2; ++i){
  //   for(int j = 0; j < Nt; ++j)
  //     cout << x[i][j] << " ";
  //   cout << endl;
  // }
  // //

  //Compute the derivative of x
  vector<vector<double> > xd(2, vector<double>(Nt, 0.));
  deriv(t, x[0], xd[0]);
  deriv(t, x[1], xd[1]);
  // //
  // for(int i = 0; i < 2; ++i){
  //   for(int j = 0; j < Nt; ++j)
  //     cout << xd[i][j] << " ";
  //   cout << endl;
  // }
  // //

  //Maximum x move range
  double Dx = (u_lim-l_lim) / 50000;

  //Maximum mu move range
  double dmu = 0.1;

  //Average prob. meassure
  //
  //
  //

  //initialize time
  timeval tim;
  gettimeofday(&tim, NULL);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);

  //Trayectory of the Monte Carlo search over the convex envelope
  vector<Point> tr;
  //To compute acceptance ratio of prob. meas. moves
  int acc_p = 0;
  int tot_p = 0;
  //To compute acceptance ratio of point moves
  int acc_x = 0;
  int tot_x = 0;

  //minimum
  double f[2] = {L_func(x[0], xd[0]), L_func(x[1], xd[1])};
  Point p, mn;
  p.m_energy = 0.;
  p.m_p.resize(Nt, 0.);
  mult(x, f, mu, p);
  mn.m_energy = 1.e6;
  mn.m_p.resize(Nt, 0.);
  if(read){
    mn.m_energy = p.m_energy;
    for(int i = 0; i < Nt; ++i)
      mn.m_p[i] = p.m_p[i]; 
  }
  //
  cout << mn.m_energy << " (";
  for(int i = 0; i < Nt-1; ++i)
    cout << mn.m_p[i] << ",";
  cout << mn.m_p[Nt-1] << ")\n";
  //
  tr.push_back(p);
    
  // Monte Carlo
  Point p_old;
  p_old.m_p.resize(Nt, 0.);
  vector<double> x_old(Nt, 0.);
  vector<double> xd_new(Nt, 0.);
  vector<double> inter(n, 0.);
  
  for(int ns = 0; ns < N_steps; ++ns){
    
    p_old.m_energy = p.m_energy;
    for(int i = 0; i < Nt; ++i)
      p_old.m_p[i] = p.m_p[i];
    
    //Modify randomly the probability meassure
    double rnd_rng = 0.;
    bool test = true;
    while(test){
      rnd_rng = 2*dmu * gsl_rng_uniform(ran) - dmu; 
      double a = mu[0] + rnd_rng;
      double b = mu[1] - rnd_rng;
      //If within range
      if((a <= 1. && a >= 0.) && (b <= 1. && b >= 0.)){
    	//sum to the first point and substract to the second
    	mu[0] = a;
    	mu[1] = b;
    	test = false;
      }
    }
    
    //Point in the convex envelope
    
    mult(x, f, mu, p);
    //cout<<mu[0]<<","<<mu[1]<<" "<<f[0]<<","<<f[1]<<" "<<p.m_energy<<"\n";
    //Accept the move if the value has decreased or temp allows it
    if((p.m_energy <= mn.m_energy) /*|| gsl_rng_uniform(ran) < exp((mn.m_energy-p.m_energy)/T)*/){
      if(ns % Pfr == 0)
    	tr.push_back(p);
      if(p.m_energy <= mn.m_energy){
    	mn.m_energy = p.m_energy;
    	for(int i = 0; i < Nt; ++i)
    	  mn.m_p[i] = p.m_p[i]; 
      }
      acc_p++;
    }else{
      mu[0] -= rnd_rng;
      mu[1] += rnd_rng;
      p.m_energy = p_old.m_energy;
      for(int i = 0; i < Nt; ++i)
    	p.m_p[i] = p_old.m_p[i];
    }
    tot_p++;
    
    //Make a random move of one of the eval. functions on the mesh
    int rnd = (gsl_rng_uniform(ran) < 0.5)? 0 : 1;
    for(int tt = 1; tt < N_mesh-1; ++tt){ //do not move the boundary
      p_old.m_energy = p.m_energy;
      for(int i = 0; i < Nt; ++i)
    	p_old.m_p[i] = p.m_p[i];
      double xm_old = xm[rnd][tt];
      //displace just one point on the mesh
      double disp = 2*Dx * gsl_rng_uniform(ran) - Dx;
      xm[rnd][tt] += disp;
      
      //If within range
      if((xm[rnd][tt]<u_lim) && (xm[rnd][tt]>l_lim)){
    	for(int i = 0; i < Nt; ++i)
    	  x_old[i] = x[rnd][i];
    	//Evaluate the interpolants
    	double p0[2] = {mesh[tt-1], xm[rnd][tt-1]};
    	double p1[2] = {mesh[tt], xm[rnd][tt]};
    	xi(p0, p1, n, inter);
    	for(int j = 1; j < n; ++j)
    	  x[rnd][(tt-1)*n+j] = inter[j];
	
    	p0[0] = mesh[tt]; p0[1] = xm[rnd][tt];
    	p1[0] = mesh[tt+1]; p1[1] = xm[rnd][tt+1];
    	xi(p0, p1, n, inter);
    	x[rnd][tt*n] = xm[rnd][tt];
    	for(int j = 1; j < n; ++j)
    	  x[rnd][tt*n+j] = inter[j];
    	//Compute the derivative of x
    	deriv(t, x[rnd], xd_new);

    	double f_old = f[rnd];
    	f[rnd] = L_func(x[rnd], xd_new);
    	//Point in the convex envelope
    	mult(x, f, mu, p);
    	//Accept the move if the value has decreased or if the temp. allows it
    	if((p.m_energy <= mn.m_energy) /*|| gsl_rng_uniform(ran) < exp((mn.m_energy-p.m_energy)/T)*/){
    	  for(int i = 0; i < Nt; ++i)
    	    xd[rnd][i] = xd_new[i];
    	  if(ns % Pfr == 0)
    	    tr.push_back(p);
    	  if(p.m_energy <= mn.m_energy){
    	    mn.m_energy = p.m_energy;
    	    for(int i = 0; i < Nt; ++i)
    	      mn.m_p[i] = p.m_p[i]; 
    	  }
    	  acc_x++;
    	}else{
    	  xm[rnd][tt] = xm_old;
    	  for(int i = 0; i < Nt; ++i)
    	    x[rnd][i] = x_old[i];
    	  f[rnd] = f_old;
    	  p.m_energy = p_old.m_energy;
    	  for(int i = 0; i < Nt; ++i)
    	    p.m_p[i] = p_old.m_p[i];
    	}
      }else{
    	xm[rnd][tt] = xm_old;
    	p.m_energy = p_old.m_energy;
    	for(int i = 0; i < Nt; ++i)
    	  p.m_p[i] = p_old.m_p[i];
      }
      tot_x++;

    }

  }


  //final time
  gettimeofday(&tim, NULL);
  double t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  double t_elap = (t2 - t1);
  cout << "Elapsed time: " << t_elap << " sec.\n";

  //acceptance ratio
  double rat_p = 0.;
  if(tot_p > 0)
    rat_p = static_cast<float>(acc_p) / tot_p;
  cout << "measure moves acc: " << rat_p << endl;
  double rat_x = 0.;
  if(tot_x > 0)
    rat_x = static_cast<float>(acc_x) / tot_x;
  cout << "point moves acc: " << rat_x << endl;
  //minimum value
  cout << "minimum: " << mn.m_energy << " (";
  for(int i = 0; i < Nt-1; ++i)
    cout << mn.m_p[i] << ",";
  cout << mn.m_p[Nt-1] << ")\n";

  //Save final function
  ofstream Out("final_f.dat");
  Out << setiosflags(ios::fixed) << setprecision(12);
  for(int i = 0; i < 2; ++i){
    for(int j = 0; j < N_mesh; ++j)
      Out << xm[i][j] << " ";
    Out << "\n";
  }
  Out << mu[0] << " " << mu[1] << "\n";
  Out.close();

  //Save function for ploting
  ofstream Plt("function.dat");
  for(int i = 0; i < Nt; ++i)
    Plt << t[i] << " " << p.m_p[i] << "\n";
  Plt.close();

  return 0;
}

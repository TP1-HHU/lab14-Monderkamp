#include "EM.hxx"
#include <cmath>
#include <iostream>
//----------------------------------
using namespace std;
//----------------------------------
class particle{
public:
	particle(){x=0; y=0; z=0;  // location at t = 0
		         px=0; py=0; pz=0;  // momentum at t = -dt
						 gam=0;  // gamma at t= -dt};
	// Print particle position, momentum u and gamma-factor
	void print_info(const double t) const;
	
  double get_x() const{return x;}
private:
	double x,y,z;
	double px,py,pz;
  double gam;
	// Calculate gamma-factor from current momentum

};
//----------------------------------
//----------------------------------
int main(){
  const double FWHM = 10;
  const double a0 = 1.44;
  const double T0 = 200;
  const double Tend = 800;
  const double dt = 0.001;
  const int N = int(Tend/dt + 0.5);

  double Ex,Ey,Ez, Bx, By, Bz;
  double t=0;

  EM em(a0, FWHM,T0);
  particle p;

  for(int i=0; i<N; i++){
	  em.fields(Ex,Ey,Ez,Bx,By,Bz, t + dt*0.5, p.get_x());
	  p.boris_push(Ex, Ey, 0, Bx, 0, Bz, dt);

	  if (i%5 == 0 ) p.print_info(t);
	  t += dt;

  }


  return 0;
}
//----------------------------------

//----------------------------------

//--------------------------------------
void particle::print_info(const double t) const
{
	cout << t << "\t" << x << "\t" << y << "\t" << z
			 << "\t" << px << "\t" << py << "\t" << pz << "\t" << gam << endl;
}

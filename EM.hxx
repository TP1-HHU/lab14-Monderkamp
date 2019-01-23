#ifndef EM_HXX_INCLUDED
#define EM_HXX_INCLUDED
//--------------------------------------------
#include <cmath>
//--------------------------------------------
class EM{
private:
    const double a0;
    const double FWHM;
    const double T0;
    double sigma;
public:
	EM(const double a0, const double FWHM, const double T0=0);
	void fields(double& Ex, double& Ey, double& Ez,
	            double& Bx, double& By, double& Bz,
	            const double t, const double x) const;
};
//--------------------------------------------
EM::EM(const double a0, const double FWHM, const double T0):
       a0(a0), FWHM(FWHM), T0(T0)
    {
        sigma = 0.5*FWHM*2*3.141/sqrt(log(2));
    }
//--------------------------------------------
void EM::fields(double& Ex, double& Ey, double& Ez,
		            double& Bx, double& By, double& Bz,
		            const double t, const double x) const
    {
        const double tau = t - T0;
        const double s2 = sigma*sigma;
        const double phi = x - tau;
        const double a = a0 * exp( -(x-tau)*(x-tau) / s2);

        Ex = 0;
        Ey = a * sin(phi);
        Ez = 0;
        Bx = 0;
        By = 0;
        Bz = a * sin(phi);
    }
//--------------------------------------------

#endif

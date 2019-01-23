#include "EM.hxx"
#include <cmath>
#include <iostream>
#include <vector>
//----------------------------------
using namespace std;
typedef vector<double> vouble;
//----------------------------------
class particle
    {
        public:
        particle()
            {
                x=0; 
                y=0; 
                z=0;    // location at t = 0
                px=0; 
                py=0; 
                pz=0;   // momentum at t = -dt
                gam=0;  // gamma at t= -dt};
                        // Print particle position, momentum u and gamma-factor
            }
        void print_info(const double t) const;
        double get_x() const{return x;}
        void boris_push(const double Ex, const double Ey, const double Ez, const double Bx, const double By, const double Bz, const double dt);
            
        private:
        double x,y,z;
        double px,py,pz;
        double gam;
        double gamma_factor();// Calculate gamma-factor from current momentum
    };
//----------------------------------
vouble crpr(vouble A, vouble B)
    {
    if (A.size() != 3 || B.size() != 3) 
        {
            exit;
        }
    vouble C(3);
    
    C[0] = A[1]*B[2] - A[2]*B[1]; 
    C[1] = A[2]*B[0] - A[0]*B[2];
    C[2] = A[0]*B[1] - A[1]*B[0];
    
    return C;
        
    }
//----------------------------------
int main()
    {
        
        /*
        vouble A(3);
        vouble B(3);
        
        A[0] = 1.0;
        A[1] = 0.0;
        A[2] = 0.0;
        
        B[0] = 0.0;
        B[1] = 1.0;
        B[2] = 0.0;
        
        vouble C(3);
        C = crpr(A,B);
        
        cout << C[0] << "   " << C[1] << "   " << C[2] << endl;
        */
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

        for(int i=0; i<N; i++)
            {
                em.fields(Ex,Ey,Ez,Bx,By,Bz, t + dt*0.5, p.get_x());
                p.boris_push(Ex, Ey, 0, Bx, 0, Bz, dt);

                if (i%5 == 0 ) p.print_info(t);
                t += dt;

            }


        return 0;
    }
//----------------------------------
void particle::boris_push(const double Ex, const double Ey, const double Ez, const double Bx, const double By, const double Bz, const double dt)
    {
        //cout << "cyka blyat it's Boris!" << endl;
        vouble E(3);
        vouble B(3);
        
        
        E[0] = Ex;
        E[1] = Ey;
        E[2] = Ez;
        
        B[0] = Bx;
        B[1] = By;
        B[2] = Bz;
        
        
        vouble t(3);
        
        double gam = gamma_factor();
        
        t[0] = 0.5*dt*B[0]/gam;
        t[1] = 0.5*dt*B[1]/gam;
        t[2] = 0.5*dt*B[2]/gam;
        
        vouble p_n_minus_half(3);
        
        p_n_minus_half[0] = px;
        p_n_minus_half[1] = py;
        p_n_minus_half[2] = pz;
        
        
        vouble p_minus(3);
        
        p_minus[0] = p_n_minus_half[0] + 0.5*dt*E[0];
        p_minus[1] = p_n_minus_half[1] + 0.5*dt*E[1];
        p_minus[2] = p_n_minus_half[2] + 0.5*dt*E[2];
        
        vouble p_strich(3);
        
        p_strich[0] = p_minus[0] + p_minus[1]*t[2] - p_minus[2]*t[1];
        p_strich[1] = p_minus[1] + p_minus[2]*t[0] - p_minus[0]*t[2];
        p_strich[2] = p_minus[2] + p_minus[0]*t[1] - p_minus[1]*t[0];
        
        vouble s(3);
        
        s[0] = 2.0*t[0]/(1.0+t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
        s[1] = 2.0*t[1]/(1.0+t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
        s[2] = 2.0*t[2]/(1.0+t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
        
        vouble p_plus(3);
        
        p_plus[0] = p_minus[0] + p_strich[1]*s[2] - p_strich[2]*s[1];
        p_plus[1] = p_minus[1] + p_strich[2]*s[0] - p_strich[0]*s[2];
        p_plus[2] = p_minus[2] + p_strich[0]*s[1] - p_strich[1]*s[0];      
        
        vouble p_n_plus_half(3);
        
        p_n_plus_half[0] = p_plus[0] + 0.5*dt*E[0];
        p_n_plus_half[1] = p_plus[1] + 0.5*dt*E[1];
        p_n_plus_half[2] = p_plus[2] + 0.5*dt*E[2];
        
        double gamma_n_plus_half = sqrt( 1.0 +  p_n_plus_half[0]*p_n_plus_half[0] 
                                             +  p_n_plus_half[1]*p_n_plus_half[1]
                                             +  p_n_plus_half[2]*p_n_plus_half[2]);
        
        x += dt* p_n_plus_half[0]/gamma_n_plus_half;
        y += dt* p_n_plus_half[1]/gamma_n_plus_half;
        z += dt* p_n_plus_half[2]/gamma_n_plus_half;
        
        px = p_n_plus_half[0];
        py = p_n_plus_half[1];
        pz = p_n_plus_half[2];
        
        
        
    }
//----------------------------------
double particle::gamma_factor()
    {
        return sqrt(1.0 + px*px+py*py+pz*pz);
    }

//--------------------------------------
void particle::print_info(const double t) const
{
	cout << t << "\t" << x << "\t" << y << "\t" << z
			 << "\t" << px << "\t" << py << "\t" << pz << "\t" << gam << endl;
}

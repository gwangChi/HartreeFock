/*
DEFINES THE OPEN SHELL BASIS BASE CLASS AND DERIVED CLASSES:
1. STO_3G
2. 4_31G
3. 6_31G*
4. 6_31G**
*/

#ifndef BASIS_H
#define BASIS_H

#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include "Eigen/Dense"
#include "NUCLEUS.h"

using namespace std;
using namespace Eigen;

int factorial(int n){return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}

//double integral_overlap(BASIS&, BASIS&);
//double integral_kinetic(BASIS&, BASIS&);
//double integral_coulombic(BASIS&, BASIS&, NUCLEI&);
//double integral_two(BASIS&, BASIS&, BASIS&, BASIS&);
class BASIS
{
public:
	//BASIS() = default;
	BASIS(NUCLEI& n, int input_label) : m_nuclei(n), m_label(input_label){}
	//virtual ~BASIS() = default;

	vector<double>& getExponent(){return m_exponents;}
	vector<double>& getCoefficient(){return m_coefficients;}

	NUCLEI& getNuclei(){return m_nuclei;}
	int getSize(){return m_size;}
	int getLabel(){return m_label;}
	/* 
	INTEGRALS TO BE OVERLOADED ON TEMPLATE BASED ON DERIVED CLASSES
	*/
	//static double integral_overlap(BASIS&, BASIS&);
	//static double integral_kinetic(BASIS&, BASIS&);
	//static double integral_coulombic(BASIS&, BASIS&, NUCLEI&);
	//static double integral_two(BASIS&, BASIS&, BASIS&, BASIS&);

protected:
	int m_size;		// SIZE OF THE GAUSSIAN CONTRACTION
	int m_label;	// THE LABEL USED BY SPECIFIC BASIS TO SPECIFY TYPE OF INTEGRATION
	
	vector<double> m_exponents;
	vector<double> m_coefficients;

	NUCLEI& m_nuclei;

private:
	static const int SET_SIZE;	// SIZE OF THE TOTAL BASIS SET
};


/*
	STO_3G LABEL FOR ATMOIC ORBITALS:
	1 - 1S
	2 - 2S
	3 - 2Px
	4 - 2Py
	5 - 2Pz
	6 - 3S
*/
class STO_3G : public BASIS
{
public:
	STO_3G(int input_label, double input_zeta, NUCLEI& n):BASIS(n, input_label), m_zeta(input_zeta)
	{
		m_size = 3;

		switch(m_label)
		{
			case 1:
				m_coefficients.push_back(0.444635);
				m_coefficients.push_back(0.535328);
				m_coefficients.push_back(0.154329);

				m_exponents.push_back(0.109818 * pow(m_zeta, 2));
				m_exponents.push_back(0.405771 * pow(m_zeta, 2));
				m_exponents.push_back(2.22766 * pow(m_zeta, 2));

				break;

			default:
				cout<<"Not yet availible.\n";
		}
	}

	double getSlaterExponent(){return m_zeta;}

private:
	double m_zeta;
};


/*
INTEGRAL DEFINITIONS
*/
double integral_overlap(BASIS& a, BASIS& b)
{
	double ans = 0;

	vector<double> &Ca = a.getCoefficient(),
		&Cb = b.getCoefficient(),
		&Ea = a.getExponent(),
		&Eb = b.getExponent();

	
	Vector3d vdist = a.getNuclei().getPosition() - b.getNuclei().getPosition();
	double sqdist = pow(vdist.norm(), 2);

	switch(a.getLabel())
	{
		case 1:
		
		switch(b.getLabel())
		{
			case 1:

			for(int i = 0; i < a.getSize(); i++)
				for(int j = 0; j < b.getSize(); j++)
					ans += Ca[i] * Cb[j] * sqrt(sqrt(pow(4*Ea[i]*Eb[j]/pow(M_PI, 2), 3))) * sqrt(pow(M_PI/(Ea[i]+Eb[j]), 3)) * exp(-1 * Ea[i] * Eb[j]/(Ea[i]+Eb[j]) * sqdist);

			cout<<endl<<ans<<endl;
			break;
		}
		break;
	}

	return ans;
}

double integral_kinetic(BASIS& a, BASIS& b)
{
	double ans = 0;

	vector<double> &Ca = a.getCoefficient(),
		&Cb = b.getCoefficient(),
		&Ea = a.getExponent(),
		&Eb = b.getExponent();

	
	Vector3d vdist = a.getNuclei().getPosition() - b.getNuclei().getPosition();
	double sqdist = pow(vdist.norm(), 2);

	switch(a.getLabel())
	{
		case 1:
		
		switch(b.getLabel())
		{
			case 1:
			for(int i = 0; i < a.getSize(); i++)
				for(int j = 0; j < b.getSize(); j++)
					ans += Ca[i] * Cb[j] * sqrt(sqrt(pow(4*Ea[i]*Eb[j]/pow(M_PI,2), 3))) *
							Ea[i]*Eb[j]/(Ea[i]+Eb[j]) * (3-2*Ea[i]*Eb[j]/(Ea[i]+Eb[j])*sqdist) *
							sqrt(pow(M_PI/(Ea[i]+Eb[j]), 3)) * exp(-1*Ea[i]*Eb[j]/(Ea[i]+Eb[j])*sqdist);
			break;
		}
		break;
	}

	return ans;
}

double integral_coulombic(BASIS& a, BASIS& b, NUCLEI& n)
{
	double ans = 0;

	vector<double> &Ca = a.getCoefficient(),
		&Cb = b.getCoefficient(),
		&Ea = a.getExponent(),
		&Eb = b.getExponent();

	
	Vector3d vdist_ab = a.getNuclei().getPosition() - b.getNuclei().getPosition();
	double sqdist_ab = pow(vdist_ab.norm(), 2);


	switch(a.getLabel())
	{
		case 1:		// 1S ORBITAL	
		switch(b.getLabel())
		{
			case 1:	// 1S ORBITAL
			for(int i = 0; i < a.getSize(); i++)
				for(int j = 0; j < b.getSize(); j++)
				{
					Vector3d Rp = (Ea[i]*a.getNuclei().getPosition()+Eb[j]*
						b.getNuclei().getPosition())/(Ea[i]+Eb[j]);
					Vector3d vdist_pn = Rp - n.getPosition();
					double sqdist_pn = pow(vdist_pn.norm(), 2);

					if(&a.getNuclei() == &n && &b.getNuclei() == &n)
					{
						ans += Ca[i] * Cb[j] * pow(4*Ea[i]*Eb[j]/pow(M_PI, 2), 3/4) * 
								(-1)*2*M_PI/(Ca[i]+Cb[j])*n.getCharge()*
								exp(-1*Ea[i]*Eb[j]/(Ea[i]+Eb[j])*sqdist_ab);
						continue;
					}

					ans += Ca[i] * Cb[j] * pow(4*Ea[i]*Eb[j]/pow(M_PI, 2), 3/4) *
								(-1)*2*M_PI/(Ca[i]+Cb[j])*n.getCharge()*
								exp(-1*Ea[i]*Eb[j]/(Ea[i]+Eb[j])*sqdist_ab)*
								1/2*sqrt(M_PI/(Ea[i]+Eb[j])/sqdist_pn)*
								erf(sqrt((Ea[i]+Eb[j])*sqdist_pn));
				}
				break;
		}
		break;
	}

	return ans;
}

double integral_two(BASIS& a, BASIS& b, BASIS& c, BASIS& d)
{
	double ans = 0;

	Vector3d vdist_ab = a.getNuclei().getPosition() - b.getNuclei().getPosition();
	double sqdist_ab = pow(vdist_ab.norm(), 2);

	Vector3d vdist_cd = c.getNuclei().getPosition() - d.getNuclei().getPosition();
	double sqdist_cd = pow(vdist_cd.norm(), 2);

	switch(a.getLabel())
	{
		case 1:
		
		switch(b.getLabel())
		{
			case 1:

			for(int i = 0; i < a.getSize(); i++)
				for(int j = 0; j < b.getSize(); j++)
					for(int k = 0; k < c.getSize(); k++)
						for(int l = 0; l < d.getSize(); l++)
						{
							double A = a.getExponent()[i], B = b.getExponent()[j], C = c.getExponent()[k], D = d.getExponent()[l];

							Vector3d Rp = (A*a.getNuclei().getPosition()+B*b.getNuclei().getPosition())/(A+B);
							Vector3d Rq = (C*c.getNuclei().getPosition()+D*d.getNuclei().getPosition())/(C+D);

							Vector3d Rpq = Rp-Rq;
							double sqdist_pq = pow(Rpq.norm(), 2);

							//FIRST CONSIDER POLE OF THE F0 FUNCTION AR ZERO
							if(sqdist_pq==0)
							{
								ans += a.getCoefficient()[i] * b.getCoefficient()[j] * c.getCoefficient()[k] * d.getCoefficient()[l]
									* sqrt(sqrt(pow(16*A*B*C*D/pow(M_PI, 4), 3)))
									* 2*sqrt(pow(M_PI, 5))/((A+B)*(C+D)*sqrt(A+B+C+D))
									* exp(-A*B/(A+B)*sqdist_ab - C*D/(C+D)*sqdist_cd);
								continue;
							}

							ans += a.getCoefficient()[i] * b.getCoefficient()[j] * c.getCoefficient()[k] * d.getCoefficient()[l]
									* sqrt(sqrt(pow(16*A*B*C*D/pow(M_PI, 4), 3)))
									* 2*sqrt(pow(M_PI, 5))/((A+B)*(C+D)*sqrt(A+B+C+D))
									* exp(-A*B/(A+B)*sqdist_ab - C*D/(C+D)*sqdist_cd)
									* 1/2*sqrt(M_PI/((A+B)*(C+D)/(A+B+C+D)*sqdist_pq))
									* erf(sqrt((A+B)*(C+D)/(A+B+C+D)*sqdist_pq));
						}
			break;
		}
		break;
	}

	return ans;
}


#endif


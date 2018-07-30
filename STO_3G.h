/*
STO_3G CLASS
PROVIDES ALL ONE-ELECTRON INTEGRAL METHODS
*/

#ifndef STO_3G_H
#define STO_3G_H

#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <complex>
#include "Eigen/Dense"
#include "Nuclei.h"

using namespace std;
using namespace Eigen;

int factorial(int n){return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}

// CURRENTLY ALL 1S ORBITALS
class STO_3G
{
public:
	STO_3G(double input_zeta, Nuclei* input_center, int input_spin);

	double getSlaterExponent();
	double* getExponent();
	double* getCoefficient();
	Nuclei* getCenter();

	// ONE ELECTRON INTEGRALS
	double integral_overlap(STO_3G* b);
	double integral_kinetic(STO_3G* b);
	double integral_coulombic(Nuclei* Rc, STO_3G* b);
private:
	// SLATER ORBITAL EXPONENT
	double zeta;
	// GAUSSIAN COEFFICIENTS
	double d[3];
	// GAUSSIAN EXPONENTS
	double a[3];

	// NUMBER OF SLATER ORBITALS USED TO APPROXIMATE
	int size_basis, m_spin;

	Nuclei* m_center;
};

STO_3G::STO_3G(double input_zeta, Nuclei* input_center, int input_spin=0)
{
	zeta = input_zeta;

	d[0] = 0.444635;
	d[1] = 0.535328;
	d[2] = 0.154329;

	a[0] = 0.109818 * zeta * zeta;
	a[1] = 0.405771 * zeta * zeta;
	a[2] = 2.22766 * zeta * zeta;

	size_basis = 1;	// NUMBER OF SLATER ORBITALS TO APPROXIMATE

	m_spin = input_spin;
	m_center = input_center;
}

double STO_3G::getSlaterExponent(){return zeta;}

double* STO_3G::getExponent(){return a;}

double* STO_3G::getCoefficient(){return d;}

Nuclei* STO_3G::getCenter(){return m_center;}

double STO_3G::integral_overlap(STO_3G* b)
{
	double ans = 0;
	Vector3d vdist = m_center->getPosition() - b->getCenter()->getPosition();
	double sqdist = pow(vdist.norm(), 2);

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			double A = a[i], B = b->getExponent()[j];
			ans += d[i] * b->getCoefficient()[j] * sqrt(sqrt(pow(4*A*B/pow(M_PI, 2), 3))) * sqrt(pow(M_PI/(A+B), 3)) * exp(-A*B/(A+B)*sqdist);
		}

	return ans;
}

double STO_3G::integral_kinetic(STO_3G* b)
{
	double ans = 0;
	Vector3d vdist = m_center->getPosition() - b->getCenter()->getPosition();
	double sqdist = pow(vdist.norm(), 2);

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			double A = a[i], B = b->getExponent()[j];
			ans += d[i] * b->getCoefficient()[j] * sqrt(sqrt(pow(4*A*B/pow(M_PI, 2), 3))) * A*B/(A+B)*(3-2*A*B/(A+B)*sqdist) * sqrt(pow(M_PI/(A+B), 3)) * exp(-A*B/(A+B)*sqdist);
		}
	return ans;
}

double STO_3G::integral_coulombic(Nuclei* c, STO_3G* b)
{

	double ans = 0;

	Vector3d vdist_ab = m_center->getPosition() - b->getCenter()->getPosition();
	double sqdist_ab = pow(vdist_ab.norm(), 2);

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			double A = a[i], B = b->getExponent()[j];
			Vector3d Rp = (A*m_center->getPosition()+B*b->getCenter()->getPosition())/(A+B);
			
			Vector3d vdist_pc = Rp - c->getPosition();
			double sqdist_pc = pow(vdist_pc.norm(), 2);
			
			//FIRST CONSIDER POLE OF THE F0 FUNCTION AR ZERO
			if(m_center==c && b->getCenter()==c)
			{
				ans += d[i] * b->getCoefficient()[j] * sqrt(sqrt(pow(4*A*B/pow(M_PI, 2), 3))) * (-1)*2*M_PI/(A+B)*c->getCharge()*exp(-A*B/(A+B)*sqdist_ab);
				continue;
			}
			ans += d[i] * b->getCoefficient()[j] * sqrt(sqrt(pow(4*A*B/pow(M_PI, 2), 3))) * (-1)*2*M_PI/(A+B)*c->getCharge()*exp(-A*B/(A+B)*sqdist_ab) * 1/2*sqrt(M_PI/(A+B)/sqdist_pc) * erf(sqrt((A+B)*sqdist_pc));
		}

	return ans;
}

#endif
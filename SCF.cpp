/*
RESTRICTED CLOSED SHELL HARTREE FOCK SELF CONSISTENT FIELD METHOD
WITH STO_3G BASIS
*/
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Nuclei.h"
#include "STO_3G.h"

using namespace std;
using namespace Eigen;

double integral_two(STO_3G* a, STO_3G* b, STO_3G* c, STO_3G* d);
MatrixXd SCF(int curr_iteration, int NUM_NUCLEI, Nuclei** n, int NUM_ELECTRON, STO_3G** basis, MatrixXd init_P, bool isConverge);

int main()
{
	/*
	1. INITIALIZATION OF THE NUCLEI LIST
	*/
	int NUM_NUCLEI = 0;
	cout << "\nInput number of nuclei: ";
	cin >> NUM_NUCLEI;

	Nuclei** n = new Nuclei*[NUM_NUCLEI];
	int input_charge = 0;
	Vector3d input_position(0, 0, 0);

	for(int i = 0; i < NUM_NUCLEI; i++)
	{
		cout << "\nInput charge of nucleus " << i+1 << ": ";
		cin >> input_charge;

		cout << "Input x coordinate of nucleus " << i+1 << ": ";
		cin >> input_position[0];
		cout << "Input y coordinate of nucleus " << i+1 << ": ";
		cin >> input_position[1];
		cout << "Input z coordinate of nucleus " << i+1 << ": ";
		cin >> input_position[2];

		n[i] = new Nuclei(input_charge, input_position);
	}
	
	/*
	2. SPECIFICATION OF THE STO_3G LIST
	*/
	int NUM_ELECTRON = 0;
	cout << "\nInput number of electrons: ";
	cin >> NUM_ELECTRON;

	STO_3G** basis = new STO_3G*[NUM_ELECTRON];
	int input_center = 0;
	double input_zeta = 0;

	for(int i = 0; i < NUM_ELECTRON; i++)
	{
		cout<<"\nInput the center nucleus of electron "<< i+1 <<" (1, 2, ...): ";
		cin>>input_center;
		cout<<"\nInput the Slater exponent of electron "<< i+1 <<": ";
		cin>>input_zeta;

		basis[i] = new STO_3G(input_zeta, n[input_center-1]);
	}

	SCF(1, NUM_NUCLEI, n, NUM_ELECTRON, basis, MatrixXd::Zero(NUM_ELECTRON, NUM_ELECTRON), false);

	/*
	RELEASE MEMORY OF DYNAMICALLY CREATED OBJECTS
	*/
	for(int i = 0; i<NUM_NUCLEI; i++)
		delete n[i];
	for(int i = 0; i<NUM_ELECTRON; i++)
		delete basis[i];
	return 0;
}

MatrixXd SCF(int curr_iteration, int NUM_NUCLEI, Nuclei** n, int NUM_ELECTRON, STO_3G** basis, MatrixXd init_P, bool isConverge)
{
	if(isConverge||curr_iteration>10)
	{
		cout<<"\nConverged! The density matrix P is:\n"<<init_P<<endl<<endl;
		return init_P;
	}

	cout<<"\n***********ITERATION NO."<<curr_iteration<<"***********\n";
	/*
	3. CONSTRUCT THE OVERLAP MATRIX S, THE KINETIC MATRIX T, AND THE COULOMBIC MATRICES Vi, 1<=i<=NUM_NUCEI
	*/

	MatrixXd S(NUM_ELECTRON, NUM_ELECTRON);
	MatrixXd T(NUM_ELECTRON, NUM_ELECTRON);
	MatrixXd** V_p = new MatrixXd*[NUM_NUCLEI];

	for(int i = 0; i < NUM_NUCLEI; i++)
		V_p[i] = new MatrixXd(NUM_ELECTRON, NUM_ELECTRON);

	for(int i = 0; i<NUM_ELECTRON; i++)
		for(int j = 0; j<NUM_ELECTRON; j++)
			{
				S(i, j) = basis[i]->integral_overlap(basis[j]);
				T(i, j) = basis[i]->integral_kinetic(basis[j]);

				for(int k = 0; k < NUM_NUCLEI; k++)
					V_p[k]->operator()(i, j) = basis[i]->integral_coulombic(n[k], basis[j]);
			}
	cout<<"\nThe overlap matrix S is:\n"<<S<<endl;
	//cout<<"\nThe kinetic matrix T is:\n"<<T<<endl;
	//for(int i = 0; i < NUM_NUCLEI; i++)
	//	cout<<"\nThe coulombic matrix V"<<i+1<<" is:\n"<<*(V_p[i])<<endl;
	
	/*
	4. COMBINE THE ONE ELECTRON INTEGRALS INTO THE CORE MATRIX
	*/
	MatrixXd Hcore = T;

	for(int i = 0; i < NUM_NUCLEI; i++)
		Hcore += *(V_p[i]);

	cout<<"\nThe core matrix Hcore is:\n"<<Hcore<<endl;

	/*
	5. CALCULATE G FROM DENSITY MATRIX P AND TWO ELECTRON INTEGRALS AND 
		ADD TO HCORE TO OBTAIN THE FOCK OPERATOR UNDER THE CURRENT BASIS
	*/
	MatrixXd G = MatrixXd::Zero(NUM_ELECTRON, NUM_ELECTRON);

	for(int i=0; i<NUM_ELECTRON; i++)
		for(int j =0; j<NUM_ELECTRON; j++)
			for(int k = 0; k<NUM_ELECTRON; k++)
				for(int l=0; l<NUM_ELECTRON; l++)
					G(i, j) += init_P(k, l)*(integral_two(basis[i], basis[j], basis[l], basis[k])
						-0.5*integral_two(basis[i], basis[k], basis[l], basis[j]));
	cout<<"\nThe two-electron G is:\n"<<G<<endl;

	MatrixXd F = Hcore + G;

	cout<<"\nThe Fock operator F is:\n"<<F<<endl;
	/*
	6. OBTAIN THE CANONICAL ORTHOGONAL TRANSFORMATION MATRIX VIA CANONICAL DIAGONALIZATION
	*/
	SelfAdjointEigenSolver<MatrixXd> S_solver(S);
	//MatrixXd s = S_solver.eigenvalues().asDiagonal();
	MatrixXd U = S_solver.eigenvectors();
	//MatrixXd s_invsqrt = U.inverse()*S_solver.operatorInverseSqrt()*U;
	MatrixXd X_canonical = S_solver.operatorInverseSqrt()*U;
	//cout<<"\nCanonical orthogonalization transformation matrix is:\n"<<X_canonical<<endl;

	/*
	7. DIAGONALIZE THE TRANSFORMED FOCK OPERATOR TO GET DENSITY MATRIX P
	*/
	SelfAdjointEigenSolver<MatrixXd> F_solver(X_canonical.transpose()*F*X_canonical);
	MatrixXd C = X_canonical*F_solver.eigenvectors();
	MatrixXd P = 2*C.leftCols(NUM_ELECTRON/2)*C.leftCols(NUM_ELECTRON/2).transpose();
	cout<<"\nThe density matrix P is:\n"<<P<<endl
		<<"\n**********END OF ITERATION**********\n";

	/*
	RELEASE MEMORY OF DYNAMICALLY CREATED OBJECTS
	*/
	for(int i = 0; i<NUM_NUCLEI; i++)
		delete V_p[i];

	/*
	8. NEXT ITERATION
	*/
	return SCF(curr_iteration+1, NUM_NUCLEI, n, NUM_ELECTRON, basis, P, P==init_P);
}
/*
STO_3G TWO ORBITAL INTEGRALS
*/
double integral_two(STO_3G* a, STO_3G* b, STO_3G* c, STO_3G* d)
{
	double ans = 0;

	Vector3d vdist_ab = a->getCenter()->getPosition() - b->getCenter()->getPosition();
	double sqdist_ab = pow(vdist_ab.norm(), 2);

	Vector3d vdist_cd = c->getCenter()->getPosition() - d->getCenter()->getPosition();
	double sqdist_cd = pow(vdist_cd.norm(), 2);

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			for(int k = 0; k < 3; k++)
				for(int l = 0; l < 3; l++)
				{
					double A = a->getExponent()[i], B = b->getExponent()[j], C = c->getExponent()[k], D = d->getExponent()[l];

					Vector3d Rp = (A*a->getCenter()->getPosition()+B*b->getCenter()->getPosition())/(A+B);
					Vector3d Rq = (C*c->getCenter()->getPosition()+D*d->getCenter()->getPosition())/(C+D);

					Vector3d Rpq = Rp-Rq;
					double sqdist_pq = pow(Rpq.norm(), 2);

					//FIRST CONSIDER POLE OF THE F0 FUNCTION AR ZERO
					if(sqdist_pq==0)
					{
						ans += a->getCoefficient()[i] * b->getCoefficient()[j] * c->getCoefficient()[k] * d->getCoefficient()[l]
							* sqrt(sqrt(pow(16*A*B*C*D/pow(M_PI, 4), 3)))
							* 2*sqrt(pow(M_PI, 5))/((A+B)*(C+D)*sqrt(A+B+C+D))
							* exp(-A*B/(A+B)*sqdist_ab - C*D/(C+D)*sqdist_cd);
						continue;
					}

					ans += a->getCoefficient()[i] * b->getCoefficient()[j] * c->getCoefficient()[k] * d->getCoefficient()[l]
							* sqrt(sqrt(pow(16*A*B*C*D/pow(M_PI, 4), 3)))
							* 2*sqrt(pow(M_PI, 5))/((A+B)*(C+D)*sqrt(A+B+C+D))
							* exp(-A*B/(A+B)*sqdist_ab - C*D/(C+D)*sqdist_cd)
							* 1/2*sqrt(M_PI/((A+B)*(C+D)/(A+B+C+D)*sqdist_pq))
							* erf(sqrt((A+B)*(C+D)/(A+B+C+D)*sqdist_pq));
				}

	return ans;
}











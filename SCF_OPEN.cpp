/*
UNRESTRICTED OPEN SHELL HARTREE FOCK SELF CONSISTENT FIELD METHOD

PARAMETER INPUT PROCEDURE
1. SPECIFY THE NUCLEI NUMBER AND CHARGE.
2. SPECIFY THE NUMBER OF ELECTRONS ON THE CURRENT NUCLEUS.
3. SPECIFY THE SPIN OF EACH ELECTRON.
4. SPECIFY THE NUMBER OF BASIS TRIAL FUNCTIONS ON THE CURRENT NUCLEUS.
	*IF SMALLER THAN THE AMOUNT OF ELECTRONS ON THAT NUCLEUS, THROW ERROR
5. INPUT THE GUESS DENSITY MATRIX.

PRE-SCF PROCEDURE TO MAXIMIZE PERFORMANCE
1. CALCULATE THE OVERLAP MATRIX.
2. STORE THE CANONICAL ORTHOGONAL TRANSFORMATION MATRIX AS X_CANONICAL.
3. STORE THE KINETIC MATRIX AND COULOMBIC MATRIX AS HCORE.


SELF CONSSITENT FIELD PROCEDURE
1. CALCULATE THE MATRIX G FOR EACH SPIN FROM THE GUESS DENSITY MATRIX.
2. CALCULATE THE FOCK MATRIX FOR EACH SPIN FROM G AND HCORE.
3. DIAGONALIZE THE TRANSFORMED FOCK MATRIX AND OBTAIN A SET OF NEW
	GUESS DENSITY MATRICES.
4. CHECK CONVERGENCE, ENTER THE NEXT ITERATION.

GUAN WANG
JULY 2018
*/

#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "NUCLEUS.h"
#include "BASIS.h"

using namespace std;
using namespace Eigen;

MatrixXd SCF(int curr_iteration, vector<NUCLEI>& n, vector<BASIS>& basis, 
				int NUM_A_ELECTRON, int NUM_B_ELECTRON, MatrixXd& init_Pa, 
				MatrixXd& init_Pb, MatrixXd& X_canonical, MatrixXd& Hcore);

int main()
{
	cout<<"*******************************************************\nOPEN SHELL HARTREE FOCK SELF CONSISTENT FIELD PROGRAM\nv1.0 by GUAN WANG\n*******************************************************\n";
	vector<NUCLEI> n;
	vector<BASIS> basis;

	int input_charge = 0;
	Vector3d input_position(0, 0, 0);
	int input_size = 0;
	double input_zeta = 0;
	char ans = 'Y';
	int NUM_A_ELECTRON = 0;
	int NUM_B_ELECTRON = 0;

	for(; ans == 'Y' || ans == 'y'; cin>>ans)
	{

		cout<<"\nInput nucleus charge: ";
		cin>>input_charge;

		cout<<"\nInput x-coordinate of current nucleus: ";
		cin>>input_position[0];

		cout<<"\nInput y-coordinate of current nucleus: ";
		cin>>input_position[1];

		cout<<"\nInput z-coordinate of current nucleus: ";
		cin>>input_position[2];

		n.push_back(*(new NUCLEI(input_charge, input_position)));

		cout<<"\nInput the size of basis set to be used on this nucleus: ";
		cin>>input_size;

		for(int curr_orbital = 1; curr_orbital < input_size + 1; curr_orbital++)
		{
			cout<<"\nInput the Slater exponent for orbital "<<curr_orbital<<": ";
			cin>>input_zeta;

			basis.push_back(*(new STO_3G(curr_orbital, input_zeta, n[n.size()-1])));
		}

		cout<<"\nAdd nucleus? (Y/N) ";
	}

	cout<<"\nInput the number of alpha-electrons: ";
	cin>>NUM_A_ELECTRON;

	cout<<"\nInput the number of beta-electrons: ";
	cin>>NUM_B_ELECTRON;

	/*
	1. ENTER THE INITIAL GUESS DENSITY MATRICES FOR RESPECTIVE SPINS
	2. SET THE CANONICAL TRANSFORMATION MATRIX AND THE CORE HAMILTONIAN MATRIX TO ZERO
	*/
	MatrixXd init_Pa = MatrixXd::Zero(basis.size(), basis.size());
	MatrixXd init_Pb = MatrixXd::Zero(basis.size(), basis.size());

	MatrixXd X_canonical = MatrixXd::Zero(basis.size(), basis.size());
	MatrixXd Hcore = MatrixXd::Zero(basis.size(), basis.size()); 
	
	SCF(1, n, basis, NUM_A_ELECTRON, NUM_B_ELECTRON, init_Pa, init_Pb, X_canonical, Hcore);
}

MatrixXd SCF(int curr_iteration, vector<NUCLEI>& n, vector<BASIS>& basis, 
				int NUM_A_ELECTRON, int NUM_B_ELECTRON, MatrixXd& init_Pa, 
				MatrixXd& init_Pb, MatrixXd& X_canonical, MatrixXd& Hcore)
{
	cout<<"\n*********************ITERATION NO."<<curr_iteration<<"*********************\n";

	if(curr_iteration == 1)
	{
		MatrixXd S(basis.size(), basis.size());
		MatrixXd T(basis.size(), basis.size());
		MatrixXd** V_p = new MatrixXd*[n.size()];

		for(int i = 0; i < n.size(); i++)
			V_p[i] = new MatrixXd(basis.size(), basis.size());

		for(int i = 0; i < basis.size(); i++)
			for(int j = 0; j < basis.size(); j++)
				{
					S(i, j) = integral_overlap(basis[i], basis[j]);
					T(i, j) = integral_kinetic(basis[i], basis[j]);

					for(int k = 0; k < n.size(); k++)
						V_p[k]->operator()(i, j) = integral_coulombic(basis[i], basis[j], n[k]);
				}
		/*
		COMBINE THE ONE ELECTRON INTEGRALS INTO THE CORE MATRIX
		*/
		Hcore = T;

		for(int i = 0; i < n.size(); i++)
			Hcore += *(V_p[i]);

		cout<<"\nThe overlap matrix S is:\n"<<S<<endl;
		cout<<"\nThe kinetic matrix T is:\n"<<T<<endl;
		cout<<"\nThe core matrix Hcore is:\n"<<Hcore<<endl;
		/*
		OBTAIN THE CANONICAL ORTHOGONAL TRANSFORMATION MATRIX VIA CANONICAL DIAGONALIZATION
		*/
		SelfAdjointEigenSolver<MatrixXd> S_solver(S);
		MatrixXd U = S_solver.eigenvectors();
		X_canonical = S_solver.operatorInverseSqrt()*U;

		cout<<"\nThe canonical transform Xcanonical is:\n"<<X_canonical<<endl;
		/*
		RELEASE MEMORY OF DYNAMICALLY CREATED OBJECTS
		*/
		for(int i = 0; i<n.size(); i++)
			delete V_p[i];
	}

	MatrixXd Ga = MatrixXd::Zero(basis.size(), basis.size());
	for(int i = 0; i < basis.size(); i++)
		for(int j = 0; j < basis.size(); j++)
			for(int k = 0; k < basis.size(); k++)
				for(int l = 0; l < basis.size(); l++)
					Ga(i, j) += (init_Pa(k, l)+init_Pb(k, l))*(integral_two(basis[i], basis[j], basis[l], basis[k])
						-init_Pa(k, l)*integral_two(basis[i], basis[k], basis[l], basis[j]));

	MatrixXd Gb = MatrixXd::Zero(basis.size(), basis.size());
	for(int i = 0; i < basis.size(); i++)
		for(int j = 0; j < basis.size(); j++)
			for(int k = 0; k < basis.size(); k++)
				for(int l = 0; l < basis.size(); l++)
					Ga(i, j) += (init_Pa(k, l)+init_Pb(k, l))*(integral_two(basis[i], basis[j], basis[l], basis[k])
						-init_Pb(k, l)*integral_two(basis[i], basis[k], basis[l], basis[j]));

	MatrixXd Fa = Hcore + Ga;
	cout<<"\nThe fock matrix for spin alpha is:\n"<<Fa<<endl;
	MatrixXd Fb = Hcore + Gb;
	cout<<"\nThe fock matrix for spin beta is:\n"<<Fb<<endl;


	SelfAdjointEigenSolver<MatrixXd> Fa_solver(X_canonical.transpose()*Fa*X_canonical);
	MatrixXd Ca = X_canonical*Fa_solver.eigenvectors();
	MatrixXd Pa = Ca.leftCols(NUM_A_ELECTRON)*Ca.leftCols(NUM_A_ELECTRON).transpose();

	SelfAdjointEigenSolver<MatrixXd> Fb_solver(X_canonical.transpose()*Fb*X_canonical);
	MatrixXd Cb = X_canonical*Fb_solver.eigenvectors();
	MatrixXd Pb = Cb.leftCols(NUM_B_ELECTRON)*Cb.leftCols(NUM_B_ELECTRON).transpose();

	cout<<"\nThe net density matrix is:\n"<<Pa+Pb<<endl
		<<"\n********************END OF ITERATION********************\n";

	if((Pa == init_Pa && Pb == init_Pb) || curr_iteration > 1)
	{
		cout<<"Converged!";
		return Pa+Pb;
	}

	return SCF(curr_iteration + 1, n, basis, NUM_A_ELECTRON, NUM_B_ELECTRON, Pa, Pb,
				X_canonical, Hcore);
}
















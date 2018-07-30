/*
NUCLEI CLASS
SPECIFIES CHARGE AND POSITION
*/

#ifndef NUCLEI_H
#define NUCLEI_H

#include "Eigen/Dense"

using namespace Eigen;

class NUCLEI
{
public:
	NUCLEI(int input_charge, Vector3d input_position)
			: m_charge(input_charge), m_position(input_position){}
	int getCharge();
	Vector3d getPosition();
private:
	int m_charge, total_center, m_orbitals;
	Vector3d m_position;
};


int NUCLEI::getCharge(){return m_charge;}

Vector3d NUCLEI::getPosition(){return m_position;}

#endif
/*
NUCLEI CLASS
SPECIFIES CHARGE AND POSITION
*/

#ifndef NUCLEI_H
#define NUCLEI_H

#include "Eigen/Dense"

using namespace Eigen;

class Nuclei
{
public:
	Nuclei(int, Vector3d);
	int getCharge();
	Vector3d getPosition();
private:
	int m_charge, total_center;
	Vector3d m_position;
};

Nuclei::Nuclei(int input_charge, Vector3d input_position)
{
	m_charge = input_charge;
	m_position = input_position;
}

int Nuclei::getCharge(){return m_charge;}

Vector3d Nuclei::getPosition(){return m_position;}

#endif
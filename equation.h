#ifndef __EQUATION__H_
#define __EQUATION__H_

#include "scikit.h"
#include "line.h"
#include "point.h"

class Equation
{
public:
	Equation(){};
	void setEquation(char*,int);
	float* getCoefficients();
	int  getDegree();
	float calculate_y_value(float );
private:
	int deg;
	float* coeff;
	char* denklem;
	void Calculate();
}; 
#endif

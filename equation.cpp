#include "equation.h"
#include <iostream>
#include "scikit.h"
#include "line.h"
#include "point.h"

using namespace std;


int Equation::getDegree()
{
	return deg;
}

void Equation::setEquation(char* denklem,int deg)
{
	this->deg = deg+1;
	string s = denklem;
	s.append("+");
	denklem = &s[0];
	this->denklem = denklem;
	coeff = new float[deg];
	for(int j=0;j<deg+1;j++)
	{
		coeff[j] = 0;
	}
	Calculate();
	
}
float* Equation::getCoefficients()
{
	float* ret = new float[deg];
	for(int i=0;i<deg;i++)
	{
		ret[i] = coeff[i];
	}
	return ret;
}
void Equation::Calculate()//x^2+2^x+1x^0 formatýnda girilecek x2+2x+1
{
	int ind =0,len = strlen(denklem),tempind =0,degind=0,tempdegind =0;
	char* temp;
	float value = 0.0;
	char* tempdeg;
	for(int i=0;i<len;i++)
	{
		if(denklem[i] == 'x')
		{
			if(denklem[i+1] == '^')
			{		
				degind = i+2;
				temp = new char[i-ind+1];
				for(int j=ind;j<i;j++)
				{
					temp[tempind] = denklem[j];
					tempind++; 	
				}
				if(denklem[ind-1] == '-')
				{
					value = atof(temp)*-1;
				}
				else
					value = atof(temp);	
				tempind =0;	
			}
			else
				return;
		}
				
		if(denklem[i] == '+' || denklem[i] == '-')
		{
			ind = i+1;
			tempdeg = new char[i-degind+1];
			for(int j=degind;j<i;j++)
			{
				tempdeg[tempdegind] = denklem[j];
				tempdegind++;
			}
		
			coeff[atoi(tempdeg)] = value;
			tempdegind=0;
		
	
		}
	
}

	if(denklem[0] == '-')
		coeff[deg-1]*=-1;
	delete [] temp;
	delete [] tempdeg;

}
float Equation::calculate_y_value(float x)
{
	float y = 0;
	for(int i=0;i<deg;i++)
	{
		y+=pow(x,i)*coeff[i];
	}
	return y;
}

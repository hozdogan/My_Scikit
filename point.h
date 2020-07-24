#ifndef __POINT__H_
#define __POINT__H_

#include "scikit.h"
#include "line.h"

#include "equation.h"
class Point
{

public:
	Point(){};
	Point(double ,double);
	void setPoints(double,double);
	double getX();
	double getY();
private:
	double x,y;
		
	
};

#endif

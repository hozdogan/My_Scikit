#ifndef __LINE__H_
#define __LINE__H_

#include "scikit.h"

#include "point.h"
#include "equation.h"

class Scikit;
class Point;
class Line
{
	public:
		Line(){};
		Line(char*);
		void setLine(Point,Point);
		char* findLine_equation();
		Point* getLine_point1();
		Point* getLine_point2();
		double get_Slope();
		
		
	private:
		double slope;
	    Point *p1,*p2;//bunlarý pointer geçmezsen hata veriyor baska bir classýn objesi olarak geçemezsin
		Scikit *sk;	
	
};


#endif

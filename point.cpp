#include "point.h"


Point* Line::getLine_point1()
{
	return p1;
}
void Line::setLine(Point p1,Point p2)
{
	this->p1 = p1;
	this->p2 = p2;	
}
Point* Line::getLine_point2()
{
	return p2;
}

Point::Point(double x,double y)
{
	this->x = x;
	this->y = y;
}
void Point::setPoints(double x,double y)
{
	this->x = x;
	this->y = y;
}
double Point::getX()
{
	return this->x;
}
double Point::getY()
{
	return this->y;
}

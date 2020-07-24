#include "line.h"

Line::Line(char* line_equation)//1000x+1000
{
	srand(time(NULL));
	int len = strlen(line_equation);
	char* c = new char[2];
	float coeff[2];
	for(int j=0;j<len;j++)
	{
		if(line_equation[j] == 'x')
		{
			c = new char[j];
			for(int k=0;k<j;k++)
			{
				c[k] = line_equation[k];
			}
			coeff[0] = atof(c);
		
		}	
		if(line_equation[j] == '+')
		{
			c = new char[len-j];
			for(int k=j;k<len;k++)
			{
				c[k-j] = line_equation[k]; 
			}
			coeff[1] = atof(c);
		}
	}
	int x1 = 1+rand()%100;
	int x2= 1+rand()%100;	
	p1.setPoints(x1,(coeff[0]*x1)+coeff[1]);
	p2.setPoints(x2,(x2*coeff[0])+coeff[1]);
	slope = (double)coeff[0];
	delete c;
}
char* Line::findLine_equation()
{
	slope = sk.slope_point(p1,p2);
	char* coeff_a = sk.int_to_str(slope);
	int b = (int)p1.getY()-(slope*p1.getX());
	char* coeff_b = sk.int_to_str(b);
	int len_a = strlen(coeff_a),len_b = strlen(coeff_b);
	char* c = new char[len_a+len_b+2];
	for(int j=0;j<len_a;j++)
	{
		c[j] = coeff_a[j];
	}
	c[len_a] = 'x';
	c[len_a+1] = '+';

	for(int j=len_a+2;j<len_a+len_b+2;j++)
	{
		c[j] = coeff_b[j-(len_b+2)];
	}
	return c;
}
double Line::get_Slope()
{
	return slope;
}

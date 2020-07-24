#include <iostream>
#include<math.h>
#include<string.h>
#include <stdlib.h>
#include <time.h>
#define PI 3.14159
#define E 2.71828
using namespace std;
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




class Line;//adamdan gormustum
class Equation;

class Scikit
{
	public:
		Scikit(){};
		double distance_two_point(Point ,Point);
		double slope_point(Point,Point);	
		char* int_to_str(int );
		bool isParallel(Line,Line);
		bool isVertical(Line,Line);
		float* equation_roots(Equation);
		float* derive_equation(int ,float*);
		float calculate_y_value(float,float*,int);
		float* derive_integrate(int ,float*);
		//matrix operations
		float** matrix_multiplecation(float** ,float**,int,int,int);
		float** matrix_addition(float**,float**,int ,int );
		float** matrix_substition(float**,float**,int ,int );
		float** matrix_transform2D_rotating(float**,double,int);
		float** matrix_transpoze(float**,int,int);
		float matrix_trace(float**,int);
		float matrix_determinant(float**,int);
		float** matrix_inverse(float**,int);
		float* gramer_rule(float**,float**,int);
		float* eigen_values(float**,int);
		float** matrix_transform3D_translate(float ,float ,float ,float ,float ,float ); 
		float** matrix_transform3D_rotate(float,float,float,int,float);
		float** matrix_transform3D_scale(float,float,float,float,float,float);
		
		float newton_raphson_square_find(float);
		float correlation_coefficient(float*,float*,int );
		
		//turev islemleri
		float tan_turev(Equation,float);
		float sin_turev(Equation,float);
		float cos_turev(Equation,float);
		float log_turev(Equation,float,int);
		float pow_turev(Equation,float,int);
		float multiple_turev(Equation,Equation,float);
		float divide_turev(Equation,Equation,float);
		float arcsin_turev(Equation,float);
		float arccos_turev(Equation,float);
		float arctan_turev(Equation,float);
		float exponantial_turev(Equation,float,float);
		
		//integral iþlemleri
	
};








/*//class line*/


class Line
{
	public:
		Line(){};
		Line(char*);
		void setLine(Point,Point);
		char* findLine_equation();
		Point getLine_point1();
		Point getLine_point2();
		double get_Slope();
		
		
	private:
		double slope;
		Point p1,p2;
		Scikit sk;	
	
};
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




Point Line::getLine_point1()
{
	return p1;
}
void Line::setLine(Point p1,Point p2)
{
	this->p1 = p1;
	this->p2 = p2;	
}
Point Line::getLine_point2()
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



float Scikit::tan_turev(Equation eq,float value)
{
	int angle = eq.getDegree();
	if(angle > 0)
	{
		angle-=1;
	}
	float t0_value = eq.calculate_y_value(value);
	float res;
	float* deg = eq.getCoefficients();
	float *t1_deg = Scikit::derive_equation(eq.getDegree(),deg);
	float t1_value = Scikit::calculate_y_value(value,t1_deg,angle);
	
	res = t1_value*(1+(tan(t0_value*PI/180)*tan(t0_value*PI/180)));
	return res;

}

float Scikit::pow_turev(Equation eq,float value,int exp)
{
	int angle = eq.getDegree();
	if(angle > 0)
	{
		angle-=1;
	}

	float t0_value = eq.calculate_y_value(value);
	float res;
	float* deg = eq.getCoefficients();
	float *t1_deg = Scikit::derive_equation(eq.getDegree(),deg);
	float t1_value = Scikit::calculate_y_value(value,t1_deg,angle);
	
	res = exp*t1_value*powf(t0_value,exp-1);
	return res;
	
}

float Scikit::multiple_turev(Equation eq1,Equation eq2,float value)
{
	int angle = eq1.getDegree();
	float res;
	if(angle > 0)
	{
		angle-=1;
	}

	float t00_value = eq1.calculate_y_value(value);
	float* deg1 = eq1.getCoefficients();
	float *t01_deg = Scikit::derive_equation(eq1.getDegree(),deg1);
	float t01_value = Scikit::calculate_y_value(value,t01_deg,angle);
	
	float t10_value = eq2.calculate_y_value(value);
	float* deg2 = eq2.getCoefficients();
	float *t11_deg = Scikit::derive_equation(eq2.getDegree(),deg2);
	float t11_value = Scikit::calculate_y_value(value,t11_deg,angle);
	
	res = (t01_value*t10_value)+(t11_value*t00_value);
	return res;
	
}

float Scikit::divide_turev(Equation eq1,Equation eq2,float value)
{
	int angle = eq1.getDegree();
	float res;
	if(angle > 0)
	{
		angle-=1;
	}

	float t00_value = eq1.calculate_y_value(value);
	float* deg1 = eq1.getCoefficients();
	float *t01_deg = Scikit::derive_equation(eq1.getDegree(),deg1);
	float t01_value = Scikit::calculate_y_value(value,t01_deg,angle);
	
	float t10_value = eq2.calculate_y_value(value);
	float* deg2 = eq2.getCoefficients();
	float *t11_deg = Scikit::derive_equation(eq2.getDegree(),deg2);
	float t11_value = Scikit::calculate_y_value(value,t11_deg,angle);
	
	res = ((t01_value*t10_value)-(t11_value*t00_value))/(t10_value*t10_value);
	return res;
}

float Scikit::arcsin_turev(Equation eq,float value)
{
	int angle = eq.getDegree();
	if(angle > 0)
	{
		angle-=1;
	}

	float t0_value = eq.calculate_y_value(value);
	float res;
	float* deg = eq.getCoefficients();
	float *t1_deg = Scikit::derive_equation(eq.getDegree(),deg);
	float t1_value = Scikit::calculate_y_value(value,t1_deg,angle);
	
	res = t1_value/(sqrt(1-(t0_value*t0_value)));
}

float Scikit::arccos_turev(Equation eq,float value)
{
	int angle = eq.getDegree();
	if(angle > 0)
	{
		angle-=1;
	}

	float t0_value = eq.calculate_y_value(value);
	float res;
	float* deg = eq.getCoefficients();
	float *t1_deg = Scikit::derive_equation(eq.getDegree(),deg);
	float t1_value = Scikit::calculate_y_value(value,t1_deg,angle);
	
	res = -t1_value/(sqrt(1-(t0_value*t0_value)));
}

float Scikit::arctan_turev(Equation eq,float value)
{
	int angle = eq.getDegree();
	if(angle > 0)
	{
		angle-=1;
	}

	float t0_value = eq.calculate_y_value(value);
	float res;
	float* deg = eq.getCoefficients();
	float *t1_deg = Scikit::derive_equation(eq.getDegree(),deg);
	float t1_value = Scikit::calculate_y_value(value,t1_deg,angle);
	
	res = t1_value/(1+(t0_value*t0_value));
}

float Scikit::cos_turev(Equation eq,float value)
{
	int angle = eq.getDegree();
	if(angle > 0)
	{
		angle-=1;
	}
	
	float t0_value = eq.calculate_y_value(value);
	float res;
	float* deg = eq.getCoefficients();
	float *t1_deg = Scikit::derive_equation(eq.getDegree(),deg);
	float t1_value = Scikit::calculate_y_value(value,t1_deg,angle);
	
	res = t1_value*cos(t0_value*PI/180);
	return res;

}
float Scikit::sin_turev(Equation eq,float value)
{
	
	int angle = eq.getDegree();
	if(angle > 0)
	{
		angle-=1;
	}
	
	float t0_value = eq.calculate_y_value(value);
	float res;
	float* deg = eq.getCoefficients();
	float *t1_deg = Scikit::derive_equation(eq.getDegree(),deg);
	float t1_value = Scikit::calculate_y_value(value,t1_deg,angle);
	
	res = -t1_value*sin(t0_value*PI/180);
	return res;
	
}

float Scikit::log_turev(Equation eq,float value,int base)
{
	int angle = eq.getDegree();
	if(angle > 0)
	{
		angle-=1;
	}

	float t0_value = eq.calculate_y_value(value);
	float res;
	float* deg = eq.getCoefficients();
	float *t1_deg = Scikit::derive_equation(eq.getDegree(),deg);
	float t1_value = Scikit::calculate_y_value(value,t1_deg,angle);
	
	if(base == 2)
	{
		res = (t1_value/t0_value)*log2f(E);
	}	
	else if(base == 10)
	{
			res = (t1_value/t0_value)*logf(E);
	}
	return res;
}


float** Scikit::matrix_transform3D_scale(float x,float y,float z,float sx,float sy,float sz)
{
	float** Vec = new float*[3];
	for(int i=0;i<3;i++)
	{
		Vec[i] = new float[1];
	}
	Vec[0][0] = x;
	Vec[1][0] = y;
	Vec[2][0] = z;
	
	float** valVec = new float*[3];
	for(int i=0;i<3;i++)
	{
		valVec[i] = new float[1];
	}
	valVec[0][0] = sx;
	valVec[1][0] = sy;
	valVec[2][0] = sz;

	
	float** res = new float*[3];
	for(int i=0;i<3;i++)
	{
		res[i] = new float[1];
	}
	
	float** transform_matrix = new float*[3];
	for(int i=0;i<3;i++)
	{
		transform_matrix[i] = new float[3];
	}
	for(int row=0;row<3;row++)
	{
		for(int col=0;col<3;col++)
		{
		   if(row == col)
		   {
		   		transform_matrix[row][col] = valVec[row][0];
		   }
		   else
		   	   transform_matrix[row][col] = 0.0f;
		}
	}
	
	res = Scikit::matrix_multiplecation(transform_matrix,Vec,3,3,1);
	return res;
	
}


float** Scikit::matrix_transform3D_rotate(float x,float y,float z,int axis,float angle)
{
	float cosT = cos(angle*PI/180);
	float sinT = sin(angle*PI/180);
	
	float** Vec = new float*[3];
	for(int i=0;i<3;i++)
	{
		Vec[i] = new float[1];
	}
	Vec[0][0] = x;
	Vec[1][0] = y;
	Vec[2][0] = z;
	
	float** res = new float*[3];
	for(int i=0;i<3;i++)
	{
		res[i] = new float[1];
	}
	
	float** transform_matrix = new float*[3];
	for(int i=0;i<3;i++)
	{
		transform_matrix[i] = new float[3];
	}
	
	for(int row=0;row<3;row++)
	{
		for(int col=0;col<3;col++)
		{
		   if(row == col)
		   {
		   		transform_matrix[row][col] = 1.0f;
		   }
		   else
		   	   transform_matrix[row][col] = 0.0f;
		}
	}
	
	if(axis == 0)//x = 0 y = 1 z = 2
	{
		transform_matrix[1][1] = cosT;
		transform_matrix[1][2] = -sinT;
		transform_matrix[2][1] = sinT;
		transform_matrix[2][2] = cosT;
	}
	else if(axis == 1)
	{
		transform_matrix[0][0] = cosT;
		transform_matrix[0][2] = sinT;
		transform_matrix[2][0] = -sinT;
		transform_matrix[2][2] = cosT;
	}
	else if(axis == 2)
	{
		transform_matrix[0][0] = cosT;
		transform_matrix[0][1] = -sinT;
		transform_matrix[1][0] = sinT;
		transform_matrix[1][1] = cosT;
	}
	
	res = Scikit::matrix_multiplecation(transform_matrix,Vec,3,3,1);
	return res;
}


float** Scikit::matrix_transform3D_translate(float x,float y,float z ,float dx,float dy,float dz)
{
	float** Vec = new float*[4];
	for(int i=0;i<4;i++)
	{
		Vec[i] = new float[1];
	}
	
	float** valVec = new float*[4];
	for(int i=0;i<4;i++)
	{
		valVec[i] = new float[1];
	}
	valVec[0][0] = dx;
	valVec[1][0] = dy;
	valVec[2][0] = dz;
	valVec[3][0] = 1.0f;
	
	float** res = new float*[4];
	for(int i=0;i<4;i++)
	{
		res[i] = new float[1];
	}
	
	float** transform_matrix = new float*[4];
	for(int i=0;i<4;i++)
	{
		transform_matrix[i] = new float[4];
	}
	Vec[0][0] = x;
	Vec[1][0] = y;
	Vec[2][0] = z;
	Vec[3][0] = 1.0f;
	
	for(int row=0;row<4;row++)
	{
		for(int col=0;col<4;col++)
		{
		   if(row == col)
		   {
		   		transform_matrix[row][col] = 1.0f;
		   }
		   else if(col == 3)
		       transform_matrix[row][col] = valVec[row][0];
		   else
		   	   transform_matrix[row][col] = 0.0f;
		}
	}
	
	res = Scikit::matrix_multiplecation(transform_matrix,Vec,4,4,1);
	
	return res;
	
}


float* Scikit::eigen_values(float** m1,int n)
{
	float teta,costeta,sinteta,temp1,temp2,temp3,temp4;
	int count,ind = 0;bool progress = true;
	float* eigen_array = new float[n];
	float** temp_mat = new float*[n];
	for(int i=0;i<n;i++)
	{
		temp_mat[i] = new float[n];
	}
			
	
	for(int row=0;row<n;row++)
	{
		for(int col=0;col<n;col++)
		{
		    temp_mat[row][col] = m1[row][col];
		}
	}
	while(progress)
	{
		for(int row=0;row<n;row++)
	{
		for(int col=0;col<n;col++)
		{
		    if(row != col && temp_mat[row][col] > 0.00001)
		    {
		    
				teta = (atan((2*temp_mat[row][col])/(temp_mat[col][col] - temp_mat[row][row]))* (180 / PI))/2+90;//arctan bole tanýmlanacak
				costeta = cos(teta*PI/180);
				sinteta = sin(teta*PI/180);
				
				temp1 = temp_mat[row][row];temp2 = temp_mat[col][col];temp3 = temp_mat[row][col];
				
		    	temp_mat[row][row] = (costeta*costeta*temp1)-(2*sinteta*costeta*temp3)+(sinteta*sinteta*temp2);
		    	temp_mat[col][col] = (sinteta*sinteta*temp1)+(2*sinteta*costeta*temp3)+(costeta*costeta*temp2);
		    	temp_mat[row][col] = (((costeta*costeta) - (sinteta*sinteta))*temp3)+(costeta*sinteta*(temp1-temp2));
		    	temp_mat[col][row] = (((costeta*costeta) - (sinteta*sinteta))*temp3)+(costeta*sinteta*(temp1-temp2));
			}
		}
	}
	for(int row=0;row<n;row++)
	{
		for(int col=0;col<n;col++)
		{
		    if(temp_mat[row][col] < 0.00001 && row != col)
		    {
		    	count++;
			}
		}
	}
	if(count == (n*n)-n)
	{
		for(int row=0;row<n;row++)
		{
			for(int col=0;col<n;col++)
			{
		    	if(row == col)
		    	{
		    		eigen_array[ind] = temp_mat[row][col];
		    		ind++;
				}
			}
		}
		progress = false;
	}
	else
		count=0;
}
	return eigen_array;
	
}



float* Scikit::gramer_rule(float** m1,float** res,int n)
{
	float* root_array = new float[n];//2 pointer diziyi eþitleme temp yapacaksan normal yap
	
	float** temp_matris = new float*[n];
	for(int i=0;i<n;i++)
	{
		temp_matris[i] = new float[n];
	}
	for(int row=0;row<n;row++)
	{
		for(int col=0;col<n;col++)
		{
		    temp_matris[row][col] = m1[row][col];
		}
	}
	
	
	
	for(int row=0;row<n;row++)
	{
		for(int col=0;col<n;col++)
		{
		    temp_matris[col][row] = res[col][0];
		}
		root_array[row] = Scikit::matrix_determinant(temp_matris,n)/Scikit::matrix_determinant(m1,n);
		for(int row=0;row<n;row++)
		{
			for(int col=0;col<n;col++)
			{
		    	temp_matris[row][col] = m1[row][col];
			}
		}
	}
	return root_array;
	
}


float** Scikit::matrix_inverse(float** m1,int n)
{
	float d,k,temp,det;
	float** unit_matris = new float*[n];
	for(int i=0;i<n;i++)
	{
		unit_matris[i] = new float[n];
	}
	if(n < 2)
	{
		exit(-1);
	}
	else if(n == 2)
	{
		det = Scikit::matrix_determinant(m1,2);
		m1[0][1]*=-1;
		m1[1][0]*=-1;
		temp = m1[0][0];
		m1[0][0] = m1[1][1];
		m1[1][1] = temp;
		for(int row=0;row<n;row++)
		{
			for(int col=0;col<n;col++)
			{
				m1[row][col]/=det;			
			}
		}
		return m1;
	}
	else
	{
		for(int row=0;row<n;row++)
		{
			for(int col=0;col<n;col++)
			{
				if(row == col)
				{
					unit_matris[row][col] = 1.0f;
				}
				else
				{
					unit_matris[row][col] = 0.0f;
				}
			}
		}
		for(int i=0;i<n;i++)
		{
			d = m1[i][i];
			for(int j=0;j<n;j++)
			{
				m1[i][j]/=d;
				unit_matris[i][j]/=d;
			}
				for(int x=0;x<n;x++)
				{
					if(x!=i)
					{
						k=m1[x][i];
						for(int l=0;l<n;l++)
						{
							m1[x][l] -= (m1[i][l]*k);
							unit_matris[x][l] -= (unit_matris[i][l]*k);
						}
					}
				}
		}
		
		return unit_matris;	
	}
	
}

float Scikit::matrix_determinant(float** m1,int n)
{
	int count = 0,temp;
	float det=0;
	float* one_dimension_temp = new float[(n-1)*(n-1)];
	float** two_dimension_temp = new float*[n-1];
	for(int j=0;j<n-1;j++)
	{
		two_dimension_temp[j] = new float[n-1];
	}
	
	if(n == 2)
	{
		det = (m1[0][0]*m1[1][1]) - (m1[0][1]*m1[1][0]);
		return det;
	}
	else if(n > 2)
	{
		for(int i=0;i<n;i++)
		{
			for(int row=0;row<n;row++)	
			{
				for(int col=0;col<n;col++)
				{
					if(row == 0 || col == i)
					{
						continue;
					}
					else
					{	
						one_dimension_temp[count] = m1[row][col];
						count++;					
					}
				}		
			}
			count = 0;
			for(int row=0;row<n-1;row++)	
			{
				for(int col=0;col<n-1;col++)
				{
					temp = row*(n-1)+col;
					two_dimension_temp[(temp-col)/(n-1)][col] = one_dimension_temp[temp];		
				}
			}
			if(i%2 == 1)
			{
				m1[0][i]*=-1;
			}
			det += m1[0][i] * matrix_determinant(two_dimension_temp,n-1);
		}
	}
	return det;
}


float Scikit::matrix_trace(float** m1,int n)
{
	float sum = 0.0f;
	for(int row=0;row<n;row++)
	{
		for(int col=0;col<n;col++)
		{
			if(row==col)
			{
				sum+=m1[row][col];
			}
		}
	}
	return sum;
}


float Scikit::newton_raphson_square_find(float sq)
{
	float p1 = 3.0f,p2;
	int count=0;
	Equation eq;
	eq.setEquation("1x^2",2);
	float* coeff = eq.getCoefficients();
	float y_value,derive_y_value;
	float* derive_array = Scikit::derive_equation(eq.getDegree(),coeff);
	
	
	y_value = eq.calculate_y_value(p1);
	derive_y_value = Scikit::calculate_y_value(p1,derive_array,eq.getDegree()-1);
	while(true)
	{
		p2 = (y_value + sq) / derive_y_value;
		if(abs(p1 - p2) < 0.001f )
		{
			return p2;
			break;
		}
		else
		{
			p1 = p2;
			y_value = eq.calculate_y_value(p1);
			derive_y_value = Scikit::calculate_y_value(p1,derive_array,eq.getDegree()-1);
		}
	}
	
}
float Scikit::correlation_coefficient(float* m1,float* m2,int len)
{
	float m1ort,m2ort,r,summ1=0.0f,summ2=0.0f,rx=0.0f,ry=0.0f;
	for(int i=0;i<len;i++)
	{
		summ1+=m1[i];
		summ2+=m2[i];
	}
	m1ort = summ1/len;
	m2ort = summ2/len;
	for(int j=0;j<len;j++)
	{
		rx+=(m1[j]-m1ort)*(m2[j]-m2ort);
		ry+=sqrt(pow((m1[j]-m1ort),2)*pow((m2[j]-m2ort),2));
	}
	r = rx/ry;
	return r;
}


float** Scikit::matrix_multiplecation(float** m1,float** m2,int m,int n,int k)
{
	float** result = new float*[m];
	for(int i=0;i<m;i++)
	{
		result[i] = new float[k];
	}
	for(int row=0;row<m;row++)
	{
		for(int col=0;col<n;col++)
		{
			result[row][col] = 0.0f;
		}
	}
	
	//multiply
	for(int row=0;row<m;row++)
	{
		for(int col=0;col<k;col++)
		{
			for(int op=0;op<n;op++)
			{
				result[row][col]+=m1[row][op]*m2[op][col];
			}
		}
	}
	return result;
}

float** Scikit::matrix_transpoze(float** m1,int m,int n)
{
	float** result = new float*[n];
	for(int i=0;i<m;i++)
	{
		result[i] = new float[m];
	}
	
	for(int row=0;row<m;row++)
	{
		for(int col=0;col<n;col++)
		{
			m1[row][col] = result[col][row];
		}
	}
	return result;
	
}

float** Scikit::matrix_addition(float** m1,float**m2,int m,int n)
{
	float** result = new float*[m];
	for(int i=0;i<m;i++)
	{
		result[i] = new float[n];
	}
	for(int row=0;row<m;row++)
	{
		for(int col=0;col<n;col++)
		{
			result[row][col]=m1[row][col]+m2[row][col];
		}
	}
	return result;
	
}


float** Scikit::matrix_substition(float** m1,float**m2,int m,int n)
{
	float** result = new float*[m];
	for(int i=0;i<m;i++)
	{
		result[i] = new float[n];
	}
	for(int row=0;row<m;row++)
	{
		for(int col=0;col<n;col++)
		{
			result[row][col]=m1[row][col]-m2[row][col];
		}
	}
	return result;
	
}

float** Scikit::matrix_transform2D_rotating(float** m1,double angle,int direction)//1 se saat yönü -1 se tersi
{
	float** result = new float*[2];
	for(int i=0;i<2;i++)
	{
		result[i] = new float[1];
	}
		for(int row=0;row<2;row++)
		{
			for(int col=0;col<1;col++)
			{
				result[row][col] = 0.0f;	
			}
		}
	
	float** rotate_matrix = new float*[2];
	for(int i=0;i<2;i++)
	{
		rotate_matrix[i] = new float[2];
	}
	if(direction == -1)
	{
		rotate_matrix[0][0] = cos(angle);
		rotate_matrix[0][1] = -1*sin(angle);	
		rotate_matrix[1][0] = sin(angle);
		rotate_matrix[1][1] = cos(angle);
	
	}
	else if(direction == 1)
	{
	  	rotate_matrix[0][0] = cos(angle);
		rotate_matrix[0][1] = sin(angle);	
		rotate_matrix[1][0] = -1*sin(angle);
		rotate_matrix[1][1] = cos(angle);
	}
	
	for(int row=0;row<2;row++)
	{
		for(int col=0;col<1;col++)
		{
			for(int op=0;op<2;op++)
			{
				result[row][col] += rotate_matrix[row][op]*m1[op][col];
			}
		}
	}

	return result;
}

double Scikit::distance_two_point(Point p1,Point p2)
{
	double valueX = pow((p2.getX() - p1.getX()),2);
	double valueY = pow((p2.getY() - p1.getY()),2);
	
	return sqrt(valueX + valueY);
	 
}
char* Scikit::int_to_str(int n)
{
    int basamak=0;
    int temp = n;
    while(temp>0)
    {
        temp = (temp-(temp%10))/10;
        basamak++;
    }
    int numbers[basamak],len = basamak;
    char* c = new char[basamak];
    temp = n;
    while(temp>0)
    {
        numbers[basamak-1] = temp % 10;
        temp = (temp-(temp%10))/10;
        basamak--;
    }
    for(int j=0;j<len;j++)
    {
        c[j] = (char)((int)numbers[j]+48);
    }
    return c;
    
}
bool Scikit::isParallel(Line l1,Line l2)
{
	if(l1.get_Slope() == l2.get_Slope())
		return true;
	return false;
}
bool Scikit::isVertical(Line l1,Line l2)
{
	if(l1.get_Slope()*l2.get_Slope() == -1)
		return true;
	return false;
}

float Scikit::calculate_y_value(float x,float* array,int deg)
{
	float y = 0;
	for(int i=0;i<deg;i++)
	{
		y+=pow(x,i)*array[i];
	}
	return y;
}

double Scikit::slope_point(Point p1,Point p2)
{
	return (p2.getY() - p1.getY())/(p2.getX()-p1.getX());
}

float* Scikit::equation_roots(Equation eq)//simdilik 2.derece köklerini veriyor
{
	float* degrees = eq.getCoefficients();
	float* ret = new float[eq.getDegree()-1];
	if(eq.getDegree() == 3)
	{
		float discriminant = pow(degrees[1],2) - (4*degrees[2]*degrees[0]);//b2-4ac
		ret[0] = (sqrt(discriminant)-degrees[1])/(2*degrees[2]);
		ret[1] = ((-1*sqrt(discriminant))-degrees[1])/(2*degrees[2]);
		return ret; 
	}
	
}
float* Scikit::derive_equation(int deg,float* coeff)
{
	float* derive_coefficient = new float[deg-1];
	for(int j=deg-1;j>=0;j--)
	{
		derive_coefficient[j-1] = coeff[j]*j; 
	}
	return derive_coefficient;
	
}
float* Scikit::derive_integrate(int deg,float* coeff)
{
	float* derive_coefficient = new float[deg+1];
	for(int j=deg;j>0;j--)
	{
		derive_coefficient[j] = coeff[j-1]/j; 
	}
	derive_coefficient[0] = 0;
	return derive_coefficient;
	
}



int main(void)
{
/*	Point p1(5.0,8.0);//first script
	Point p2;
	p2.setPoints(1.0,5.0);
	Scikit sc;
	double deger = sc.distance_two_point(p1,p2);
	cout<<"Distance = "<<deger<<endl;
	cout<<"Egim = "<<sc.slope_point(p1,p2)<<endl;
	*/
	
	/*
	Point p1(2.0,500.0);
	Point p2(4.0,900.0);
	Line ln;
	ln.setLine(p1,p2);//setli noktalardan doðru denklemi çýkarýyor
	char* ret = ln.findLine_equation();//aradaki karakterleri temizleriz tool da
	cout<<ret<<endl;
	*/
	
	
	///calisti///
	Line ln("0.75x+19");
	Point p3=ln.getLine_point1();
	Point p4 =ln.getLine_point2();
	cout<<"x1 = "<<p3.getX()<<" y1 = "<<p3.getY()<<endl;
	cout<<"x2 = "<<p4.getX()<<" y2 = "<<p4.getY()<<endl;
	cout<<"Egim = "<<ln.get_Slope()<<endl;
	
	
	Equation denk;
	denk.setEquation("1x^2-4x^1-5x^0",2);//x10 a kadar buluyor yakýnda çözeriz
	float* degree = denk.getCoefficients();
	for(int i=0;i<3;i++)
	{
		cout<<i<<".kok = "<<degree[i]<<endl;
	}
	Scikit sc;
	float* root = sc.equation_roots(denk);
	cout<<"1.kok: "<<root[0]<<"\t"<<"2.kok: "<<root[1]<<endl; 
	cout<<"x = 2 icin y = "<<denk.calculate_y_value(2)<<endl;
	
	
	
	Equation eq1;
	eq1.setEquation("1x^12+1x^11-1x^10",12);
	float* deg = eq1.getCoefficients();
	float eski = eq1.calculate_y_value(2);
	float* derive_array = sc.derive_equation(eq1.getDegree(),deg);
	cout<<"Turevsiz x = 2 icin y = "<<eski<<endl;
	cout<<"Turevli x = 2 icin y = "<< sc.calculate_y_value(2,derive_array,eq1.getDegree()-1)<<endl;
	
	/*cout<<endl;
	Equation eq1;
	eq1.setEquation("1x^12+1x^11-1x^10",12);
	float* deg = eq1.getCoefficients();
	float eski = eq1.calculate_y_value(2);
	float* integrate_array = sc.derive_integrate(eq1.getDegree(),deg);
	cout<<"Turevsiz x = 2 icin y = "<<eski<<endl;
	cout<<"Turevli x = 2 icin y = "<< sc.calculate_y_value(2,integrate_array,eq1.getDegree()+1)<<endl;//derece +1 olur*/
	
	
	/*
	Scikit sc;
	Line l1("1x+5");
	Line l2("-1x+5");
	cout<<sc.isVertical(l1,l2)<<endl;
	*/
	/*
    float** m1 = new float*[3];//matrix multiple test
	for(int i=0;i<3;i++)
	{
		m1[i] = new float[3];
	}
	 
	float** m2 = new float*[3];
	for(int i=0;i<3;i++)
	{
		m2[i] = new float[3];
	} 
    for(int row=0;row<3;row++)
    {
       for(int col=0;col<3;col++)
       {
       		m1[row][col] = row*3+col;
       		m2[row][col] = row*3+col;
       }
       
    }
    
   float** res = sc.matrix_multiplecation(m1,m2,3,3,3); 
   for(int row=0;row<3;row++)
   {
       for(int col=0;col<3;col++)
       {
            cout<<res[row][col]<<endl;       
       }
       cout<<"\n"<<endl;
   }
	*/
	/*
	float** result = new float*[2];
	for(int i=0;i<2;i++)
	{
		result[i] = new float[1];
	}
	
	float** m1 = new float*[2];
	for(int i=0;i<2;i++)
	{
		m1[i] = new float[1];
	}
	m1[0][0] = 6;
	m1[1][0] = 8;
	result = sc.matrix_transform2D_rotating(m1,45,1);
	cout<<"saat yonu = "<<endl;
	for(int row=0;row<2;row++)
   {
       for(int col=0;col<1;col++)
       {
            cout<<result[row][col]<<endl;       
       }
       cout<<"\n"<<endl;
   }
	result = sc.matrix_transform2D_rotating(m1,45,-1);
	cout<<"saat yonu tersi = "<<endl;
	for(int row=0;row<2;row++)
   {
       for(int col=0;col<1;col++)
       {
            cout<<result[row][col]<<endl;       
       }
      cout<<"\n"<<endl;
   }
  
	//newton raphson test hata payý 0.001  
   float rrres = sc.newton_raphson_square_find(1000);
   cout<<rrres<<endl;
   
    float* corx = new float[10];
    float* cory = new float[10];
    for(int i=0;i<10;i++)
    {
    	corx[i] = i*10+1;
    	cory[i] = i*10+1;
	}
	float corres = sc.correlation_coefficient(corx,cory,10);
	cout<<"correlation = "<<corres<<endl;
	*/


	float** matrix_determinant = new float*[3];//determinant hesabi
	for(int i=0;i<3;i++)
	{
		matrix_determinant[i] = new float[3];
	}

   matrix_determinant[0][0] = 5;
   matrix_determinant[0][1] = 3;
   matrix_determinant[0][2] = 7;
   matrix_determinant[1][0] = 2;
   matrix_determinant[1][1] = 4;
   matrix_determinant[1][2] = 9;
   matrix_determinant[2][0] = 3;
   matrix_determinant[2][1] = 6;
   matrix_determinant[2][2] = 4;
   cout<<sc.matrix_determinant(matrix_determinant,3)<<endl;

	float** mat = new float*[3];//tersi matris hesabý
	for(int i=0;i<3;i++)
	{
		mat[i] = new float[3];
	}

   mat[0][0] = 3;
   mat[0][1] = 2;
   mat[0][2] = 2;
   mat[1][0] = 4;
   mat[1][1] = 3;
   mat[1][2] = 0;
   mat[2][0] = 1;
   mat[2][1] = 0;
   mat[2][2] = 5;
   
   float** matres = new float*[3];//determinant hesabi
	for(int i=0;i<3;i++)
	{
		matres[i] = new float[3];
	}
   matres = sc.matrix_inverse(mat,3);
	
   for(int row=0;row<3;row++)
   {
       for(int col=0;col<3;col++)
       {
            cout<<matres[row][col]<<"\t";       
       }
      cout<<"\n"<<endl;
   }
	
	/*
	//det hesabý
	float** m1 = new float*[2];
	for(int i=0;i<2;i++)
	{
		m1[i] = new float[2];
	}
	
	float** res = new float*[2];
	for(int i=0;i<2;i++)
	{
		res[i] = new float[1];
	}
	
	m1[0][0] = 8.0f;
	m1[0][1] = 6.0f;
	m1[1][0] = 5.0f;
	m1[1][1] = 7.0f;
	
	res[0][0] = 6.0f;
	res[1][0] = 5.0f;
	
	float* sonuc = sc.gramer_rule(m1,res,2);
	for(int i=0;i<2;i++)
	{
		cout<<sonuc[i]<<"\t";
	}
*/
	//ozdeger hesabý jacobi ile
	
	float** m1 = new float*[3];
	for(int i=0;i<3;i++)
	{
		m1[i] = new float[3];
	}	
	
	m1[0][0] = 1.0f;
	m1[0][1] = 0.0f;
	m1[0][2] = 0.0f;
	m1[1][0] = 0.0f;
	m1[1][1] = 2.0f;
	m1[1][2] = 0.0f;
	m1[2][0] = 0.0f;
	m1[2][1] = 0.0f;
	m1[2][2] = 3.0f;
	
	float* sonuc = sc.eigen_values(m1,3);//simetrik matrislerde calýsýyor jacobi yontemi
	for(int i=0;i<3;i++)
	{
		cout<<sonuc[i]<<"\t";
	}
	
	
	
	//turev test 
	Equation eq11;
	eq11.setEquation("2x^1+5x^0",1);
	Equation eq22;
	eq22.setEquation("3x^1+2x^0",1);
	Equation eq33;
	eq33.setEquation("10x^1",1);
	//carpýmýn turevi
	cout<<"Turevsiz x = 2 icin y = "<<eq11.calculate_y_value(2)*eq22.calculate_y_value(2)<<endl;
	cout<<"Turevli x = 2 icin y = "<< sc.multiple_turev(eq11,eq22,2.0f)<<endl;
	//bolumun turevi
	cout<<"Turevsiz x = 2 icin y = "<<eq11.calculate_y_value(2)/eq22.calculate_y_value(2)<<endl;
	cout<<"Turevli x = 2 icin y = "<< sc.divide_turev(eq11,eq22,2.0f)<<endl;
	
	cout<<"Turevli x = 2 icin y = "<< sc.log_turev(eq33,2.0f,10)<<endl;	//suspect
	
	cout<<"Turevsiz x = 19 icin y = "<<cos(eq11.calculate_y_value(19)*PI/180)<<endl;
	cout<<"Turevli x = 19 icin y = "<< sc.cos_turev(eq22,19)<<endl;
	
	system("PAUSE");
	return 0;
}

#include <iostream>
#include<math.h>
#include<string.h>
#include <stdlib.h>
#include <time.h>

#include "scikit.h"
#include "line.h"
#include "point.h"
#include "equation.h"

using namespace std;
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
	Point *p3=ln.getLine_point1();
	Point *p4 =ln.getLine_point2();
	cout<<"x1 = "<<p3->getX()<<" y1 = "<<p3->getY()<<endl;
	cout<<"x2 = "<<p4->getX()<<" y2 = "<<p4->getY()<<endl;
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
    float** m1 = new float*[3];
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
	
	
	
	return 0;
}

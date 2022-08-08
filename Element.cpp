#include "Element.h"
#include <iostream>
#include <cmath>

/// quadratures for integration over a triangular element; exact of polynomials of order up to 5, 
double xgaus[] = {-0.3333333, -0.0597158717,-0.0597158717, -0.8805682564, -0.7974269853, -0.7974269853, 0.5948539707};
double ygaus[] = {-0.3333333, -0.0597158717,-0.8805682564, -0.0597158717, -0.7974269853,  0.5948539707,-0.7974269853};
double wgaus[] = {0.45, 0.2647883055, 0.2647883055, 0.2647883055, 0.251878361, 0.251878361, 0.251878361};

Element::Element(const double x[NO_NODES], const double y[NO_NODES], int i1, int i2, int i3){
	x_lok[0] = x[0];
	x_lok[1] = x[1];
	x_lok[2] = x[2];

	y_lok[0] = y[0];
	y_lok[1] = y[1];
	y_lok[2] = y[2];

	idxs[0] = i1;	
	idxs[1] = i2;	
	idxs[2] = i3;	
}


double Element::phi_i(double ksi1, double ksi2, int i){
	if(i == 0){
		return -0.5 * (ksi1 + ksi2);
	}
	else if(i==1){
		return 0.5 * (1 + ksi1);
	}
	else if(i==2){
		return 0.5 * (ksi2 + 1);
	}
	else{
		std::cout << "Cos nie tak z funkcjami ksztaltu! " << i << ", " << std::endl;
		return -100;
	}
}


double Element::dphi_ksi1(int i){
	double fun=-10;
	switch(i){
		case 0:	
			fun = -0.5;
			break;
		case 1:
			fun = 0.5;
			break;
		case 2:
			fun = 0;
			break;
		default:
			printf("Cod nir tak z indeksami do mapowania!\n");
			break;
	}
	return fun;
}

double Element::dphi_ksi2(int i){
	double fun=-10;
	switch(i){
		case 0:	
			fun = -0.5;
			break;
		case 1:
			fun = 0;
			break;
		case 2:
			fun = 0.5;
			break;
		default:
			printf("Cod nir tak z indeksami do mapowania!\n");
			break;
	}
	return fun;
}

double area(double* x, double *y){
	double x1 = x[0];
	double x2 = x[1];
	double x3 = x[2];

	double y1 = y[0];
	double y2 = y[1];
	double y3 = y[2];

	double S = (-x1*(y3-y2)+x2*y3-x3*y2+(x3-x2)*y1 );
	return 0.5*S;
}

bool Element::point_inside(double x, double y){

	//calkowite pole
	double S = area( x_lok, y_lok );
	// liczymy pola trojkatow zlozonych z 2 wierzcholkow i naszego pktu
	double sumArea = 0;
	double xx[NO_NODES], yy[NO_NODES];
	for (int i = 0; i < NO_NODES; ++i)
	{
		xx[0]=x_lok[0];
		xx[1]=x_lok[1];
		xx[2]=x_lok[2];

		yy[0]=y_lok[0];
		yy[1]=y_lok[1];
		yy[2]=y_lok[2];

		xx[i]=x;
		yy[i]=y;
		sumArea += fabs(area(xx,yy));
	}
	return fabs(sumArea-S)<1e-8;
}


void Element::zetaeta(double x, double y, double &zeta, double &eta){
  double x1 = x_lok[0];
  double x2 = x_lok[1];
  double x3 = x_lok[2];

  double y1 = y_lok[0];
  double y2 = y_lok[1];
  double y3 = y_lok[2];

  eta=-1/(-x1*y3-y1*x2+x2*y3+y1*x3+x1*y2-y2*x3)*(2*x1*y-x1*y2-x1*y3+y1*x3-2*x2*y+x2*y3-y2*x3-2*y1*x+y1*x2+2*y2*x);
  zeta=(2*x1*y-x1*y2-x1*y3+y1*x2+y2*x3+y1*x3-2*y1*x-2*x3*y+2*x*y3-x2*y3)/(-x1*y3-y1*x2+x2*y3+y1*x3+x1*y2-y2*x3);
}


void Element::fill_E(){
	for (int i = 0; i < NO_NODES; ++i){
		double dphi_i_1 = dphi_ksi1(i) * dksi1_dx() + dphi_ksi2(i) * dksi2_dx();
		double dphi_i_2 = dphi_ksi1(i) * dksi1_dy() + dphi_ksi2(i) * dksi2_dy();
		for (int j = 0; j < NO_NODES; ++j){
			E[i][j] = 0;

			double dphi_j_1 = dphi_ksi1(j) * dksi1_dx() + dphi_ksi2(j) * dksi2_dx();
			double dphi_j_2 = dphi_ksi1(j) * dksi1_dy() + dphi_ksi2(j) * dksi2_dy();
			E[i][j] = dphi_i_1 * dphi_j_1 + dphi_i_2 * dphi_j_2;
			E[i][j] *= jm * 2; // 2 comes from the triangle area (because we integrate over a constant independent on ksi1, ksi2)
			E[i][j] *= 0.5; // 1/2 from the schrodingera equation in atomic units
		}
	}
}


void Element::fill_vector(){
	for (int i = 0; i < NO_NODES; ++i){
		F[i] = 0;
	}

}

void Element::fill_V(std::function<double(double x, double y)> potential){
	for (int i = 0; i < NO_NODES; ++i){
		for (int j = 0; j < NO_NODES; ++j){
			V[i][j] = 0;
			for (int k = 0; k < 7; ++k){
				double ksi1 = xgaus[k];
				double ksi2 = ygaus[k];
				double x, y;
				xy(ksi1, ksi2, x, y);
				V[i][j] += jm * phi_i(ksi1, ksi2, i) * phi_i(ksi1, ksi2, j) * potential(x, y) * wgaus[k]; 
			} 
		}
	}

}


void Element::fill_O(){
	for (int i = 0; i < NO_NODES; ++i){
		for (int j = 0; j < NO_NODES; ++j){
			O[i][j] = 0;
			for (int k = 0; k < 7; ++k){
				double ksi1 = xgaus[k];
				double ksi2 = ygaus[k];
				double x, y;
				xy(ksi1, ksi2, x, y);

				O[i][j] += jm * phi_i(ksi1, ksi2, i) * phi_i(ksi1, ksi2, j) * wgaus[k];
			}
		}
	}

}

void Element::calc_jm(){
	jm = 0;
	jm = dx_dksi1() * dy_dksi2() - dy_dksi1() * dx_dksi2();
}

void Element::xy(double ksi1, double ksi2, double &x, double &y){
	x=0;
	y=0;
	for(int l=0; l<3; l++){
		x += x_lok[l] * phi_i(ksi1, ksi2, l);
		y += y_lok[l] * phi_i(ksi1, ksi2, l);
	}
}

double Element::dx_dksi1(){
	double res=0;
	for(int l=0; l<3; l++){
		res += x_lok[l] * dphi_ksi1(l);
	}	
	return res;
}
double Element::dx_dksi2(){
	double res=0;
	for(int l=0; l<3; l++){
		res += x_lok[l] * dphi_ksi2(l);
	}	
	return res;
}
double Element::dy_dksi1(){
	double res=0;
	for(int l=0; l<3; l++){
		res += y_lok[l] * dphi_ksi1(l);
	}
	return res;
}
double Element::dy_dksi2(){
	double res=0;
	for(int l=0; l<3; l++){
		res += y_lok[l] * dphi_ksi2(l);
	}
	return res;
}

double Element::dksi1_dx(){
	return dy_dksi2()/jm;
}
double Element::dksi2_dx(){
	return -dy_dksi1() / jm;
}
double Element::dksi1_dy(){
	return -dx_dksi2() / jm;
}
double Element::dksi2_dy(){
	return dx_dksi1() / jm;
}

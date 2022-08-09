#include "Lead.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <functional> // for invoke
#include </home/alina/Alinkowe/Eigsolvers/eigen-3.4.0/Eigen/Eigenvalues> 

void Lead::add_node(double x, double y, int idx_glob){
	lok.push_back({x,y,idx_glob});
	no_nodes++;
}

// for the nodes of the lead to be ordered from lowest y or x (to make integration easier)
void Lead::sort_nodes(std::function<bool(node &A, node &B)> cmp){
	std::sort(lok.begin(), lok.end(), cmp);
	
	for(int i = 0; i < no_nodes; i++){
		ksi.push_back( std::sqrt(std::pow(lok[0].x-lok[i].x, 2) + std::pow(lok[0].y-lok[i].y, 2)) );
	}
	d = std::sqrt(std::pow(lok[0].x-lok[no_nodes-1].x, 2) + std::pow(lok[0].y-lok[no_nodes-1].y, 2));
}

std::complex<double> Lead::mode_analytical(double xi, int m){
	return std::sqrt(2 / d) * std::sin(m * M_PI / d * xi);
}

// fake fucntion that assumes linear approximation between nodes of the lead
std::complex<double> Lead::mode_fem(double xi, int m){
	// finding in which element the mode is
	int i=1;
	for(int j = 1; j < no_nodes; ++j){
		if (xi > ksi[j-1] && xi < ksi[j]+1e-4){
			i = j;
			break;
		}
	}
	double funsin = (std::sin(m * M_PI / d * ksi[i]) - std::sin(m * M_PI / d * ksi[i-1]))/(ksi[i] - ksi[i-1]) * (xi - ksi[i-1]) + std::sin(m * M_PI / d * ksi[i-1]);
	return std::sqrt(2 / d) * funsin;
}

std::complex<double> Lead::mode_fdm(double xi, int m){
	// finding in which element the mode is
	int i=1;
	int mode_idx = m - 1; // because numbering is from 0, while m is from 1
	for(int j = 1; j < no_nodes; ++j){
		if (xi > xi_num[j-1] && xi < xi_num[j]+1e-4){
			i = j;
			break;
		}
	}
	std::complex<double> funsin = (modes_num[mode_idx][i] - modes_num[mode_idx][i-1])/(xi_num[i] - xi_num[i-1]) * (xi - xi_num[i-1]) + modes_num[mode_idx][i-1];
	return funsin;
}

void Lead::calc_dispersion_relation(std::function<double(double x, double y)> potential){
	int N = no_nodes;	
	Eigen::MatrixXcd A = Eigen::MatrixXcd(N, N);
	double dx = (ksi[no_nodes-1] - ksi[0]) / (N - 1);
	double alpha = 1 / (2 * mass * dx * dx);
	V = std::vector<double>(N, 0);
	// filling the potential energy
	for(int i = 0; i < N; i++){
		//~ double ksi = ksi[0] + dx * i;
		// watch out! it works only when N = no_nodes
		double x = lok[i].x;
		double y = lok[i].y;
		V[i] = potential(x, y);
	}
	std::complex<double> j(0, 1);
	 
	std::ofstream file; 
  	file.open("Ek.txt");
	for(double k =- M_PI/dx; k < M_PI/dx; k += M_PI/dx/100.){
		 //~ filling the A matrix of eigenproblem
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				A(i, j) = 0;
			}
		}
		for(int i = 0; i < N; i++){
			A(i, i) = (4 * alpha + V[i]) - alpha * (std::exp(j * k * dx) + std::exp(-j * k * dx));
		}
		A(0, 0) = 100;
		A(N - 1, N - 1) = 100; // to ensure zeros at the edges
		
		for(int i = 2; i < N - 1; i++){
			A(i - 1, i) = -alpha;
			A(i, i - 1) = -alpha;
		}
		
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces;
		ces.compute(A);
		file << k << "\t";
		for(int i=0; i<N; ++i){
			file << ces.eigenvalues()(i) << "\t";
		}
		file << std::endl;
	}
	file.close();
}
	
void Lead::calculate_modes_fdm(double E, std::function<double(double x, double y)> potential){
	
	//~ calc_dispersion_relation(potential);
	
	ki.clear();
	xi_num.clear();
	modes_num.clear();
	// matrices for the generalized eigenproblem
	int N = no_nodes;	
	//~ int N = 25; //no_nodes;	
	Eigen::MatrixXcd A = Eigen::MatrixXcd(N * 2, N * 2);
	Eigen::MatrixXcd B = Eigen::MatrixXcd(N * 2, N * 2);
	double dx = (ksi[no_nodes-1] - ksi[0]) / (N - 1);
	
	// filing the nodes 
	xi_num.reserve(N);
	for(int i = 0; i < N; i++){
		xi_num.push_back(ksi[0] + dx * i);
	}	
	
	double alpha = 1 / (2 * mass * dx * dx);
	// V vector should have been filled in MES::calc_Hamiltonian, here temporarily zeros
	V = std::vector<double>(N, 0);
	for(int i = 0; i < N; i++){
		// watch out! it works only when N = no_nodes
		double x = lok[i].x;
		double y = lok[i].y;
		V[i] = potential(x, y);
	}
	
	// filling the A nad B matrix of generalized eigenproblem
	for(int i = 0; i < 2 * N; i++){
		for(int j = 0; j < 2 * N; j++){
			A(i, j) = 0;
			B(i, j) = 0;
		}
	}
	for(int i = 0; i < N; i++){
		B(i, i) = 1;
		B(i + N, i + N) = -alpha;
	}
	for(int i = 0; i < N; i++){
		A(i, i + N) = 1;
		A(i + N, i) = alpha;
		A(i + N, i + N) = E - (4 * alpha + V[i]);
	}
	//~ for(int i = 1; i < N; i++){
	for(int i = 2; i < N-1; i++){
		A(i - 1 + N, i + N) = alpha;
		A(i + N, i - 1 + N) = alpha;
	}
	A(N, N) = 100;
	A(2 * N - 1, 2 * N - 1) = 100; // to ensure zeros at the edges
	
	// solving the eigenproblem
	// Eigen does not have a generalized eigensolver for complex matrices, so we multiply by inverted B matrix
	A = B.inverse() * A;
	
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
	ces.compute(A);
	// filing the numerical modes and k's
	for(int i=0; i < 2 * N; ++i){
		if( std::abs(1 - std::abs(ces.eigenvalues()(i))) < 1e-4 && std::arg(ces.eigenvalues()(i)) > 0){
			ki.push_back( std::arg(ces.eigenvalues()(i)) / dx);
			no_modes++;
			
			std::vector<std::complex<double>> wf;
			wf.reserve(N);
			for(int j = 0; j < N; j++){
				wf.push_back(ces.eigenvectors()(j, i) / std::sqrt(dx/2));
			}
			modes_num.push_back(wf);
		}
	}
	no_all_modes = no_modes * 1; // later we can implement the evanescent mdoes for FDM as well
}

	
void Lead::calc_Ni(double E, std::function<double(double x, double y)> potential){
	double dx;
	double xmin;
	double xmax;
	int num = 400; // the denser the better the result for the analytical modes ;) 
	
	double Ewell0 = M_PI*M_PI / 2 / mass / d / d;
	no_modes = 0;
	if(calculation_mode == Lead_mode::ANALYTICAL || calculation_mode == Lead_mode::FEM){
		bool propagating = true;
		no_modes = 0;
		int m = 1;
		
		ki.clear();
		while(propagating){
			double Ei = Ewell0 * m * m;
			if(Ei > E) break;
			double k = std::sqrt( 2 * mass * (E - Ei) ); // in 1/a0
			//~ std::cout << "Ei = " << Ei << ", k= " << k << std::endl;
			ki.push_back(k);
			no_modes++;
			m++;
		}
		no_all_modes = no_modes * 1;
		for(int i = no_modes; i < no_all_modes; ++i){
			double Ei = Ewell0 * m * m;
			double k = std::sqrt( -2 * mass * (E - Ei) ); // in 1/a0
			ki.push_back(k);
			m++;
		}
		
	}
	else if(calculation_mode == Lead_mode::FDM){
		calculate_modes_fdm(E, potential);
	}
	
	Nim.clear();
	for(int m = 1; m <= no_all_modes; ++m){
		std::vector<std::complex<double>> Ni;
		for(int i = 0; i < no_nodes; ++i){
			if(i > 0 && i < no_nodes - 1){
				xmin = ksi[i-1];
				xmax = ksi[i+1];
			}
			else if(i == 0){
				xmin = ksi[i];
				xmax = ksi[i+1];
			}
			else if(i == no_nodes - 1){
				xmin = ksi[i-1];
				xmax = ksi[i];
			}
			dx = (xmax - xmin) / (num - 1);
			// calkowanie ze wzou trapezow 
			std::complex<double> integral = 0;
			for(int j = 1; j<num - 1; j++){
				double xi = xmin + dx * j;
				integral += std::invoke(mode_function, this, xi, m) * phi_i(xi, i); // invoking the mode method we chose at the beginning
				
			}
			integral += 0.5 * std::invoke(mode_function, this, xmin, m) * phi_i(xmin, i);
			integral += 0.5 * std::invoke(mode_function, this, xmax, m) * phi_i(xmax, i);
			Ni.push_back(integral * dx);
		}
		Nim.push_back(Ni);
	}
	
	B = Eigen::MatrixXcd(no_modes, no_modes);
	// overlap integrals
	for(int m = 1; m <= no_modes; ++m){
		for(int n = 1; n <= no_modes; ++n){
			double dx;
			
			std::complex<double> integral = 0;
			for(int i = 1; i < no_nodes; ++i){
				xmin = ksi[i-1];
				xmax = ksi[i];
				int num = 100; 
			
				dx = (xmax - xmin) / (num - 1);
				// integrate using trapeze rule
				for(int j = 1; j<num - 1; j++){
					double xi = xmin + dx * j;
					integral += std::conj(std::invoke(mode_function, this, xi, m)) * std::invoke(mode_function, this, xi, n);
				}
				integral += 0.5 * std::conj(std::invoke(mode_function, this, xmin, m)) * std::invoke(mode_function, this, xmin, n);
				integral += 0.5 * std::conj(std::invoke(mode_function, this, xmax, m)) * std::invoke(mode_function, this, xmax, n);
			}	
			B(m-1, n-1) = integral * dx;
		}
	}
	
	// inverted S
	Sinv = B.inverse();
		
	// test shape functions i calek: zapis do pliku
	//~ std::ofstream file;
	//~ file.open("fcjelead.txt");
	//~ int n=100;
	//~ dx = d/(n-1);
	//~ for (int i = 0; i < n; ++i){
	 //~ double x = dx*i;
	 //~ file << x << "\t" << "\t";
	   
	 //~ for (int j = 0; j < no_nodes; ++j){
		//~ for (int m1 = 1; m1 <= no_modes; ++m1){
			//~ file << std::invoke(mode_function, this, x, m1).real() << "\t"; 
			//~ file << std::invoke(mode_function, this, x, m1).imag() << "\t"; 
		//~ }
		//~ file << Nim[0][j] << "\t"; 
	 //~ }
	 //~ file << std::sin(M_PI / d * x) << std::endl; 
	//~ }
   //~ file.close();
     
}


void Lead::fill_P(int mode_id){
	
	P.clear();
	for(int i = 0; i < no_nodes; i++){
		std::complex<double> p_val = 0;
		p_val = -0.5 * 2 * ki[mode_id] * Nim[mode_id][i] * std::complex<double>(0, 1.);
		//~ p_val = -0.5 * 1 * ki[mode_id] * Nim[mode_id][i] * std::complex<double>(0, 1.) - ; // may be generalized if left and right-propagating modes are different
		P.push_back(p_val);
	}
	
}
	
void Lead::fill_C(){
	
	C.clear();
	for(int i = 0; i < no_nodes; i++){
		std::vector<std::complex<double>> tmp(no_nodes, 0);
		C.push_back(tmp);
	}
	for(int i = 0; i < no_nodes; i++){
		for(int j = 0; j < no_nodes; j++){
			std::complex<double> tmp = 0;
			for(int m = 0; m < no_modes; ++m){
				tmp += -0.5 * ki[m] * std::conj(Nim[m][i]) * Nim[m][j] * std::complex<double>(0, 1.); 
				//~ tmp += -0.5 * ki[m] * Nim[m][i] * Nim[m][j] * std::complex<double>(0, 1.); 
			}
			for(int m = no_modes; m < no_all_modes; ++m){
				tmp += 0.5 * ki[m] * std::conj(Nim[m][i]) * Nim[m][j]; 
				//~ tmp += 0.5 * ki[m] * Nim[m][i] * Nim[m][j]; 
			}
			C[i][j] = tmp;
		}
	}
	
}

double Lead::phi_i(double x1, int i){
	double xmin;
	double xmax;
	double x0;
	double val;
	if(i > 0 && i < no_nodes - 1){
		xmin = ksi[i-1];
		xmax = ksi[i+1];
	}
	else if(i == 0){
		xmin = ksi[i];
		xmax = ksi[i+1];
	}
	else if(i == no_nodes - 1){
		xmin = ksi[i-1];
		xmax = ksi[i];
	}
	else{
		std::cout << "zly indeks w phi leadu!\n";
	}
	
	x0 = ksi[i];
	if(x1 < xmin+1e-6 || x1 > xmax-1e-6){
		val = 0;
	}
	else if(x1 < x0){
		val = (x1 - xmin) / (x0 - xmin);
	}
	else{
		val = (xmax - x1) / (xmax - x0);
	}
	return val;
}

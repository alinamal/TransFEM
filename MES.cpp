#include "MES.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "Msh_reader.h"

 
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::MatrixXcd;
 
Mes::Mes(std::string filename){

    no_leads = 0;
    
    double **vnodes;
    int *vflags;
    int **velems;
    int no_elements;
    int no_nodes;
    
    std::cout << "Reading: " << std::endl;
   
    msh_read(filename, vnodes, vflags, velems, no_nodes, no_elements);
    
    std::cout << "Number of nodes: " << no_nodes << std::endl;
    std::cout << "Number of elements: " << no_elements << std::endl;
    
    Np = no_nodes;
    x_n = new double[Np];
    y_n = new double[Np];
	flags = new int[Np];
    
    // copying the read nodes
    std::cout <<"Saving nodes:" << std::endl;
    for(int i=0;i<no_nodes;i++){
		  x_n[i] = vnodes[i][0] * 100;
		  y_n[i] = vnodes[i][1] * 100; // test wiekszego ukladu		
		  flags[i] = vflags[i];	
    }     
    
    // write the nodes to file
	//~ std::ofstream file;
	//~ file.open("nodes.txt");
	//~ for (int m = 0; m < Np; ++m){
		//~ file << x_n[m] << "\t" << y_n[m] << std::endl;
	//~ }
	//~ file.close();
    
    // copying the read elements 
     for(int i = 0 ; i< no_elements ; i++){ 
		double xs[3];
		double ys[3];
		int idx1 = velems[i][0];
		int idx2 = velems[i][1];
		int idx3 = velems[i][2];

		xs[0] = x_n[idx1];
		xs[1] = x_n[idx2];
		xs[2] = x_n[idx3];

		ys[0] = y_n[idx1];
		ys[1] = y_n[idx2];
		ys[2] = y_n[idx3];

		// check if the nodes are clockwise 
		double S = area( xs, ys);		
		if( fabs(S)>1e-6 ){
			if( S<0 ){ // if not, swap the 2 nodes
				xs[1] = x_n[idx3];
				xs[2] = x_n[idx2];
				ys[1] = y_n[idx3];
				ys[2] = y_n[idx2];
				// and swap the indices
				int tmp = idx2;
				idx2 = idx3;
				idx3 = tmp;
			}
		}
		
		Element tmp(xs, ys, idx1, idx2, idx3);
		tmp.calc_jm();
		elements.push_back(tmp);       
     }
    // dealocation of the data used by the reader
    for(int i=0;i < Np;i++){
        delete [] vnodes[i];
    }
    delete [] vnodes;
    delete [] vflags;
    for(int i=0;i<no_elements;i++){
        delete [] velems[i];
    }
    delete [] velems;
    
	M = no_elements;    

	// define Eigen matrices
	Am = SparseMatrix<std::complex<double>>(Np, Np);
    X = VectorXd(Np), 
    bm = VectorXd(Np);

	// alocating all the matrices  
	S = new std::complex<double>[Np * Np];
	A = new std::complex<double>[Np * Np];
	O = new std::complex<double>[Np * Np];
	V = new std::complex<double>[Np * Np];
	F = new std::complex<double>[Np];
	psi = new std::complex<double>[Np];
  
	// zeroing the matrices
	for (int k = 0; k < Np * Np; ++k) A[k] = 0;
	for (int k = 0; k < Np * Np; ++k) S[k] = 0;
	for (int k = 0; k < Np * Np; ++k) O[k] = 0;
	for (int k = 0; k < Np * Np; ++k) V[k] = 0;
	for (int k = 0; k < Np; ++k) F[k] = 0;
}


Mes::~Mes(){
  delete [] x_n;
  delete [] y_n;
  delete [] A;
  delete [] S;
  delete [] O;
  delete [] V;
  delete [] F;
  delete [] psi;
  delete [] flags;
}


int Mes::S_idx(int i, int j){
  return (Np) * j + i; 
}

double diff1(double x, int i, double (*op)(double, int) ){
  double dx = 0.001;
  return  ( op(x+dx, i) - op(x-dx, i) )/(2.*dx);
}

double diff2(double x, int i, int alfa, double (*op)(double ksi, int i, int alfa) ){
  double dx = 0.001;
  return  ( op(x+dx, i, alfa) + op(x-dx, i, alfa) - op(x, i, alfa)*2. )/(dx*dx);
}



void Mes::solve(double energy, std::vector<int> lead_in_list){
	std::streamsize ssize = std::cout.precision();
	//~ std::ofstream file2;
	//~ file2.open("fcje.txt");

	//~ for(int i = 0; i < no_nodes; ++i){
		//~ if(i > 0 && i < no_nodes - 1){
			//~ xmin = ksi[i-1];
			//~ xmax = ksi[i+1];
		//~ }
		//~ else if(i == 0){
			//~ xmin = ksi[i];
			//~ xmax = ksi[i+1];
		//~ }
		//~ else if(i == no_nodes - 1){
			//~ xmin = ksi[i-1];
			//~ xmax = ksi[i];
		//~ }
		//~ dx = (xmax - xmin) / (num - 1);

		//~ file2 << xmin << "\t" << std::sin(m * M_PI / d * xmin) * phi_i(xmin, i) << "\n";
		//~ for(int j = 1; j<num - 1; j++){
			//~ double xi = xmin + dx * j;
			//~ double funsin1 = (std::sin(m * M_PI / d * ksi[i]) - std::sin(m * M_PI / d * ksi[i-1]))/(ksi[i] - ksi[i-1]) * (ksi[i] - xi) + std::sin(m * M_PI / d * ksi[i-1]);
			//~ file2 << xmin << "\t" << funsin1  << "\n";
		//~ }
		//~ file2 << xmax << "\t" << std::sin(m * M_PI / d * xmax) * phi_i(xmin, i) << "\n";

	//~ }
	//~ file2.close();
   
  
//   // // sprawdzenie mapowania ksi->x
   //~ file2.open("test_glob.txt");
   //~ int n=14;
   //~ double dx = 2./(n-1);
   //~ for (int m = 0; m < M; ++m){
     //~ for( double ksi1=-1; ksi1<1.01; ksi1+=dx ){
       //~ for( double ksi2=-1; ksi2<-ksi1; ksi2+=dx ){
         //~ double x, y;
         //~ elements[m].xy(ksi1, ksi2, x, y);
         //~ file2 << x << "\t" << y << "\t" << ksi1 << "\t" << ksi2 << std::endl; 

       //~ }
       //~ file2 << std::endl; 
     //~ }
     //~ file2 << std::endl; 
     //~ file2 << std::endl; 
   //~ }
   //~ file2.close();
  
	fill_S(energy);

	typedef Eigen::Triplet<std::complex<double>> T;
	std::vector<T> tripletList;
	tripletList.reserve(elements.size() * NO_NODES);

	for(int i = 0; i < Np; ++i){
		for (int j = 0; j < Np; ++j){
			std::complex<double> v_ij = A[S_idx(i, j)];
			//~ if(std::abs(v_ij) > 1e-3) Am.coeffRef(j, i) += v_ij; // this way was slower 
			if(std::abs(v_ij) > 1e-3) tripletList.push_back(T(j, i, v_ij));
		}
	}
	Am.setFromTriplets(tripletList.begin(), tripletList.end());
	Eigen::SparseLU<SparseMatrix<std::complex<double>>, Eigen::COLAMDOrdering<int> > solver;
	// Compute the ordering permutation vector from the structural pattern of A
	solver.analyzePattern(Am); 
	// Compute the numerical factorization 
	solver.factorize(Am); 

	double Ttot = 0;
	double Rtot = 0;
	
	// create the scattering matrices
	SM = std::vector<std::vector<smatrix>>(no_leads, std::vector<smatrix>(no_leads));
	// zero the scattering matrices
	for(int i = 0; i < no_leads; ++i){
		for (int j = 0; j < no_leads; ++j){
			SM[i][j].Tnm = std::vector< std::vector<std::complex<double>> >(leads[i].no_modes, std::vector<std::complex<double>>(leads[j].no_modes, 0));
		}
	}

	//~ for(int i=0; i < lead_in_list.size(); ++i)
	for(int lead_idx: lead_in_list){
		//~ int lead_idx = 1;
		for(int l = 0; l < no_leads; l++){
			if(l != lead_idx){
				for(int idx_in = 0; idx_in < leads[lead_idx].no_modes; ++idx_in){
					
					//~ std::cout << ", no_modes= ";
					//~ for(int li = 0; li < no_leads; li++) std::cout << leads[li].no_modes << ", ";
					//~ std::cout << std::endl;
					
					fill_F(idx_in, lead_idx);
					for(int i = 0; i < Np; ++i){
						bm(i) = F[i];
					}

					//Use the factors to solve the linear system 
					X = solver.solve(bm);

					for (int i = 0; i < Np; ++i){
						psi[i] = X(i);
					}

					// calculate T, R
					//~ std::cout << "T: " << std::endl;
					double Tn = 0;
					for (int m = 1; m <= leads[l].no_modes; ++m){
						std::complex<double> tm =0;
						int nop = 200;
						for(int i = 1; i < leads[l].no_nodes; ++i){
							double dx = (leads[l].ksi[i] - leads[l].ksi[i-1]) / (nop - 1);
							for(int j = 0; j < nop; ++j){
								double x = j * dx + leads[l].ksi[i-1];
								std::complex<double> funpsi = (psi[leads[l].lok[i].idx] - psi[leads[l].lok[i-1].idx])/(leads[l].ksi[i] - leads[l].ksi[i-1]) * (x - leads[l].ksi[i-1]) + psi[leads[l].lok[i-1].idx];
								tm += std::conj(funpsi) * std::invoke(leads[l].mode_function, leads[l], x, m) * dx;
							}
							std::complex<double> funpsi = psi[leads[l].lok[0].idx];
							tm += 0.5 * std::conj(funpsi) * std::invoke(leads[l].mode_function, leads[l], leads[l].ksi[0], m) * dx;
							funpsi = psi[leads[l].lok[leads[l].no_nodes-1].idx];
							tm += 0.5 * std::conj(funpsi) * std::invoke(leads[l].mode_function, leads[l], leads[l].ksi[leads[l].no_nodes-1], m) * dx;
						}
						Tn += std::pow(std::abs(tm), 2) * std::abs(leads[l].ki[m-1] / leads[lead_idx].ki[idx_in]);
						
						SM[lead_idx][l].Tnm[idx_in][m-1] = tm * std::sqrt(std::abs(leads[l].ki[m-1] / leads[lead_idx].ki[idx_in]));
					}
					//~ std::cout << "Tn= " << Tn << std::endl;
					Ttot += Tn;

					//~ std::cout << "R: " << std::endl;
					double Rn = 0;
					for (int m = 1; m <= leads[lead_idx].no_modes; ++m){
						std::complex<double> rm =0;
						int nop = 200;
						for(int i = 1; i < leads[lead_idx].no_nodes; ++i){
							double dx = (leads[lead_idx].ksi[i] - leads[lead_idx].ksi[i-1]) / (nop - 1);
							for(int j = 0; j < nop; ++j){
								double x = j * dx + leads[lead_idx].ksi[i-1];
								// linear approximation of the wave function between the nodes
								std::complex<double> funpsi = (psi[leads[lead_idx].lok[i].idx] - psi[leads[lead_idx].lok[i-1].idx])/(leads[lead_idx].ksi[i] - leads[lead_idx].ksi[i-1]) * (x - leads[lead_idx].ksi[i-1]) + psi[leads[lead_idx].lok[i-1].idx];
								rm += std::conj(funpsi) * std::invoke(leads[lead_idx].mode_function, leads[lead_idx], x, m) * dx;
							}
							std::complex<double> funpsi = psi[leads[lead_idx].lok[0].idx];
							rm += 0.5 * std::conj(funpsi) * std::invoke(leads[lead_idx].mode_function, leads[lead_idx], leads[lead_idx].ksi[0], m) * dx;
							funpsi = psi[leads[lead_idx].lok[leads[lead_idx].no_nodes-1].idx];
							rm += 0.5 * std::conj(funpsi) * std::invoke(leads[lead_idx].mode_function, leads[lead_idx], leads[lead_idx].ksi[leads[lead_idx].no_nodes-1], m) * dx;
						}
						if(m == idx_in+1) rm -= 1;
						Rn += std::pow(std::abs(rm), 2) * std::abs(leads[lead_idx].ki[m-1] / leads[lead_idx].ki[idx_in]);
						
						SM[lead_idx][lead_idx].Tnm[idx_in][m-1] = rm * std::sqrt(std::abs(leads[lead_idx].ki[m-1] / leads[lead_idx].ki[idx_in]));
					}
					//~ std::cout << "Rn= " << Rn << ", Tn+Rn= " << Tn + Rn << std::endl;
					Rtot += Rn;
				}
				std::cout << "Ttot= " << Ttot << std::endl;
				std::cout << "Rtot= " << Rtot << std::endl;
			}
		}
	}
}



void Mes::write_density(std::string filename, double xmin, double xmax, double ymin, double ymax){
	std::ofstream file; 
  	file.open(filename);
	double dx = 4.91; // some value different than the node spacing 
	for(double x = xmin; x<xmax+dx/2; x+=dx ){
	for(double y = ymin; y<ymax+dx/2; y+=dx ){
	  // check if (x,y) is in the mth element
	  for (int m = 0; m < M; ++m){
	    if( elements[m].point_inside(x,y) ){
	      double ze,ee;
	      elements[m].zetaeta(x, y, ze, ee);

	      file << x << "\t" << y << "\t" ;

			std::complex<double> u=0;
			double uabs=0;

			for(int i=0; i<NO_NODES; i++){
			int idx1 = elements[m].idxs[i];
	        
	        //~ u += F[idx1] * elements[m].phi_i(ze, ee, i); 
	        u += psi[idx1] * elements[m].phi_i(ze, ee, i); 
	        uabs += std::abs(psi[idx1]) * elements[m].phi_i(ze, ee, i); 
			}
		  
		  file << u.real() << "\t" << u.imag() << "\t" << uabs;
		  file << std::endl;
	    }
	  }
	}
	file << std::endl;

	}
	file.close();
	
}


void Mes::add_lead(std::function<bool(double x, double y)> lead_position, std::function<bool(node &A, node &B)> cmp){
	leads.push_back(Lead());
	for (int i = 0; i < Np; ++i){
		if(lead_position(x_n[i], y_n[i])){
		  leads[no_leads].add_node(x_n[i], y_n[i], i);
		}
	}
	leads[no_leads].sort_nodes(cmp);
	no_leads++;
}


//~ void Mes::fill_leads(){
  //~ // wypelniamy leady
  //~ // na razie 2 leady
  //~ for (int l = 0; l < no_leads; ++l){ 
	//~ leads.push_back(Lead());
	  //~ for (int i = 0; i < Np; ++i){
		//~ //if( fabs(x_n[i]-xmin)<1e-1 || fabs(x_n[i]-xmax)<1e-1 || fabs(y_n[i]-xmin)<1e-1 || fabs(y_n[i]-xmax)<1e-1 ){
		//~ //if( fabs(y_n[i]-ymin)<1e-2 || fabs(y_n[i]-ymax)<1e-2 ){
		//~ //if(l==0 && fabs(y_n[i]-ymin)<1e-2 && fabs(x_n[i])<ymax+1e-3){
		//~ if(l==0 && fabs(y_n[i]-ymin)<1e-2){
		//~ //if(l==0 && fabs(y_n[i]-ymin)<1e-2 && fabs(x_n[i])<50){
		//~ //if(l==0 && fabs(y_n[i]-ymin)<1e-2 && fabs(x_n[i])<99.9){	
		  //~ leads[l].add_node(x_n[i], y_n[i], i);
		//~ }
		//~ //else if(l==1 && fabs(y_n[i]-ymax)<1e-2 && fabs(x_n[i])<ymax+1e-3){
		//~ else if(l==1 && fabs(y_n[i]-ymax)<1e-2){
		//~ //else if(l==1 && fabs(y_n[i]-ymax)<1e-2 && fabs(x_n[i])<50){
		//~ //else if(l==1 && fabs(y_n[i]-ymax)<1e-2 && fabs(x_n[i])<99.9){
		  //~ leads[l].add_node(x_n[i], y_n[i], i);
		//~ }
	  //~ }
	  //~ //leads[l].sort_nodes();
  //~ }
  
  //~ // proba posortowania wezlow tak zeby orientacja leadu byla w obu taka sama wzgledem obszaru rozpraszania. Ale nie pomoglo :/ 
  //~ leads[0].sort_nodes([](node &A, node &B) -> bool {
			//~ return A.x < B.x;}
			//~ );
  //~ leads[1].sort_nodes([](node &A, node &B) -> bool {
			//~ return A.x > B.x;}
			//~ );
			
  //~ //// wypisanie na probe :
  //~ //for (int l = 0; l < 2; ++l){ 
	//~ //  std::cout << "lead " << l << ": \n";
	 //~ // for (int i = 0; i < leads[l].no_nodes; ++i){
		//~ //  std::cout << i << "= " << leads[l].ksi[i] << ", " << leads[l].lok[i].x << ", " << leads[l].lok[i].idx << "\n";
//~ //	  }	  
//~ //	  std::cout << std::endl;
 //~ // }
  
//~ }

void Mes::save_mesh(){
	std::ofstream file;
	
	file.open("mesh.txt");
	for (int m = 0; m < elements.size(); ++m){
		for (int i = 0; i < NO_NODES; ++i)	{
			file << elements[m].x_lok[i] << "\t" << elements[m].y_lok[i] << "\t" << i << "\t" << elements[m].idxs[i] << std::endl;
		} 
    file << elements[m].x_lok[0] << "\t" << elements[m].y_lok[0] << "\t" << 0 << "\t" << elements[m].idxs[0] << std::endl << std::endl << std::endl;
	}
	std::cout << "saved mesh " << std::endl;
	file.close();
}

void Mes::fill_Hamiltonian(std::function<double(double x, double y)> potential){
	this->potential = potential;

	// //zerowanie macierzy na poczatek
	for (int k = 0; k < Np * Np; ++k) S[k] = 0;
	for (int k = 0; k < Np * Np; ++k) O[k] = 0;
	for (int k = 0; k < Np * Np; ++k) V[k] = 0;
	for (int k = 0; k < Np; ++k) F[k] = 0;

	for (int m = 0; m < M; ++m){
		elements[m].fill_E();
		elements[m].fill_O();
		elements[m].fill_V(potential);
		elements[m].fill_vector();
	}

	// the matrices are assembled here, to be used for the system of equations
	for (int m = 0; m < elements.size(); ++m){
		for (int i = 0; i < NO_NODES; ++i){
			int p = elements[m].idxs[i];
			for (int j = 0; j < NO_NODES; ++j){
				int q = elements[m].idxs[j];
				S[S_idx(p, q)] += elements[m].E[i][j];
				O[S_idx(p, q)] += elements[m].O[i][j];
				V[S_idx(p, q)] += elements[m].V[i][j];
			}
		}
	}
}

void Mes::fill_Hamiltonian(std::function<double(double x, double y)> potential, std::function<bool(double x, double y)> Dirichlet){
	Dirichlet_boundary = Dirichlet;
	//~ this->potential = potential;

	//~ // //zerowanie macierzy na poczatek
	//~ for (int k = 0; k < Np * Np; ++k) S[k] = 0;
	//~ for (int k = 0; k < Np * Np; ++k) O[k] = 0;
	//~ for (int k = 0; k < Np * Np; ++k) V[k] = 0;
	//~ for (int k = 0; k < Np; ++k) F[k] = 0;

	//~ for (int m = 0; m < M; ++m){
		//~ elements[m].fill_E();
		//~ elements[m].fill_O();
		//~ elements[m].fill_V(potential);
		//~ elements[m].fill_vector();
	//~ }

	//~ // the matrices are assembled here, to be used for the system of equations
	//~ for (int m = 0; m < elements.size(); ++m){
		//~ for (int i = 0; i < NO_NODES; ++i){
			//~ int p = elements[m].idxs[i];
			//~ for (int j = 0; j < NO_NODES; ++j){
				//~ int q = elements[m].idxs[j];
				//~ S[S_idx(p, q)] += elements[m].E[i][j];
				//~ O[S_idx(p, q)] += elements[m].O[i][j];
				//~ V[S_idx(p, q)] += elements[m].V[i][j];
			//~ }
		//~ }
	//~ }
	fill_Hamiltonian(potential);
}
	
void Mes::fill_S(double energy){
	for (int k = 0; k < Np * Np; ++k) A[k] = S[k];

	// the incident carriers for transport: filling the leads' matrices C and P (see Lent's paper)
	for(int l = 0; l < no_leads; ++l){
		leads[l].calc_Ni(energy, potential);
		leads[l].fill_C();
	}

	for (int l = 0; l < no_leads; ++l){
		for (int i = 0; i < leads[l].no_nodes; ++i){
			for (int j = 0; j < leads[l].no_nodes; ++j){
				A[S_idx(leads[l].lok[j].idx, leads[l].lok[i].idx)] += leads[l].C[i][j]; 
			}
		}
	}

	for (int i = 0; i < Np; ++i){
		for (int j = 0; j < Np; ++j){
			A[S_idx(i, j)] += (V[S_idx(i, j)] - energy) * O[S_idx(i, j)];
		}  
	}

	// warunki brzegowe: zero na krawedzi, jak do rownania Laplace'a
	for (int i = 0; i < Np; ++i){
		if(Dirichlet_boundary != nullptr){
			if(Dirichlet_boundary(x_n[i], y_n[i])){
				for (int j = 0; j < Np; ++j){
					A[S_idx(j, i)] = 0; // the entire row is set to zero, only on the diagonal we set 1
				}
				A[S_idx(i, i)] = 1;
			}
		}
		else{
			if(flags[i] == flag_d1){
				for (int j = 0; j < Np; ++j){
					A[S_idx(j, i)] = 0; // the entire row is set to zero, only on the diagonal we set 1
				}
				A[S_idx(i, i)] = 1;
			}
		}
	}

}


void Mes::fill_F(int mode_in, int lead_idx){

	for (int k = 0; k < Np; ++k) F[k] = 0;

	leads[lead_idx].fill_P(mode_in);
	

	for (int i = 0; i < leads[lead_idx].no_nodes; ++i){
		F[leads[lead_idx].lok[i].idx] += leads[lead_idx].P[i];  
	}

	// boundary conditions
	for (int i = 0; i < Np; ++i){ 
		if(Dirichlet_boundary != nullptr){
			if(Dirichlet_boundary(x_n[i], y_n[i])){
				F[i] = 0;
			}
		}
		else{
			if(flags[i] == flag_d1){
				F[i] = 0;// + 1 * y_n[i];// * y_n[i] * x_n[i] * x_n[i]/2000; // tu wejdzie wartosc napiecia na brzegu. Chwilowo dajemy po prostu kondensator plaski 
			}
		}
	}
}


void Mes::choose_lead_calc_mode(Lead_mode mod){
	for(int i=0; i < no_leads; ++i){
		leads[i].calculation_mode = mod;
		if(mod == Lead_mode::ANALYTICAL){
			leads[i].mode_function = &Lead::mode_analytical; // obtaining a pointer to member function of Lead
		}
		else if(mod == Lead_mode::FEM){
			leads[i].mode_function = &Lead::mode_fem;
		}
		else if(mod == Lead_mode::FDM){
			leads[i].mode_function = &Lead::mode_fdm;
		}
		else{
			std::cout << "not implemented" << std::endl;
		}
	}
}


double determinant( double x11, double x12, double x13, double x21, double x22, double x23, double x31, double x32, double x33 ){
	return x13*( x21*x32-x31*x22 ) - x23*( x11*x32-x31*x12 ) + x33*( x11*x22-x21*x12 );
}
double dist2( double x1, double y1, double x2, double y2 ){
	return pow(x1-x2,2) + pow(y1-y2,2);
}


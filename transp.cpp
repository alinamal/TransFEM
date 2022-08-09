#include <iostream>
#include <fstream>
#include <cmath>
#include "MES.h"


// How to imrpove? 
// - evanescent modes! :) 
// - detect errors like trying to run solve without choosing the mode function, itd.

void run_2_leads();
void run_3_leads();
void run_4_leads();

int main(){

	run_2_leads();
	//~ run_3_leads();
	//~ run_4_leads();
	
}

void run_2_leads(){
	double ymax = 100; 
	double ymin = -100;
	double xmax = 100; 
	double xmin = -100;

	Mes syst("Meshes/KwadLow20.msh");
	//~ Mes syst("Meshes/KwadLow10.msh"); //this is a denser mesh
	syst.save_mesh();
	syst.add_lead( 
		[&](double x, double y) -> bool{
			return fabs(y-ymin)<1e-2; // the lead in the bottom
		},
		[](node &A, node &B) -> bool {
			return A.x < B.x;
		}
	);
		
	syst.add_lead( 
		[&](double x, double y) -> bool{
			return fabs(y-ymax)<1e-2; // the lead in the top
		},
		[](node &A, node &B) -> bool {
			return A.x > B.x;
		}
	);
	
	syst.choose_lead_calc_mode(Lead_mode::FDM); // ANALYTICAL or FDM are implemented, but FEM not implemented
	double energy = 0.005;
	std::ofstream file; 
  	file.open("T2.txt");
  	syst.fill_Hamiltonian(
  		[&](double x, double y) -> double { 
			// here we define a QPC potential
			return 0.1 * std::exp(-(std::pow((x - xmax)/20,2)+std::pow(y/(20),2)))
				 + 0.1 * std::exp(-(std::pow((x - xmin)/20,2)+std::pow(y/(20),2))); },
			// here we define a finite potential well
			//~ return 0.01 * (std::tanh(-(x + xmax/2)/5) + 
						  //~ std::tanh((x - xmax/2)/5) + 2); },
		[&](double x, double y){return fabs(x-xmin)<1e-2 || fabs(x-xmax)<1e-2;} // the side edges have Dirichlet b.c.
	);
	
	// scan the incident electron energy
  	for(energy = 0.0001; energy < 0.005; energy += 0.0001){
  	//~ for(energy = 0.002; energy < 0.01; energy += 0.0002){
		syst.solve(energy, {0});
		
		double Ttot01 = 0;
		double Rtot00 = 0;
		// we can access the scattering matrix elements of choice
		for (auto& v : syst.SM[0][1].Tnm)
			for (auto& n: v)
				Ttot01 += std::pow(std::abs(n), 2);
		for (auto& v : syst.SM[0][0].Tnm)
			for (auto& n: v)
				Rtot00 += std::pow(std::abs(n), 2);
		file << energy << "\t" << Ttot01 << "\t" << Rtot00 << "\t" << Ttot01 + Rtot00 << std::endl;
		
	}
	file.close();
	
	// saving the wave function to file
	syst.write_density("rozw2.txt", xmin, xmax, ymin, ymax);
}


void run_3_leads(){
	double ymax = 100; 
	double ymin = -100;
	double xmax = 100; 
	double xmin = -100;

	Mes syst("Meshes/KwadLow20.msh");
	//~ Mes syst("Meshes/KwadLow10.msh"); //this is a denser mesh
	//~ syst.save_mesh();
	
	// now we divide it into 2 leads on the left
	syst.add_lead( 
		[&](double x, double y) -> bool{
			return fabs(y-ymin)<1e-2 && x<-1e-2;
		},
		[](node &A, node &B) -> bool {
			return A.x < B.x;
		}
	);
	syst.add_lead( 
		[&](double x, double y) -> bool{
			return fabs(y-ymin)<1e-2 && x>1e-2;
		},
		[](node &A, node &B) -> bool {
			return A.x < B.x;
		}
	);
	
	syst.add_lead( 
		[&](double x, double y) -> bool{
			return fabs(y-ymax)<1e-2;
		},
		[](node &A, node &B) -> bool {
			return A.x > B.x;
		}
	);
	int no_leads = 3;
	
	syst.choose_lead_calc_mode(Lead_mode::ANALYTICAL); // ANALYTICAL or FDM are implemented, but FEM not implemented
	double energy = 0.00045;
	std::ofstream file; 
  	file.open("T.txt");
  	syst.fill_Hamiltonian(
  		[](double x, double y) -> double {	return 0.; },
		[&](double x, double y){return fabs(x-xmin)<1e-2 || fabs(x-xmax)<1e-2 || fabs(x)<20+1e-2 ;}
	);
	
  	//~ for(energy = 0.0001; energy < 0.0025; energy += 0.00004){
  	for(energy = 0.001; energy < 0.01; energy += 0.0002){
		syst.solve(energy, {0, 1, 2});
		
		file << energy << "\t";
		std::vector<std::vector<double> > Tkl = std::vector<std::vector<double> >(no_leads, std::vector<double>(no_leads, 0) );
		for(int i = 0; i < no_leads; i++){
			for(int j = 0; j < no_leads; j++){
				for (auto& v : syst.SM[i][j].Tnm)
					for (auto& n: v)
						Tkl[i][j] += std::pow(std::abs(n), 2);
				file << Tkl[i][j] << "\t";
			}
		}
		file << std::endl;
	}
	file.close();
	
	syst.write_density("rozw3.txt", xmin, xmax, ymin, ymax);
}



void run_4_leads(){
	double ymax = 150; 
	double ymin = -150;
	double xmax = 200; 
	double xmin = -200;

	Mes syst("Meshes/Rectangle_with_leads.msh");
	syst.save_mesh();
	
	// now we divide it into 2 leads on the left
	syst.add_lead( 
		[&](double x, double y) -> bool{
			return fabs(y-ymin)<1e-2;
		},
		[](node &A, node &B) -> bool {
			return A.x < B.x;
		}
	);
	
	syst.add_lead( 
		[&](double x, double y) -> bool{
			return fabs(x-xmin)<1e-2;
		},
		[](node &A, node &B) -> bool {
			return A.y < B.y;
		}
	);
	
	syst.add_lead( 
		[&](double x, double y) -> bool{
			return fabs(y-ymax)<1e-2;
		},
		[](node &A, node &B) -> bool {
			return A.x > B.x;
		}
	);
	
	syst.add_lead( 
		[&](double x, double y) -> bool{
			return fabs(x-xmax)<1e-2;
		},
		[](node &A, node &B) -> bool {
			return A.y < B.y;
		}
	);
	
	int no_leads = 4;
	
	syst.choose_lead_calc_mode(Lead_mode::ANALYTICAL); // ANALYTICAL or FDM are implemented, but FEM not implemented
	double energy = 0.00045;
	std::ofstream file; 
  	file.open("T4.txt");
  	syst.fill_Hamiltonian(
  		[](double x, double y) -> double {	return 0.; },
		[&](double x, double y){
			return ((fabs(y-100)<1e-2 || fabs(y+100)<1e-2) && fabs(x)>100 ) 
				|| ((fabs(x-100)<1e-2 || fabs(x+100)<1e-2) && fabs(y)>100);
		}
	);
	
  	for(energy = 0.001; energy < 0.005; energy += 0.0002){
		syst.solve(energy, {0, 1, 2, 3});
		
		file << energy << "\t";
		std::vector<std::vector<double> > Tkl = std::vector<std::vector<double> >(no_leads, std::vector<double>(no_leads, 0) );
		for(int i = 0; i < no_leads; i++){
			for(int j = 0; j < no_leads; j++){
				for (auto& v : syst.SM[i][j].Tnm)
					for (auto& n: v)
						Tkl[i][j] += std::pow(std::abs(n), 2);
				file << Tkl[i][j] << "\t";
			}
		}
		file << std::endl;
	}
	file.close();
	
	syst.write_density("rozw4.txt", xmin, xmax, ymin, ymax);
}

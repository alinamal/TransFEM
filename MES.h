#pragma once
#include <fstream>
#include "Element.h"
#include "Lead.h"
#include <vector>
#include <complex>
#include <string>
#include <functional> // for std::function
#include "/home/alina/Alinkowe/Eigsolvers/eigen-3.4.0/Eigen/Sparse"
#include "/home/alina/Alinkowe/Eigsolvers/eigen-3.4.0/Eigen/StdVector"

/// Used to store the scattering matrices Tnm between the leads
struct smatrix{
	std::vector<std::vector<std::complex<double>>> Tnm;
};

class Mes{
public:
	std::vector<std::vector<smatrix>> SM; ///< scattering matrices between all the leads
	
	Mes(std::string filename);
	~Mes();

	/// Maps 2 indices of the matrix to 1 index
	int S_idx(int i, int j); 
	
	/** 
	 * Adds lead to the system
	 * 
	 * @param lead_position function describing coordinates of the lead nodes
	 * @param cmp function that tells how to sort the nodes in the lead 
	 */
	void add_lead(std::function<bool(double x, double y)> lead_position, std::function<bool(node &A, node &B)> cmp); 
	
	/** 
	 * Fills all the elements' matrices E, O, V, and global S, O, V
	 * needs to be called when defining the system, changing Hamiltonian parameters, e.g. potential
	 * 
	 * @param potential function describing the potential profile
	 * @param Dirichlet function describing coordinates of the nodes with Dirichlet hard wall b.c.
	 */
	void fill_Hamiltonian(std::function<double(double x, double y)> potential); 
	
	/** 
	 * Variant where nodes with Dirichlet b.c. are defined by hand; can be used when one wants to overwrite the b.c. defined in gmsh.
	 * 
	 * @param potential function describing the potential profile
	 * @param Dirichlet function describing coordinates of the nodes with Dirichlet hard wall b.c.
	 */
	void fill_Hamiltonian(std::function<double(double x, double y)> potential, std::function<bool(double x, double y)> Dirichlet); 
	
	/**
	 * Assembles the system of equations (with the leads), needs to be called whenever energy changes.
	 * Sets the Dirichlet boundary condition (wave function = 0) at nodes 
	 * (i) using user-defined logical function Dirichlet passed to fill_Hamiltonian.
	 * (ii) or defined in gmsh .msh file; this option is used when no logical function Dirichlet is passed to fill_Hamiltonian.
	 * 
	 * @param energy Fermi energy
	 */	
	void fill_S(double energy); 
	
	/**
	 *  Solves the problem, and fills the scattering matrix, 
	 * 
	 * @param energy The incident carrier energy
	 * @param lead_in_list the indices of the leads in which we inject carriers
	 */
	void solve(double energy, std::vector<int> lead_in_list); 
	
	/// Save the nodes and elements for visualisation
	void save_mesh(); 
	
	/// Pick the modes out of all the enum class Lead_mode options
	void choose_lead_calc_mode(Lead_mode mod); 

	/// Save the wave function modulus squared for the latest solution (only no_modes in 1st lead)
	void write_density(std::string filename, double xmin, double xmax, double ymin, double ymax); 
private:
	/// number of all the nodes 
	int Np;
	/// number of elements 
	int M; 
	int no_leads; 
	/// all the elements
	std::vector<Element> elements; 
	/// list of the leads
	std::vector<Lead> leads; 
	/// coordinates of the nodes
	double *x_n; 
	double *y_n;
	int *flags;
	
	std::complex<double> *S; /// global stiffness matrix 
	std::complex<double> *O; /// global overlap matrix
	std::complex<double> *A; /// matrix for the system of equations (to be copied to Eigen sparse matrix Am)
	std::complex<double> *F; /// free terms
	std::complex<double> *V; /// matrix of the integrals with the potential  
	std::complex<double> *psi; /// solution vector
	
    Eigen::SparseMatrix<std::complex<double>> Am; 
    Eigen::VectorXcd bm; /// free variables for Eigen
    Eigen::VectorXcd X; /// solutions from Eigen (copied to psi after solving)
  
	/// describes the entire edge with Dirichlet b.c.
    std::function<bool(double x, double y)> Dirichlet_boundary; 
  
	/// describes the potential profile; used e.g. in Lead::calc_Ni
    std::function<double(double x, double y)> potential; 
    
    /**
     * Assembles the free (source) terms. Must be called for each incident mode (it is done automatically inside fill_S)
     * 
     * @param mode_in index of the mode to which carriers are injected
     * @param lead_idx index of the lead to which carriers are injected
     */
    void fill_F(int mode_in, int lead_idx);	
};

/// Utils
/// Numerical derivative of function op
double diff1(double x, int i, double (*op)(int, double ) );	
/// 3x3 matrix determinant
double determinant( double x11, double x12, double x13, double x21, double x22, double x23, double x31, double x32, double x33); 
/// Distance between points in 2D
double dist2( double x1, double y1, double x2, double y2 ); 

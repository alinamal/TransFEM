#pragma once
#include <vector>
#include <complex>
#include <functional> // for std::function
#include "/home/alina/Alinkowe/Eigsolvers/eigen-3.4.0/Eigen/Dense"

class Mes;

/// choosing the calculation method for modes
enum class Lead_mode{
	ANALYTICAL, ///< analytical for an infinite well; potential inside well is zero
	FDM, ///< finite difference method
	FEM  ///< finite element method (not implemented)
};

/// Used to store the coordinates and indices of the nodes in the lead
struct node{
	double x;
	double y;
	int idx;
};

class Lead{
	friend class Mes;
public:
	/**
	 *  Linear 1D shape functions within the lead 
	 * 
	 * @param x1 coordinate in the lead (between 0 and d -- width of the lead)
	 * @param i index of the node where the shape function is 1 (counted from 0)
	 */
	double phi_i(double x1, int i); 
private:
	int no_nodes = 0;
	/// number of propagating modes
	int no_modes = 0;
	/// number of all modes including both propagating and evanescent. Can equal no_modes for only propagataing modes
	int no_all_modes = 0;
	/// lead's width
	double d;
	/// effective mass
	double mass = 1;
	/// coordinates and global indices of nodes included in the lead
	std::vector<node> lok;
	/// coordinates of the lead's nodes mapped to the interval (0, d)
	std::vector<double> ksi;
	/// overlap matrix of the nodes
	Eigen::MatrixXcd B;
	/// inverse of overlap matrix
	Eigen::MatrixXcd Sinv;
	
	/// integral of psi_i * mode over the lead's elements
	std::vector<std::vector<std::complex<double>>> Nim;
	/// the wave numbers of the modes. First no_modes are propagating, the rest up tp no_all_modes are evanescent
	std::vector<double> ki;
	/// the matrix built up of Nim.T * Nim
	std::vector<std::vector<std::complex<double>>> C;
	/// the sourse terms
	std::vector<std::complex<double>> P;
	/// How the modes are calculated. Default is analytical unless anything else is chosen
	Lead_mode calculation_mode = Lead_mode::ANALYTICAL; 
	
	/// Numerically obtained propagating modes
	std::vector<std::vector<std::complex<double>>> modes_num; 
	/// Nodes for numerical modes
	std::vector<double> xi_num; 
	/// Potential in the lead (just in 1D)
	std::vector<double> V; 

	/**
	 * Adds node to the lead
	 * 
	 * @param x global x coordinate
	 * @param y global y coordinate
	 * @param idx_glob global index of the node
	 */
	void add_node(double x, double y, int idx_glob);
	
	/**
	 * Sort nodes added to the lead. THe global array of nodes in not affected.
	 * 
	 * @param cmp function that accepts two nodes and returns bool value indicating the sorting condition for nodes
	 */
	void sort_nodes(std::function<bool(node &A, node &B)> cmp);
	
    /**
     * Calculates the free terms for lead (matrix P from Lent's paper). 
     * 
     * @param mode_in index of the mode to which carriers are injected
     */
	void fill_P(int mode_id);
	
	/// Calculates the matrix C from Lent's paper
	void fill_C();
	/// Calculates the modes at energy E. For Leadmode::ANALYTICAL is works only  with V=0 in the lead
	void calc_modes(double E); 
	/// Calculates the integrals of phi_i * mode 
	void calc_Ni(double E); 
	
	/// The calculation method for modes
    std::complex<double>(Lead::*mode_function)(double, int);
	std::complex<double> mode_analytical(double xi, int m); ///< analytical for an infinite well; potential inside well is zero
	std::complex<double> mode_fem(double xi, int m); ///< finite element method (not implemented; returns sin, linearly interpolated between nodes
	std::complex<double> mode_fdm(double xi, int m); ///< finite different method
	
	/// Calculates modes with finite difference method
	void calculate_modes_fdm(double E);
	/// Calculates dispersion relation of the lead in the 1D Brillouin zone with finite difference method
	void calc_dispersion_relation();
};


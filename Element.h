#pragma once
#define NO_NODES 3
#include <functional> // for std::function

class Mes;

class Element{
	friend class Mes;
public:
	Element(const double x[NO_NODES], const double y[NO_NODES], int i1, int i2, int i3);
	/**
	 * (linear) shape function
	 * 
	 * @param ksi1 local coordinate 'x'
	 * @param ksi2 local coordinate 'y'
	 * @param i index of the shape function (according to the order of nodes in (x_lok, y_lo)
	 */
	static double phi_i(double ksi1, double ksi2, int i);
	/**
	 * Numerical derivative of the shape functions w.r.t local coordinates. 
	 * Here shape functions are linear so the derivative does not depend on x, y
	 */
	static double dphi_ksi1(int i);
	static double dphi_ksi2(int i);
	/// mapping from local coordinates ksi1, ksi2 to global x, y 
	void xy(double ksi1, double ksi2, double &x, double &y);
	/// Derivatives used for chain rule
	double dx_dksi1();
	double dx_dksi2();
	double dy_dksi1();
	double dy_dksi2();

	double dksi1_dx();
	double dksi2_dx();
	double dksi1_dy();
	double dksi2_dy();

	// Fills the stiffness matrix of the element
	void fill_E();
	// Fills the overlap matrix of the element (overlap integrals of the shape functions)
	void fill_O(); 
	/**
	 *  Fills the matrix of integrals over the potential given by the function provided
	 * 
	 * @param potential function that accepts point coordinates and returns potential 
	 */
	void fill_V(std::function<double(double x, double y)> potential); 
	/// fills the local load vector (here with zeros)
	void fill_vector();
	/// Checks is a given point lies within the element
	bool point_inside(double, double);
	/**
	 * Maps the global point to the local coordinates of the element
	 * 
	 * @param x coordinate x of the point
	 * @param y coordinate y of the point
	 * @param[out] zeta local coordinate 'x' of the point in the element
	 * @param[out] eta local coordinate 'y' of the point in the element
	 */
	void zetaeta(double x, double y, double &zeta, double &eta);
private:
	/// Coordinates of the triangular element nodes
	double x_lok[NO_NODES];
	double y_lok[NO_NODES];
	/// global indices of the nodes
	int idxs[NO_NODES];
	/// local stiffness matrix of the element
	double E[NO_NODES][NO_NODES];
	/// local overlap matrix of the element
	double O[NO_NODES][NO_NODES];
	/// local potential matrix of the element
	double V[NO_NODES][NO_NODES];
	/// local load vector of the element
	double F[NO_NODES];
	/// Jacobian of the mapping from global to local element's coordinates
	double jm;
	/// Calculated the Jacobian
	void calc_jm();
};

/// Calculates the area of a triangle with nodes at coordinated in arrays x, y
double area(double* x, double *y);

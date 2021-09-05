/********************************************************************************
 *     Functions for sxrdcalc							*
 * 										*
 * Calculate theoretical structure factor, patterson map, chi^2...		*
 *										*
 * calc_structure_factor: Calculates structure factors from model structure.	*
 * patterson_map:  Calculates Patterson map and writes it to file.		*
 * calc_chi_square: Calculates chi^2						*
 * set_symmetry_matrix: Sets the symmetry matrix				*
 *										*
 ********************************************************************************/
#include <gsl/gsl_siman.h>



/********************************************************************************************************
 * double calc_structure_factor(const struct atom_structure  atoms[MAX_ATOMS], const struct atom_structure atoms_bulk[MAX_ATOMS] , const struct latt_par_structure lattice_par, double penetration_depth, 
 *           double h, double k, double l, size_t nr_atoms, const double domain_occ[],                 *
 *           const size_t nr_domains, size_t nr_atoms_bulk)                                             *
 * 			  										*
 * 													*
 * Calculates the real and imaginary parts of the structure factor for the (h,k,l) reflection           *
 * from the model structure. Structure factors for different domains are added incoherently.            *
 * Returns the absolute value of the structure factor.							*
 *													*
 * atoms: 												*
 *  struct atom_structure {	 Structure containing the data of atoms in unit cell			*
 *  char name[3]; 		 Abbreviated name for atom (e.g. Si, C) 				*
 *  double x;			 atom position in unit cell x (in fractional coordinates)		*
 *  double y;			 atom position in unit cell y (in fractional coordinates)		*
 *  double z;			 atom position in unit cell z (in fractional coordinates)		*
 *  double scatt_f;		 atomic scattering factor 						*
 *  double dw_par;		 Debye-Waller factor of the atoms  					*
 *  double occ;			 Occupancy of atom							*
 *  } ;													*
 *													*
 * h, k, l: Indices of the reflection                                                            	*
 * nr_atoms: Number of atoms in the structure file	                                                *
 * domain_occ: Occupation of domain                                                                     *
 * nr_domains: number of domains                                                                        *
 * nr_atoms_bulk: Number of atoms in the bulk structure                                                 *
 ********************************************************************************************************/

double calc_structure_factor(struct atom_structure atoms[][MAX_ATOMS], const struct atom_structure atoms_bulk[MAX_ATOMS],const struct latt_par_structure lattice_par, double penetration_depth,  const double h, const double k, const double l, const size_t nr_atoms[], const double domain_occ[], const size_t nr_domains, size_t nr_atoms_bulk);

/********************************************************************************************************
 * void patterson_map(filename[])
 * 
 * Calculate Patterson map and write it to file.
 *
 * Uses calc_patterson to calculate the value of the Patterson map at each position.
 *
 * filename: Name of file  to write the Patterson map to.
 * step_r1: Number of steps in the r1 direction.
 * step_r2: Number of steps in the r2 direction.
 * nr_struc_f: Number of structure factors. 
 *
 * Output format is:
 * 1st line :comment
 * Other lines: r1 r2 value , 
 * where r1 and r2 are the fractional coordinates in the real-space coordinate system corresponding to 
 * the coordinate system of the structure factor data.
 *******************************************************************************************************/
extern void patterson_map(char filename[], unsigned long int step_r1, unsigned long int step_r2, size_t nr_struc_f, char comment[]);


/********************************************************************************************************
 *  extern double calc_electron_diff_map(double r1, double r2, struct structure_factor exp_struc_f[],   *
 * size_t nr_f,  struct atom_structure  atoms[][MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth,  size_t nr_atoms[], size_t nr_domains,         *
 * double domain_occ[], double scale);                                                                  *
 *                             										*
 * Calculates and returns the value of the electron difference map at r1, r2. 				*
 *													*
 *
 * !!!!!!!!!!!! CHECK  SCALE => REALLY IN UNITS OF e????????????????????
 *													
 * struct_f: Matrix of structures containing the structure factors for h, k, l				*
 * 	.h												*
 *	.k												*
 *	.l												*
 *	.value												*
 *	.sigma												*
 * r1, r2: r1 and r2 positions of the electron difference map (reduced coordinates)			*
 * nr_f: Number of structure factors									*
 * atoms: 												*
 *  struct atom_structure {	 Structure containing the data of atoms in unit cell			*
 *  char name[3]; 		 Abbreviated name for atom (e.g. Si, C) 				*
 *  double x;			 atom position in unit cell x (in fractional coordinates)		*
 *  double y;			 atom position in unit cell y (in fractional coordinates)		*
 *  double z;			 atom position in unit cell z (in fractional coordinates)		*
 *  double scatt_f;		 atomic scattering factor 						*
 *  double dw_par;		 Debye-Waller factor of the atoms  					*
 *  } ;													*
 * nr_atoms: Number of atoms in the structure file							*
 * scale: Scale factor of experimental and calculated structure factors					* 
 ********************************************************************************************************/                             										
/* extern double calc_electron_diff_map(double r1, double r2, struct structure_factor exp_struc_f[], unsigned int nr_f, const struct atom_structure  atoms[MAX_ATOMS], unsigned int nr_atoms, double scale); */
 extern double calc_electron_diff_map(double r1, double r2, struct structure_factor exp_struc_f[], size_t nr_f,  struct atom_structure  atoms[][MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth, size_t nr_atoms[], size_t nr_domains, double domain_occ[], double scale); 






/********************************************************************************************************
 * calc_chi_square(struc structure_factor struct_f_exp, struc structure_factor struc_f_calc, _          *
 * _ unsigned int nr_struc_f, unsigned int nr_par)                                                                  	*
 *												        *
 * Calculates chi^2 to check the agreement of experimental and calculated structure factors.   		*
 * A scale factor is varied for mimimum chi^2.                                                          *
 * Assumes that h,k,l values for the same position in the matrix are same.		                *
 *													*
 * struct_f_exp and struc_f_calc: Matrices with structures containing the structure factors for h, k, l	*
 * 	.h												*
 *	.k												*
 *	.l												*
 *	.value												*
 *	.sigma	(only needed for struct_f_exp )								*
 * 											      		*
 * nr_struc_f: number of structure factors.                                                             *
 * nr_par: Number of parameters								  		*
 * initial_scale: The initial scale factor of the two sets  of structure factors.                       *
 ********************************************************************************************************/

/* extern double calc_chi_square(unsigned int nr_struc_f, unsigned int nr_par, double initial_scale); */


/********************************************************************************
 * void set_symmetry_matrix(enum SURFACE_SYMMETRY symmetry_to_set)     		*
 *                                                                     		*
 * Set symmetry matrices                                               		*
 * This sets the following global variables:					*
 * - the global symmetry matrix symmetry_matrix,                   		*
 * - the number of matrices in the symmetry nr_symmetry_matrices.              	*
 *                                                                     		*
 * symmetry_to_set gives the symmetry the symmetry matrix should be set to.	*
 * It can be at present:				               		*
 * - p6 or									*
 * - p6mm									*
 ********************************************************************************/

extern void set_symmetry_matrix(enum SURFACE_SYMMETRY symmetry_to_set);


/****************************************************************************************************************
 * void do_ls_fit(struct atom_structure atoms[][MAX_ATOMS], , struct atom_structure atoms_bulk[MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth,  struct displacements displ[][MAX_ATOMS], int fixed_par[], double start_par[], double domain_occ[], size_t nr_atoms[], size_t nr_domains, size_t nr_atoms_bulk, struct structure_factor exp_struc_f[], size_t nr_struc_f, unsigned int max_iteration, double delta_abs, double delta_rel, struct atom_structure end_positions[][MAX_ATOMS],  unsigned long int rng_seed, enum output print_intermediate, int n_fit, double par_var_width);
 *														*
 * Main function for Levenberg-Marquardt least-squares structure refinement. 					*
 * The structure is refined by fitting to the absolute square of the experimental structure factors.		*
 * Uses the Nonlinear fitting routines from the GNU scientific library						*
 *														*
 * atoms[domain nr][atom nr]: Position of the atom (atom nr) in the unit cell in the domain (domain nr).	*
 * lattice_par: 
 *   struct latt_par_structure {          Lattice parameters a,b,c and                          *
 *   double a;                             angles alpha, beta, gamma                            *
 *   double b;
 *   double c;
 *   double alpha;
 *   double beta;
 *   double gamma;
 *   };
 * penetration_depth: Penetration depth of X-rays into crystal (reciprocal units of lattice parameter)  *
 * displ[domain nr][atom nr]: Displacements for atom (atom nr) in domain (domain nr).				*
 * fixed_parameters[]: Contains the numbers of the fixed parameters						*
 * start_par[nr_param]: Starting values for the fit parameters.							*
 * nr_atoms: Number of atoms in the unit cell.									*
 * nr_atoms_bulk: Number of atoms in the bulk structure                                                 *
 * exp_struc_f[f nr]: Structure factor h,k,l of structure factor (f nr).					*
 * nr_struc_f: Number of structure factors.									*
 * max_iteration: Maximum number of iterations.									*
 * delta_abs: Absolute convergence criterion 									*
 * delta_rel: Relative convergence criterion (not very useful)							*
 * *end_positions: Pointer to array where final positions should be stored.					*
 ****************************************************************************************************************/

extern void do_ls_fit(struct atom_structure atoms[][MAX_ATOMS], struct atom_structure atoms_bulk[], const struct latt_par_structure lattice_par, double penetration_depth,  struct displacements displ[][MAX_ATOMS], int fixed_par[], double start_par[], double domain_occ[], size_t nr_atoms[], size_t nr_domains, size_t nr_atoms_bulk, struct structure_factor exp_struc_f[], size_t nr_struc_f, unsigned long int max_iteration, double delta_abs, double delta_rel, struct atom_structure end_positions[][MAX_ATOMS],  unsigned long int rng_seed, enum output print_intermediate, long unsigned int n_fit, double par_var_width);


/****************************************************************************************************************
 * void do_sa_fit(struct atom_structure atoms[][],  struct atom_structure atoms_bulk[MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth, struct displacements displ[][], double start_par[],		*
 * double domain_occ[], size_t nr_atoms[], size_t nr_domains, size_t nr_atoms_bulk, struct structure_factor exp_struc_f[],	*
 * size_t nr_struc_f, gsl_siman_params_t sa_par, struct atom_structure (*end_positions)[MAX_ATOMS],	*
 *  unsigned long int rng_seed, enum output  print_intermediate) 								*
 *														*
 * Main function for simulated annealing least-squares structure refinement. 					*
 * The structure is refined by fitting to the absolute square of the experimental structure factors.		*
 * Uses the simulated annealing routines from the GNU scientific library					*
 *														*
 * atoms[domain nr][atom nr]: Position of the atom (atom nr) in the unit cell in the domain (domain nr).	*
 * lattice_par: 
 *   struct latt_par_structure {          Lattice parameters a,b,c and                          *
 *   double a;                             angles alpha, beta, gamma                            *
 *   double b;
 *   double c;
 *   double alpha;
 *   double beta;
 *   double gamma;
 *   };
 * penetration_depth: Penetration depth of X-rays into crystal (reciprocal units of lattice parameter)  * * displ[domain nr][atom nr]: Displacements for atom (atom nr) in domain (domain nr).				*
 * fixed_parameters[]: Contains the numbers of the fixed parameters						*
 * start_par[nr_param]: Starting values for the fit parameters. 						*
 * nr_atoms: Number of atoms in the unit cell.									*
 * nr_atoms_bulk: Number of atoms in the bulk structure                                                 *
 * exp_struc_f[f nr]: Structure factor h,k,l of structure factor (f nr).					*
 * nr_struc_f: Number of structure factors.									*
 * max_iteration: Maximum number of iterations.									*
 * delta_abs: Absolute convergence criterion 									*
 * delta_rel: Relative convergence criterion (not very useful)							*
 * *end_positions: Pointer to array where final positions should be stored.					*
 * rng_seed: Seed for random number generator.
 * print_intermediate: Should intermediate information be printed?
 ****************************************************************************************************************/

extern void do_sa_fit(struct atom_structure atoms[][MAX_ATOMS], struct atom_structure atoms_bulk[MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth, struct displacements displ[][MAX_ATOMS], int fixed_par[], double start_par[], double domain_occ[], size_t nr_atoms[], size_t nr_atoms_bulk, size_t nr_domains, struct structure_factor exp_struc_f[], size_t nr_struc_f, gsl_siman_params_t sa_par, struct atom_structure (*end_positions)[MAX_ATOMS], unsigned long int rng_seed, enum output print_intermediate, unsigned long int n_fit, double par_var_width);



/************************************************************************************************
 * double find_scale(const struct structure_factor exp_struc_f[], const size_t nr_f,  struct atom_structure  atoms[MAX_DOMAINS][MAX_ATOMS], const struct atom_structure atoms_bulk[MAX_ATOMS], const size_t nr_atoms[], const double domain_occ[], const size_t nr_domains, size_t nr_atoms_bulk) 
 * 												*
 * Find the scale factor of theoretical structure factors relative to experimental ones.	* 
* this is linear least-squares minimization.							*
 *												*
 * struct_f: Matrix of structures containing the structure factors for h, k, l			*
 * 	.h											*
 *	.k											*
 *	.l											*
 *	.value											*
 *	.sigma											*
 * nr_f: Number of structure factors								*
 * atoms: 											*
 *  struct atom_structure {	 Structure containing the data of atoms in unit cell		*
 *  char name[3]; 		 Abbreviated name for atom (e.g. Si, C) 			*
 *  double x;			 atom position in unit cell x (in fractional coordinates)	*
 *  double y;			 atom position in unit cell y (in fractional coordinates)	*
 *  double z;			 atom position in unit cell z (in fractional coordinates)	*
 *  double scatt_f;		 atomic scattering factor 					*
 *  double dw_par;		 Debye-Waller factor of the atoms  				*
 *  } ;												*
 * nr_atoms: Number of atoms in the structure file.						*
 * nr_domains: Number of domains.								*
 ************************************************************************************************/
extern double find_scale(const struct structure_factor exp_struc_f[], const size_t nr_f,  struct atom_structure  atoms[MAX_DOMAINS][MAX_ATOMS], const struct atom_structure atoms_bulk[MAX_ATOMS],const struct latt_par_structure lattice_par, double penetration_depth,  const size_t nr_atoms[], const double domain_occ[], const size_t nr_domains, size_t nr_atoms_bulk);




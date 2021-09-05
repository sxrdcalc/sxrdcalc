/********************************************************************************
 *     Functions for sxrdcalc							*
 * 										*
 * Calculate theoretical structure factor, patterson map, chi^2...		*
 * Least-squares and simulated annealing fit                                    *
 *										*
 * calc_structure_factor: Calculates structure factors from model structure.	*
 * patterson_map:  Calculates Patterson map and writes it to file.		*
 * calc_chi_square: Calculates chi^2						*
 * set_symmetry_matrix: Sets the symmetry matrix				*
 *										*
 * Internal:									*
 * calc_patterson: Calculates value of the Patterson map for position r1,r2.	*
 ********************************************************************************/

#include "sxrddefs.h"
#include "calcs.h"
/* GNU scientific library */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_siman.h>

/****** Common data structures ******/
/* Structure containing  atom data, displacement vectors and other parameters for structure factor calculation in fit */
static struct fit_all {					
    struct atom_structure (*atoms_ptr)[MAX_ATOMS];	/* Pointer to atom data */
    struct atom_structure (*atoms_bulk_ptr); /* Pointer to atoms of bulk */
    struct displacements (*displ_ptr)[MAX_ATOMS];	/* Pointer to displacement vectors */
    struct structure_factor (*struc_f_ptr);		/* Pointer to structure factor data */
    struct latt_par_structure lattice_par;	/* Structure containing lattice parameters */
    double (*bulk_structure_factor_real_ptr);   /* Pointer to array with real part of bulk structurefactors */
    double (*bulk_structure_factor_im_ptr);   /* Pointer to array with imaginary part of bulk structurefactors */
    double penetration_depth;           /* Penetration depth of X-rays into crystal [reciprocal units of lattice vectors] */
    double (*start_par_ptr);		/* Pointer to starting parameters */
    int (*fixed_par_ptr);		/* Pointer to array with flag for fixed parameters */
    double (*domain_occ_ptr);		/* Pointer to Occupancies of domains */
    size_t (*displtox);		/* Pointer to array for translating from x to displacement numbers */
    gsl_vector (*fit_par_ptr);		/* Pointer to parameters of the fit (only for simulated annealing, not used in least-squares fit */
    size_t (*nr_atoms_ptr);
    size_t nr_atoms_bulk;
    size_t nr_domains;
    size_t nr_struc_f;
    size_t nr_parameters; 
};

/********************************************************
 * Constants 						*
 ********************************************************/

const struct cromer_mann_par cm_par[21][4] = { 			/* Cromer-Mann parameters */
    { /* Nothing 0*/
	{0.0, 0.0 , 0.0},	/* a1, b1, c */
	{0.0, 0.0, 0.0},	/* a2, b2, 0.0 */
	{0.0, 0.0, 0.0},	/* a3, b3, 0.0 */
	{0.0, 0.0 , 0.0}	/* a4, b4, 0.0 */
    },
    { /* H 1 */
	{0.0, 0.0 , 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0 , 0.0}
    },
    { /* He 2 */
	{0.0, 0.0 , 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0 , 0.0}
    },
    { /* Li 3 */
	{0.0, 0.0 , 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0 , 0.0}
    },
    { /* Be 4 */
	{0.0, 0.0 , 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0 , 0.0}
    },
    { /* B 5 */
	{0.0, 0.0 , 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0 , 0.0}
    },
    { /* C 6 */
	{2.310, 20.844, 0.216},		/* a1 b1 c */
	{1.020, 10.208, 0.0},		/* a2 b2  */
	{1.589, 0.569,  0.0},		/* a3 b3  */
	{0.865, 51.651, 0.0}		/* a4 b4  */
    },
    { /* N 7 */
	{12.2126, 0.005700, -11.529},
	{3.13220, 9.89330, 0.0},
	{2.01250, 28.9975, 0.0},
	{1.16630, 0.582600 , 0.0}
    },
    { /* O 8 !!! CHECK !!! */
	{3.04856, 13.2771, 0.2508},
	{2.28680, 5.70110, 0.0},
	{1.5463, 0.3239, 0.0},
	{0.867, 32.9089, 0.0}
    },
/*    { /\* F *\/
	{0.0, 0.0 , 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0 , 0.0}
    },
    { /\* Ne *\/
	{0.0, 0.0 , 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0 , 0.0}
    },
    { /\* Na *\/
	{0.0, 0.0 , 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0 , 0.0}
    },
    { /\* Mg *\/
	{0.0, 0.0 , 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0 , 0.0}
    },
    { /\* Al *\/
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0}
    },							*/
/*     { /\* Cromer-Mann parameters for silicon *\/ */
/* 	{6.292, 2.439, 1.141},					/\* a1 b1 c *\/ */
/* 	{3.035, 32.334, 0.0},					/\* a2 b2  *\/ */
/* 	{1.989, 0.678, 0.0},					/\* a3 b3  *\/ */
/* 	{1.541, 81.694, 0.0}					/\* a4 b4  *\/ */
/*     } */
    { /* Cromer-Mann parameters for silicon 9 */
	{6.2915, 2.4386, 1.1407},					/* a1 b1 c */
	{3.0353, 32.3337, 0.0},					/* a2 b2  */
	{1.9891, 0.6785, 0.0},					/* a3 b3  */
	{1.541, 81.694, 0.0}					/* a4 b4  */
    },
    { /* Si val 10 */
	{5.66269, 2.66520, 1.24707},
	{3.07164, 38.6634, 0.0},
	{2.62446, 0.916946, 0.0},
	{1.39320, 93.5458, 0.0}
    },
    { /* Si4+ CHECK 11 */
	{4.43918, 1.64167, 0.746297},
	{3.20345, 3.43757, 0.0},
	{1.19453, 0.214900, 0.0},
	{0.416530, 6.65365, 0.0}
    },
    {  /*  O- CHECK! 12 */
	{4.19160, 12.8573, 21.9412},
	{1.63969, 4.17236, 0.0},
	{1.52673, 47.0179, 0.0},
	{-20.307, -0.01404, 0.0}
    },
    {  /*  P CHECK! 13  */
	{6.43450, 1.90670, 1.11490},
	{4.17910, 27.1570, 0.0},
	{1.78000,  0.5260, 0.0},
	{1.49080, 68.1645, 0.0}
    },
    {  /*  In CHECK! 14 */
	{19.1624, 0.54760, 4.93910},
	{18.5596, 6.37760, 0.0},
	{4.29480, 25.8499, 0.0},
	{2.03960, 92.8029, 0.0}
    },
    {  /*  In3+ CHECK! 15  */
	{19.1045, 0.551522, 4.99635},
	{18.1108, 6.32470, 0.0},
	{3.78897, 17.3595, 0.0},
	{    0.0,     0.0, 0.0}
    },
    {  /* Au 16 */
      {16.8819, 0.4611, 12.0658},
      {18.5913, 8.6216, 0.0},
      {25.5582, 1.4826, 0.0},
      {5.86, 36.3956, 0.0}
    },

    {  /* Mn 17 */
      {11.2819, 5.3409,1.0896},
      {7.3573, 0.3432, 0.0},
      {3.0193, 17.8674, 0.0},
      {2.2441, 83.7543, 0.0}
    }, 

    {   /* Fe 18 */
      {11.7695, 4.76110, 1.03690},
      {7.35730, 0.30720, 0.0},
      {3.52220, 15.3535, 0.0},
      {2.30450, 76.8805, 0.0}
    },

    {   /* Pb 19 */
      {31.0617, 0.690200, 13.4118},
      {13.0637, 2.35760, 0.0},
      {18.4420, 8.61800, 0.0},
      {5.96960, 47.2579, 0.0}
    }
    

};


/* Function prototypes */
double atom_scatt_f(int element_nr, double scattering_length_sq);

/****************************************************************************************
 *  double calc_scattering_length_sq(const struct latt_par_structure lattice_par, 	*
 * const double h, const double k, const double l)				       	*
 *											*
 * Calculates and returns the square of the scattering length for h, k, l reflection	*
 * (used for calculating the atomic scattering factor)				       	*
 * This is (sin(theta)/lambda)^2 = q^2 / (16*pi^2) 					*
 *											*
 * !!!!! ONLY FOR HEXAGONAL OR CUBIC SYSTEM !!!!					*
 * 											*
 * lattice_par:   Structure with the lattice parameters	a,b,c,alpha,beta,gamma          *
 * h, k, l:	  Miller indizes of reflection						*
 ****************************************************************************************/

double calc_scattering_length_sq(const struct latt_par_structure lattice_param, const double h, const double k, const double l)
{
    double scattering_length_sq;
    double a, b, c, alpha, beta, gamma;
    double sina,sinb,sing,cosa,cosb,cosg;
    
    /* Cubic structure */
    if (lattice_param.alpha == 90.0 && lattice_param.beta == 90.0 && lattice_param.gamma == 90.0 && lattice_param.a == lattice_param.b)	{
      scattering_length_sq = (h*h + k*k) / (lattice_param.a*lattice_param.a)/ 4.0 + l*l/ (lattice_param.c*lattice_param.c)/ 4.0;
    }
    else { 
        /* Hexagonal structure */ 
        if (lattice_param.alpha == 90.0 && lattice_param.beta == 90.0 && lattice_param.gamma == 120.0 && lattice_param.a == lattice_param.b) {
	    scattering_length_sq = (4.0/3.0/(lattice_param.a*lattice_param.a) * (h*h + k*k + h*k) + l*l / (lattice_param.c*lattice_param.c)) / 4.0;
	}
                                  
	/* General structure */
	else { 
   	    a = lattice_param.a;
	    b = lattice_param.b;
	    c = lattice_param.c;
	    alpha = lattice_param.alpha;
	    beta = lattice_param.beta;
	    gamma = lattice_param.gamma;
	    sina = sin(alpha/360.0*2.0*PI);
	    sinb = sin(beta/360.0*2.0*PI);
	    sing = sin(gamma/360.0*2.0*PI);
	    cosa = cos(alpha/360.0*2.0*PI);
	    cosb = cos(beta/360.0*2.0*PI);
	    cosg = cos(gamma/360.0*2.0*PI);

	    scattering_length_sq = 1.0/4.0/(1.0 + 2.0*cosa*cosb*cosg - cosa*cosa - cosb*cosb - cosg*cosg) * (h*h*sina*sina/(a*a) + k*k*sinb*sinb/(b*b) + l*l*sing*sing/(c*c) + 2.0*h*k/(a*b)*(cosa*cosb-cosg) + 2.0*k*l/(b*c)*(cosb*cosg-cosa) + 2.0*l*h/(a*c)*(cosg*cosa - cosb));

	};
    };
    return scattering_length_sq;
}


/********************************************************************************************************
 * void calc_structure_factor_domain(const struct atom_structure  atoms[MAX_ATOMS],const struct latt_par_structure lattice_par, double h, double k,	*
 *                               double l, size_t n_atoms_f, double *f_real, double *f_im)              *
 * 			  										*
 * 													*
 * Calculates the real and imaginary parts of the structure factor for the (h,k,l) reflection           *
 * from the model structure for one domain		                                                        *
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
 * lattice_par: 
 *   struct latt_par_structure {          Lattice parameters a,b,c and                          *
 *   double a;                             angles alpha, beta, gamma                            *
 *   double b;
 *   double c;
 *   double alpha;
 *   double beta;
 *   double gamma;
 *   };
 *													*
 * h, k, l: Indices of the reflection                                                            	*
 * n_atoms_f: Number of atoms in the structure file	                                                *
 * struc_f_real_ptr, *struc_f_im_ptr: pointer to variables where real and imaginary part of             *
 *                                    structure factor should be stored.				*
 ********************************************************************************************************/

void calc_structure_factor_domain(const struct atom_structure atoms[MAX_ATOMS],const struct latt_par_structure lattice_par, const double h, const double k, const double l, const size_t n_atoms_f, double *struc_f_real_ptr, double *struc_f_im_ptr)  
{ 
    double struc_f_real = 0.0;   		/* Real part of the structure factor */
    double struc_f_im = 0.0;   			/* Complex part of the structure factor */
    size_t atom_counter;  		/* Atom number used in loop */
    double dw_factor;   		        /* debye-waller factor: exp(-1*dw_par*Q^2) */
    double aniso_dw;   		        /* anisotropic debye-waller factor: */
    double atom_scattering_f = 0.0;		/* Atomic scattering factor */
    double scattering_length_sq;		/* The length of the scattering vector (used for
						   calculating the atomic scattering factor) */
  
    /* Square of the scattering length  */
    scattering_length_sq = calc_scattering_length_sq(lattice_par, h, k, l);

    /* Calculate the structure factor by summing the terms. */
    for (atom_counter = 0; atom_counter < n_atoms_f; ++atom_counter) {
	/* Add the contribution of atom j to the structure factor. */
	atom_scattering_f = atom_scatt_f(atoms[atom_counter].el_nr, scattering_length_sq); /* Atomic scattering factor */
	dw_factor = exp(-1.0 * atoms[atom_counter].dw_par * scattering_length_sq);	/* Debye-Waller factor  */
	aniso_dw = exp(-1.0 * (atoms[atom_counter].aniso_dw[0]*h*h + atoms[atom_counter].aniso_dw[1]*k*k + atoms[atom_counter].aniso_dw[2]*l*l + 2.0*atoms[atom_counter].aniso_dw[3]*h*k + 2.0*atoms[atom_counter].aniso_dw[4]*h*l + 2.0*atoms[atom_counter].aniso_dw[5]*k*l));
	struc_f_real += (atoms[atom_counter].occ * atom_scattering_f * dw_factor * aniso_dw *cos(2.0*PI*(h*atoms[atom_counter].x+k*atoms[atom_counter].y+l*atoms[atom_counter].z)));
	struc_f_im +=  (atoms[atom_counter].occ * atom_scattering_f * dw_factor * aniso_dw * sin(2.0*PI*(h*atoms[atom_counter].x+k*atoms[atom_counter].y+l*atoms[atom_counter].z)));
    };

    /* Assign structure factor to return */
    *struc_f_real_ptr = struc_f_real;
    *struc_f_im_ptr =struc_f_im ;
  
    return;
}


/************************************************************************************************
 * void calc_structure_factor_bulk(const struct atom_structure atoms_bulk[MAX_ATOMS], size_t  nr_atoms_bulk, const struct latt_par_structure lattice_par, double penetration_depth, double h, double k, double l, double *struc_f_real_ptr, double *struc_f_im_ptr)
 *                                                                                              *
 * Calculate the structure factor of the truncated bulk crystal                                 *
 *                                                                                              *
 * F is multiplied with 1/(1 - exp(-2*PI*i*l) * exp(lattice_par.c/penetration_depth))           *
 * (The formula used is more complex, because it is divided into real and imaginary part        *
 *                                                                                              *
 * atoms_bulk: 										        *
 *  struct atom_structure {	 Structure containing the data of atoms in unit cell	        *
 *  char name[3]; 		 Abbreviated name for atom (e.g. Si, C) 			*
 *  double x;			 atom position in unit cell x (in fractional coordinates)	*
 *  double y;			 atom position in unit cell y (in fractional coordinates)	*
 *  double z;			 atom position in unit cell z (in fractional coordinates)	*
 *  double scatt_f;		 atomic scattering factor 					*
 *  double dw_par;		 Debye-Waller factor of the atoms  				*
 *  double occ;			 Occupancy of atom						*
 *  } ;												*
 *												*
 * h, k, l: Indices of the reflection                                                           *
 * nr_atoms: Number of atoms in the structure file	                                        *
 * lattice_par: 
 *   struct latt_par_structure {          Lattice parameters a,b,c and                          *
 *   double a;                             angles alpha, beta, gamma                            *
 *   double b;
 *   double c;
 *   double alpha;
 *   double beta;
 *   double gamma;
 *   };
 * penetration_depth: Penetration depth of X-rays into crystal (reciprocal units of lattice parameter)
 * struc_f_real_ptr, *struc_f_im_ptr: pointer to variables where real and imaginary part of     *
 *                                    structure factor should be stored.			*
 ************************************************************************************************/

void calc_structure_factor_bulk(const struct atom_structure atoms_bulk[MAX_ATOMS], size_t  nr_atoms_bulk, const struct latt_par_structure lattice_par, double penetration_depth, double h, double k, double l, double *struc_f_real_ptr, double *struc_f_im_ptr)
{
    double struc_f_real = 0.0;   			/* Real part of the structure factor */
    double struc_f_im = 0.0;   				/* Complex part of the structure factor */
    double abs_exp;             /* Exponential containing absorption effect, e.g. exp(-lattice_par.c/penetration_depth) */
    double divisor;              /* */
    double dividend_real;
    double dividend_im;
    double tpil;                /* 2.0 * PI * l */
    double eef_real;
    double eef_im;
    double epsilon=1.0E-10;   /* Set bulk F to 0 if smaller than this value */

    /* Calculation of CTR at integer positions can result in 0/0 division.
       Simply change l little bit to avoid this. Introduced Version 3.2. dl changed from 1E-8 to 1E-6 Version 3.32*/
    if (fmod(l, 1.0) < 1E-6) {
      l += 1E-6;
      };
    
    calc_structure_factor_domain(&(atoms_bulk[0]), lattice_par, h, k, l, nr_atoms_bulk, &struc_f_real, &struc_f_im);

/*     if (struc_f_real > epsilon ||struc_f_real < -1.0*epsilon || struc_f_im > epsilon || struc_f_im < -1.0*epsilon) {   /\* Only calculate following when bulk is not 0, to avoid nans. *\/  */

    if (penetration_depth >= 0) {
      abs_exp = exp(-1.0 * lattice_par.c / penetration_depth); 
    }
    else {
      abs_exp = 1.0; 
    };
    
    tpil = 2.0*PI*l;
    /*    divisor = 1.0 - 2.0 * cos(tpil)*abs_exp + cos(tpil)*cos(tpil)*abs_exp*abs_exp + sin(tpil)*sin(tpil)*abs_exp*abs_exp; */
    /*    divisor = (1.0 - cos(tpil)*abs_exp) * (1.0 - cos(tpil)*abs_exp) + sin(tpil) *sin(tpil)*abs_exp*abs_exp; */
    divisor = 1.0 - 2.0*cos(tpil)*abs_exp + abs_exp*abs_exp;
    
    /*    if (divisor < 1E-10) divisor = 1E-10;  */
    dividend_real = (1.0 - cos(tpil) *abs_exp); 
    eef_real =  dividend_real / divisor;
    dividend_im = -1.0 * sin(tpil) * abs_exp;
    eef_im = dividend_im / divisor;  
    /*printf(" h=%f k=%f l=%f F=%f %f EEF=%.20f %.20f DIV=%.20f \n", h,k,l,struc_f_real, struc_f_im, dividend_real,dividend_im,divisor ); */

    /* Structure factor to return: bulk structure factor * enhancement factor */
    *struc_f_real_ptr = struc_f_real * eef_real - struc_f_im * eef_im;
    *struc_f_im_ptr = struc_f_im * eef_real + struc_f_real * eef_im;

    
/*     else {     /\* Bulk structure factor set to 0 if struc_f_real and struc_f_im are 0  *\/ */
/*     *struc_f_real_ptr = 0.0; */
/*     *struc_f_im_ptr = 0.0; */
/*     }; */

    return;
}




/********************************************************************************************************
 * double calc_structure_factor(const struct atom_structure  atoms[MAX_ATOMS], const struct atom_structure atoms_bulk[MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth
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
 * h, k, l: Indices of the reflection                                                            	*
 * nr_atoms: Number of atoms in the structure file	                                                *
 * domain_occ: Occupation of domain                                                                     *
 * nr_domains: number of domains                                                                        *
 * nr_atoms_bulk: Number of atoms in the bulk structure                                                 *
 ********************************************************************************************************/

double calc_structure_factor(struct atom_structure atoms[][MAX_ATOMS], const struct atom_structure atoms_bulk[MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth, const double h, const double k, const double l, const size_t nr_atoms[], const double domain_occ[], const size_t nr_domains, size_t nr_atoms_bulk)
{ 
    unsigned int domain_counter;		/* Counter for loop over domains */
    double struc_f_square = 0.0;		/* For summing different domains */ 
    double struc_f_real = 0.0;   			/* Real part of the structure factor */
    double struc_f_im = 0.0;   				/* Complex part of the structure factor */
    double struc_f_bulk_real = 0.0;   			/* Real part of the bulk structure factor */
    double struc_f_bulk_im = 0.0;   				/* Complex part of the bulk structure factor */
    
    /* Calculate bulk CRT */
    calc_structure_factor_bulk(atoms_bulk, nr_atoms_bulk, lattice_par, penetration_depth, h, k, l, &struc_f_bulk_real, &struc_f_bulk_im);

    /* Loop over domains and add the structure factor incoherently (i.e. absolute values)*/
    for (domain_counter = 0; domain_counter < nr_domains; ++domain_counter) {
        struc_f_real = 0.0;
	struc_f_im = 0.0;
	calc_structure_factor_domain(&(atoms[domain_counter][0]), lattice_par, h, k, l, nr_atoms[domain_counter], &struc_f_real, &struc_f_im);
	struc_f_real = struc_f_real + struc_f_bulk_real;
	struc_f_im = struc_f_im + struc_f_bulk_im;
	struc_f_square += domain_occ[domain_counter] * (struc_f_real*struc_f_real + struc_f_im*struc_f_im); 	/* Sum over domains */
    };
    return sqrt(struc_f_square);
}

/****************************************************************************************************************
 * void get_pos_dw_occ(const gsl_vector * x, void *params, double *r1_ptr, double *r2_ptr, double *r3_ptr,	*
 *  double *dw_ptr, double *occ_ptr, size_t atom_counter, size_t domain_counter)			*
 *														*
 * For fit - calculate the current positions, DW factor and occupation from starting positions etc, parameters  *
 * and starting parameters											*
 *
 * x: Vector with displacement parameters used in fit.								*
 * params: Structure with all information 									*
 * r1_ptr, r2_ptr, r3_ptr, dw_ptr, aniso_dw_ptr, occ_ptr: Pointers to to variables to change					*
 * atom_counter: Number of current atom										*
 * domain_counter: Number of current domain									*
 ****************************************************************************************************************/

void get_pos_dw_occ(const gsl_vector * x, void *params, double *r1_ptr, double *r2_ptr, double *r3_ptr, double *dw_ptr, double *aniso_dw_ptr, double *occ_ptr, size_t atom_counter, size_t domain_counter)
{   
    struct atom_structure (*atom_ptr)[MAX_ATOMS];	/* Pointer to array with atomic data */
    struct displacements (*displ_ptr)[MAX_ATOMS];	/* Pointer to array with atomic displacements vectors */
    double (*start_par_ptr); 				/* Pointer to starting parameters */
    int (*fixed_par_ptr);				/* Pointer to array with flag for fixed parameters */
    int i;

    atom_ptr = ((struct fit_all *) params)->atoms_ptr;
    displ_ptr = ((struct fit_all *) params)->displ_ptr;
    fixed_par_ptr = ((struct fit_all *) params)->fixed_par_ptr;
    start_par_ptr = ((struct fit_all *) params)->start_par_ptr;

#define PAR ((struct fit_all *) params)
#define CURR_ATOM  ((*(atom_ptr + domain_counter))[atom_counter]) 
#define CURR_DISPL  ((*(displ_ptr + domain_counter))[atom_counter]) 
#define CURR_X1 gsl_vector_get(x, *((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr1)))
#define CURR_X2 gsl_vector_get(x, *((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr2)))
#define CURR_X3 gsl_vector_get(x, *((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr3)))
#define CURR_X_DW gsl_vector_get(x, *((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr_dw)))
#define CURR_X_ANISO_DW0 gsl_vector_get(x, *((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr_aniso_dw[0])))
#define CURR_X_ANISO_DW1 gsl_vector_get(x, *((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr_aniso_dw[1])))
#define CURR_X_ANISO_DW2 gsl_vector_get(x, *((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr_aniso_dw[2])))
#define CURR_X_ANISO_DW3 gsl_vector_get(x, *((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr_aniso_dw[3])))
#define CURR_X_ANISO_DW4 gsl_vector_get(x, *((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr_aniso_dw[4])))
#define CURR_X_ANISO_DW5 gsl_vector_get(x, *((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr_aniso_dw[5])))
#define CURR_X_OCC gsl_vector_get(x, *((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr_occ)))

    /* Calculate current position of atom */

    *r1_ptr = CURR_ATOM.x; 	/* Starting positions */
    *r2_ptr = CURR_ATOM.y;
    *r3_ptr = CURR_ATOM.z;
    *dw_ptr = CURR_ATOM.dw_par;
    for (i=0; i<6; i++) {
      *(aniso_dw_ptr+i) = CURR_ATOM.aniso_dw[i];
    };
    *occ_ptr = CURR_ATOM.occ;
	
    if (CURR_DISPL.nr1 > 0) {		/* First displacement, if defined */
	if (*(fixed_par_ptr + CURR_DISPL.nr1) == FALSE) {	/* Parameter not fixed */
	    *r1_ptr +=  CURR_X1 * CURR_DISPL.dx1; 
	    *r2_ptr +=  CURR_X1 * CURR_DISPL.dy1; 
	    *r3_ptr +=  CURR_X1 * CURR_DISPL.dz1; 
	};
	if (*(fixed_par_ptr + CURR_DISPL.nr1) == TRUE) {	/* Parameter fixed */
	    *r1_ptr += *(start_par_ptr + CURR_DISPL.nr1) *  CURR_DISPL.dx1;
	    *r2_ptr+= *(start_par_ptr + CURR_DISPL.nr1) *  CURR_DISPL.dy1;
	    *r3_ptr += *(start_par_ptr + CURR_DISPL.nr1) *  CURR_DISPL.dz1;
	};
    };

    if (CURR_DISPL.nr2 > 0) {		/* Second displacement, if defined */
	if (*(fixed_par_ptr + CURR_DISPL.nr2) == FALSE) {	/* Parameter not fixed */
	    *r1_ptr +=  CURR_X2 * CURR_DISPL.dx2; 
	    *r2_ptr +=  CURR_X2 * CURR_DISPL.dy2; 
	    *r3_ptr +=  CURR_X2 * CURR_DISPL.dz2; 
	};
	if (*(fixed_par_ptr + CURR_DISPL.nr2) == TRUE) {	/* Parameter fixed */
	    *r1_ptr += *(start_par_ptr + CURR_DISPL.nr2) *  CURR_DISPL.dx2;
	    *r2_ptr += *(start_par_ptr + CURR_DISPL.nr2) *  CURR_DISPL.dy2;
	    *r3_ptr += *(start_par_ptr + CURR_DISPL.nr2) *  CURR_DISPL.dz2;
	};
    };
	

    if (CURR_DISPL.nr3 > 0) {		/* Third displacement, if defined */
	if (*(fixed_par_ptr + CURR_DISPL.nr3) == FALSE) {	/* Parameter not fixed */
	    *r1_ptr +=  CURR_X3 * CURR_DISPL.dx3; 
	    *r2_ptr+=  CURR_X3 * CURR_DISPL.dy3; 
	    *r3_ptr +=  CURR_X3 * CURR_DISPL.dz3; 
	};
	if (*(fixed_par_ptr + CURR_DISPL.nr3) == TRUE) {	/* Parameter fixed */
	    *r1_ptr += *(start_par_ptr + CURR_DISPL.nr3) *  CURR_DISPL.dx3;
	    *r2_ptr += *(start_par_ptr + CURR_DISPL.nr3) *  CURR_DISPL.dy3;
	    *r3_ptr += *(start_par_ptr + CURR_DISPL.nr3) *  CURR_DISPL.dz3;
	};
    };
	
    if (CURR_DISPL.nr_dw > 0) {		/* Debye-Waller parameter, if defined */
	if (*(fixed_par_ptr + CURR_DISPL.nr_dw) == FALSE) {	/* Parameter not fixed */
	    *dw_ptr += CURR_X_DW * CURR_DISPL.dw_scale;
	};
	if (*(fixed_par_ptr + CURR_DISPL.nr_dw) == TRUE) {	/* Parameter fixed */
	    *dw_ptr += *(start_par_ptr + CURR_DISPL.nr_dw) *  CURR_DISPL.dw_scale;
	};
    };

    /* Parameters for anisotropic Debye-Waller factor */
    if (CURR_DISPL.nr_aniso_dw[0] > 0) {
      if (*(fixed_par_ptr + CURR_DISPL.nr_aniso_dw[0]) == FALSE) {  /* Parameter not fixed */
	for (i=0; i<6; i++) {           
	  *(aniso_dw_ptr+i) += CURR_X_ANISO_DW0 * CURR_DISPL.aniso_dw_scale[0][i];
	};
      };
      if (*(fixed_par_ptr + CURR_DISPL.nr_aniso_dw[0]) == TRUE) {	/* Parameter fixed */
	for (i=0; i<6; i++) {           
	  *(aniso_dw_ptr+i) += *(start_par_ptr + CURR_DISPL.nr_aniso_dw[0]) *  CURR_DISPL.aniso_dw_scale[0][i];
	};
      };
    };
    if (CURR_DISPL.nr_aniso_dw[1] > 0) {
      if (*(fixed_par_ptr + CURR_DISPL.nr_aniso_dw[1]) == FALSE) {  /* Parameter not fixed */
	for (i=0; i<6; i++) {           
	  *(aniso_dw_ptr+i) += CURR_X_ANISO_DW1 * CURR_DISPL.aniso_dw_scale[1][i];
	};
      };
      if (*(fixed_par_ptr + CURR_DISPL.nr_aniso_dw[1]) == TRUE) {	/* Parameter fixed */
	for (i=0; i<6; i++) {           
	  *(aniso_dw_ptr+i) += *(start_par_ptr + CURR_DISPL.nr_aniso_dw[1]) *  CURR_DISPL.aniso_dw_scale[1][i];
	};
      };
    };
    if (CURR_DISPL.nr_aniso_dw[2] > 0) {
      if (*(fixed_par_ptr + CURR_DISPL.nr_aniso_dw[2]) == FALSE) {  /* Parameter not fixed */
	for (i=0; i<6; i++) {           
	  *(aniso_dw_ptr+i) += CURR_X_ANISO_DW2 * CURR_DISPL.aniso_dw_scale[2][i];
	};
      };
      if (*(fixed_par_ptr + CURR_DISPL.nr_aniso_dw[2]) == TRUE) {	/* Parameter fixed */
	for (i=0; i<6; i++) {           
	  *(aniso_dw_ptr+i) += *(start_par_ptr + CURR_DISPL.nr_aniso_dw[2]) *  CURR_DISPL.aniso_dw_scale[2][i];
	};
      };
    };
    if (CURR_DISPL.nr_aniso_dw[3] > 0) {
      if (*(fixed_par_ptr + CURR_DISPL.nr_aniso_dw[3]) == FALSE) {  /* Parameter not fixed */
	for (i=0; i<6; i++) {           
	  *(aniso_dw_ptr+i) += CURR_X_ANISO_DW3 * CURR_DISPL.aniso_dw_scale[3][i];
	};
      };
      if (*(fixed_par_ptr + CURR_DISPL.nr_aniso_dw[3]) == TRUE) {	/* Parameter fixed */
	for (i=0; i<6; i++) {           
	  *(aniso_dw_ptr+i) += *(start_par_ptr + CURR_DISPL.nr_aniso_dw[3]) *  CURR_DISPL.aniso_dw_scale[3][i];
	};
      };
    };
    if (CURR_DISPL.nr_aniso_dw[4] > 0 ){
      if (*(fixed_par_ptr + CURR_DISPL.nr_aniso_dw[4]) == FALSE) {  /* Parameter not fixed */
	for (i=0; i<6; i++) {           
	  *(aniso_dw_ptr+i) += CURR_X_ANISO_DW4 * CURR_DISPL.aniso_dw_scale[4][i];
	};
      };
      if (*(fixed_par_ptr + CURR_DISPL.nr_aniso_dw[4]) == TRUE) {	/* Parameter fixed */
	for (i=0; i<6; i++) {           
	  *(aniso_dw_ptr+i) += *(start_par_ptr + CURR_DISPL.nr_aniso_dw[4]) *  CURR_DISPL.aniso_dw_scale[4][i];
	};
      };
    };
    if (CURR_DISPL.nr_aniso_dw[5] > 0 ){
      if (*(fixed_par_ptr + CURR_DISPL.nr_aniso_dw[5]) == FALSE) {  /* Parameter not fixed */
	for (i=0; i<6; i++) {           
	  *(aniso_dw_ptr+i) += CURR_X_ANISO_DW5 * CURR_DISPL.aniso_dw_scale[5][i];
	};
      };
      if (*(fixed_par_ptr + CURR_DISPL.nr_aniso_dw[5]) == TRUE) {	/* Parameter fixed */
	for (i=0; i<6; i++) {           
	  *(aniso_dw_ptr+i) += *(start_par_ptr + CURR_DISPL.nr_aniso_dw[5]) *  CURR_DISPL.aniso_dw_scale[5][i];
	};
      };
    };
	
    if (CURR_DISPL.nr_occ > 0) {		/* Occupancy, if defined */
      if (*(fixed_par_ptr + CURR_DISPL.nr_occ) == FALSE) {	/* Parameter not fixed */
	*occ_ptr += CURR_X_OCC * 1.0;
      };
      if (*(fixed_par_ptr + CURR_DISPL.nr_occ) == TRUE) {	/* Parameter fixed */
	*occ_ptr += *(start_par_ptr + CURR_DISPL.nr_occ) * 1.0;
      };
    };
    return;
}


/****************************************************************************************************************
 * void calc_f_real_im(const gsl_vector * x, void *params, double *f_real, double *f_im)			*
 * 	                                                                                                        *
 * USED IN FIT                                                         						*
 * Calculate the real and imaginary part of the structure factor (h k l) (returned in *struc_f_real_ptr and 	*
 * struc_f_im_ptr. This function is used by the fit functions. The structure used is the one given in params,	*
 * modified by the displacements given in x and the displacement vectors in params. 				*
 * 														*
 * x: Vector with displacement parameters used in fit.								*
 * params: Structure with following members:									*
 * static struct fit_all {											*
 *	 struct atom_structure (*atoms_ptr)[MAX_ATOMS];	 Pointer to atom data 					*
 *	 struct displacements (*displ_ptr)[MAX_ATOMS];	 Pointer to displacement vectors 			*
 *	 struct structure_factor (*struc_f_ptr);	 Pointer to structure factor data (not used here)	* 
 *	 struct latt_par_structure lattice_par;		 Structure containing lattice parameters (not used here)*
 *	 double (*domain_occ_ptr);			 Pointer to Occupancies of domains (not used here)	*
 *	 size_t (*displtox);			 Pointer to array for translating from x to displacement numbers 
 *	 gsl_vector (*fit_par_ptr);			 Pointer to parameters of the fit (only for simulated annealing, not used in least-squares fit (not used here)
 *	 size_t nr_atoms[MAX_DOMAINS];												*
 *	 size_t nr_domains;	(not used here)									*
 *       size_t nr_struc_f;	(not used here)									*
 *	 size_t nr_parameters; 	(not used here)									*
 * } ;														*
 * domain_counter: Number of current domain 									*
 *  h, k, l: Indices of the reflection 										*
 * scattering_length_sq: Square of the scattering length sin(theta)/ lambda					*
 * domain_counter: Number of current domain 									*
 * *struc_f_real_ptr, *struc_f_im_ptr: pointer to variables where real and imaginary part of structure factor	*
 * 				       should be stored.							*
 ****************************************************************************************************************/

static void calc_real_im_struc_f(const gsl_vector * x, void *params, double h, double k, double l, double scattering_length_sq, size_t domain_counter,  double *struc_f_real_ptr, double *struc_f_im_ptr)
{
    size_t atom_counter; 				/* Atom number used in loop */
    double dw_factor;   			        /* debye-waller factor: exp(-1*dw_par*Q^2) */
    double aniso_dw_factor;   		        /* anisotropic debye-waller factor: */
    double atom_scattering_f = 0.0;	    	        /* Atomic scattering factor */
    double f_real = 0.0;   			/* Real part of the structure factor */
    double f_im = 0.0;   			/* Complex part of the structure factor */
    double r1;					/* r1 coordinate of current atom */
    double r2;					/* r2 coordinate of current atom */
    double r3;					/* r3 coordinate of current atom */
    double dw;					/* Debye-Waller parameter of current atom */
    double aniso_dw[6];                         /* Components of anisotropic Debye-Waller factor (beta11,beta22,beta33,beta12,beta13,beta23) */
    double occ;					/* Occupancy of current atom */

    struct atom_structure (*atom_ptr)[MAX_ATOMS];	/* Pointer to array with atomic data */
    struct displacements (*displ_ptr)[MAX_ATOMS];	/* Pointer to array with atomic displacements vectors */
    double (*domain_occ_ptr);				/* Pointer to array with domain occupancies */
    double (*start_par_ptr); 				/* Pointer to starting parameters */
    int (*fixed_par_ptr);				/* Pointer to array with flag for fixed parameters */

    atom_ptr = ((struct fit_all *) params)->atoms_ptr;
    displ_ptr = ((struct fit_all *) params)->displ_ptr;
    domain_occ_ptr = ((struct fit_all *) params)->domain_occ_ptr; 
    fixed_par_ptr = ((struct fit_all *) params)->fixed_par_ptr;
    start_par_ptr = ((struct fit_all *) params)->start_par_ptr;



#define CURR_ATOM  ((*(atom_ptr + domain_counter))[atom_counter]) 


/* Loop over all atoms in unit cell and sum the structure factor */
    for (atom_counter = 0; atom_counter < *((PAR->nr_atoms_ptr)+domain_counter); ++atom_counter) {
	
	/* Calculate current position of atom */
	get_pos_dw_occ (x, params, &r1, &r2, &r3, &dw, &aniso_dw[0], &occ, atom_counter, domain_counter);

	/* Add the contribution of atom j to the structure factor. */
	atom_scattering_f = atom_scatt_f(CURR_ATOM.el_nr, scattering_length_sq); /* Atomic scattering factor */
	dw_factor = exp(-1.0 * dw * scattering_length_sq);	/* Debye-Waller factor  */
	aniso_dw_factor = exp(-1.0 * (aniso_dw[0]*h*h + aniso_dw[1]*k*k + aniso_dw[2]*l*l + 2.0*aniso_dw[3]*h*k + 2.0*aniso_dw[4]*h*l + 2.0*aniso_dw[5]*k*l)); /* Anisotropic Debye-Waller factor */
	f_real += (occ * atom_scattering_f * dw_factor * aniso_dw_factor * cos(2.0*PI*(h*r1+k*r2+l*r3))) ;
	f_im += (occ * atom_scattering_f * dw_factor * aniso_dw_factor * sin(2.0*PI*(h*r1+k*r2+l*r3)));
    };	/* End atom loop */

    /* Assign structure factor to return */
    *struc_f_real_ptr = f_real;
    *struc_f_im_ptr = f_im;
    
    return;
}

/********************************************************************************************************
 * double calc_patterson(struct structure_factor struct_f[], double r1, double r2, size_t nr_f)	*
 * 													*
 * Calculates the Patterson map for a set of (h,k) reflections for position r1,r2.	 		*
 * (No scale factor for unit cell size)									*
 * Returns the Patterson value at position r1,r2							*
 *													*
 * struct_f: Matrix of structures containing the structure factors for h, k, l				*
 * 	.h												*
 *	.k												*
 *	.l												*
 *	.value												*
 *	.sigma												*
 * r1, r2: r1 and r2 position of the Patterson map (reduced coordinates)				*
 * nr_f: Number of structure factors									*
 ********************************************************************************************************/

double calc_patterson(double r1, double r2, struct structure_factor struc_f[], size_t nr_f)
{
    double patterson_value = 0.0;		/* Value of the Patterson map at r1,r2 */
    size_t f_counter;			/* Counter used to sum all structure factors */
    size_t matrix_counter;        	/* Counter used for the symmetry matrices */
    double h;                             /* Current h */
    double k;                             /* Current k */
    double value = 0.0;                     /* Current structure factor value */



    /**** Calculate the value of the Patterson map at r1, r2 *****/

    /* Loop over all reflections */
    for (f_counter = 0; f_counter < nr_f; ++f_counter) {
 
	/* generate symmetrical equivalent reflections */
	for (matrix_counter = 0; matrix_counter < nr_symmetry_matrices; ++matrix_counter) {
	    h = struc_f[f_counter].h * symmetry_matrix[matrix_counter][0][0] + struc_f[f_counter].k * symmetry_matrix[matrix_counter][0][1];  
	    k = struc_f[f_counter].h * symmetry_matrix[matrix_counter][1][0] + struc_f[f_counter].k * symmetry_matrix[matrix_counter][1][1];
	    value = struc_f[f_counter].value;
      
	    /* Sum all reflections for value of patterson map */ 
	    patterson_value = patterson_value + value*value * cos(2.0*PI * (h * r1 + k * r2));

#ifdef DEBUG1
	    printf("r1: %lf. r2: %lf h: %lf, k: %lf, patterson value: %lf, struct f: %lf, cos: %lf\n", r1, r2, h, k,patterson_value, value, cos(2.0*PI * (h * r1 + k * r2)));
	    fflush;
#endif /*DEBUG1 */

	};

    };
    return patterson_value;

}



/********************************************************************************************************
 * void patterson_map(char filename[], unsigned int step_r1, unsigned int step_r2, size_t nr_struc_f, char comment[])
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

void patterson_map(char filename[], unsigned long int step_r1, unsigned long int step_r2, size_t nr_struc_f, char comment[])
{
    unsigned int patt_r1_counter;	        /* Counts steps of the Patterson map in the r1 direction */
    unsigned int patt_r2_counter;		/* Counts steps of the Patterson map in the r2 direction */
  
    double patt_r1;			/* Current r1 value */
    double patt_r2;			/* Current r2 value */
    double patt_value;                  /* Current value of the Patterson function at r1,r2 */
    int number;                     	/* Number of charcters written */

    FILE *patterson_file;		/* Output file for the Patterson function */    
  
  
    /* Open file for Patterson map */
    patterson_file = fopen(filename,"w");
    if (patterson_file == NULL) { 			/* Check for error */
	fprintf(stderr, "Error: Could not open file %s for output of Patterson map!\n", filename);
	fprintf(stdout, "Error: Could not open file %s for output of Patterson map!\n", filename);
	exit(8);
    };
  
    number = fprintf(patterson_file, comment);					/* Write comment in first line */

    /* Loops for calculating the Patterson map and writing it to output file */
    for (patt_r1_counter=0; patt_r1_counter < step_r1; ++patt_r1_counter) {   /* Loop over the r1 direction */
	patt_r1 = (double)patt_r1_counter / (double)step_r1;				  /* Current r1 value */
	for (patt_r2_counter=0; patt_r2_counter < step_r2; ++patt_r2_counter) {      /* Loop over the r2 direction */
	    patt_r2 =  (double)patt_r2_counter / (double)step_r2;				  	/* Current r2 value */
      
	    /* Calculate value of the patterson map at (patt_r1. patt_r2) */
	    patt_value = calc_patterson(patt_r1, patt_r2, struc_f_exp, nr_struc_f);
      
	    /* Write current value of patterson map to file */
	    number = fprintf(patterson_file, "%f %f %f\n", patt_r1 , patt_r2, patt_value);
	    if (number ==-1) {							/* Check for error */
		fprintf (stderr, "Error: Problem writing to file for Patterson map: %s !", filename);
		fprintf (stdout, "Error: Problem writing to file for Patterson map: %s !", filename);
		exit(8);
	    };
	};
    };
    /* Close file for Patterson map */
    fclose(patterson_file);
    return;
}


/********************************************************************************************************
 * calc_chi_square(struc structure_factor struct_f_exp, struc structure_factor struc_f_calc, _          *
 * _ size_t nr_struc_f, size_t nr_par)                                                    	*
 *												        *
 * Calculates chi^2 to check the agreement of experimental and calculated structure factors.   		*
 * A scale factor is varied for mimimum chi^2.                                                          *
 * Assumes that h,k,l values for the same position in the matrix are same.		                *
 *			
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!  DOES NOT WORK!    !!!!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!
 * 
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

/* double calc_chi_square(size_t nr_struc_f, size_t nr_par, double initial_scale)
{
    double chi_square=0.0;	*/		/* Chi^2    */
/*    size_t struc_f_counter; */			/* Counter for the summation over all structure factors */
    /*  double scale;				

    scale = initial_scale; 
*/
      /* Calculate Chi^2 */
/*   for (struc_f_counter = 0; struc_f_counter < nr_struc_f; ++struc_f_counter) {  
	chi_square = chi_square + ((struc_f_exp[struc_f_counter].value*struc_f_exp[struc_f_counter].value -  struc_f_calc[struc_f_counter].value*struc_f_calc[struc_f_counter].value)   * (struc_f_exp[struc_f_counter].value*struc_f_exp[struc_f_counter].value-struc_f_calc[struc_f_counter].value*struc_f_calc[struc_f_counter].value))/struc_f_exp[struc_f_counter].sigma  ;
    };
  
    chi_square = chi_square / (nr_struc_f - nr_par);
    return chi_square;
}  */

/********************************************************************************************************
 * double calc_chi_square(const struct structure_factor exp_struc_f[], const size_t nr_f,		*
 *  const struct atom_structure atoms[][MAX_ATOMS],  const size_t nr_atoms[],			*
 *  const  double domain_occ[], unsigned  int nr_domains)				*
 *													*
 * Calculates and returns chi, with chi^2  = sum ((f_exp - |f_calc|)^2 / sigma_exp^2)*
 *													*
 * !!! This is chi for structure factors, not  intensities !!!!						*
 * 													*
 *  struct_f: Matrix of structures containing the structure factors for h, k, l				*
 * 	.h												*
 *	.k												*
 *	.l												*
 *	.value												*
 *	.sigma												*
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
 ********************************************************************************************************/                        
double calc_chi_square(const struct structure_factor exp_struc_f[], const size_t nr_f, struct atom_structure atoms[][MAX_ATOMS], const struct atom_structure atoms_bulk[MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth, const size_t nr_atoms[], const  double domain_occ[], size_t nr_domains, size_t nr_atoms_bulk) 
{
    size_t f_counter;		
    double struc_factor;		/* Absolute value of structure factor */ 
    double chi_square = 0.0;		/* chi_square */
    double scale;			/* Scale between experimental and calculated structure factors */

    /* Find scale */
    scale = find_scale(exp_struc_f, nr_f, atoms, atoms_bulk, lattice_par, penetration_depth, nr_atoms, domain_occ, nr_domains, nr_atoms_bulk);

    for (f_counter = 0; f_counter < nr_f; f_counter++) {	/* Loop over all structure factors */

        struc_factor = calc_structure_factor(atoms, atoms_bulk, lattice_par, penetration_depth, exp_struc_f[f_counter].h, exp_struc_f[f_counter].k, exp_struc_f[f_counter].l, nr_atoms, domain_occ, nr_domains, nr_atoms_bulk);

	/* Sum for getting chi^2 */
	chi_square += (exp_struc_f[f_counter].value - struc_factor * scale) *  (exp_struc_f[f_counter].value - struc_factor * scale) / (exp_struc_f[f_counter].sigma * exp_struc_f[f_counter].sigma);
    };
    

    return chi_square;
}



/********************************************************************************************************
 * double calc_r_factor(const struct structure_factor exp_struc_f[], const size_t nr_f,	        	*
 *  const struct atom_structure atoms[][MAX_ATOMS],  const size_t nr_atoms[],			*
 *  const  double domain_occ[], unsigned  int nr_domains) 						*
 *													*
 * Calculates and returns the value of the R-factor R = (sum |f_exp - |f_calc||) / (sum f_exp)		*
 *													*
 *  struct_f: Matrix of structures containing the structure factors for h, k, l				*
 * 	.h												*
 *	.k												*
 *	.l												*
 *	.value												*
 *	.sigma												*
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
 ********************************************************************************************************/                             										
double calc_r_factor(const struct structure_factor exp_struc_f[], const size_t nr_f, struct atom_structure atoms[][MAX_ATOMS], const struct atom_structure atoms_bulk[MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth, const size_t nr_atoms[], const  double domain_occ[], size_t nr_domains, size_t nr_atoms_bulk) 
{
    size_t f_counter;		
    double sum_exp_f = 0.0;		/* Sum of all experimental structure factors */
    double sum_diff_f = 0.0;		/* Sum of F_exp - F_calc for all structure factors */
    double r_factor;			/* R-factor */
    double scale;			/* Scale between experimental and calculated structure factors */
    double struc_factor;		/* Absolute value of structure factor */ 
 
    /* Find scale */
    scale = find_scale(exp_struc_f, nr_f, atoms, atoms_bulk, lattice_par, penetration_depth, nr_atoms, domain_occ, nr_domains, nr_atoms_bulk);

    for (f_counter = 0; f_counter < nr_f; f_counter++) {	/* Loop over all structure factors */
      
        struc_factor = calc_structure_factor(atoms, atoms_bulk, lattice_par, penetration_depth, exp_struc_f[f_counter].h, exp_struc_f[f_counter].k, exp_struc_f[f_counter].l, nr_atoms, domain_occ, nr_domains, nr_atoms_bulk);

	sum_exp_f += exp_struc_f[f_counter].value;	/* Sum experimental structure factors */ 
	sum_diff_f += fabs(exp_struc_f[f_counter].value - struc_factor * scale); /* Sum F_exp - F_calc */
    };
    
    /* R-factor */
    r_factor = sum_diff_f / sum_exp_f;

    return r_factor;
}
	

/********************************************************************************
 * void set_symmetry_matrix(enum SURFACE_SYMMETRY symmetry_to_set)     		*
 *                                                                     		*
 * Set symmetry matrices                                               		*
 *                                                                     		*
 * This sets the following global variables:					*
 * - the global symmetry matrix symmetry_matrix,                   		*
 * - the number of matrices in the symmetry nr_symmetry_matrices.              	*
 *                                                                     		*
 * symmetry_to_set gives the symmetry the symmetry matrix should be set to.	*
 * It can be at present:				               		*
 * - p6 or									*
 * - p6mm									*
 ********************************************************************************/

void set_symmetry_matrix(enum SURFACE_SYMMETRY symmetry_to_set)
{


    switch (symmetry_to_set) 
	{
	case p1:
		nr_symmetry_matrices = 1;
		double sym_matrixp1[1][2][2] = {{{1.0, 0.0}, /* 1 */
	        				          {0.0, 1.0}}};
		
                {
		unsigned int counterp1;
      
		for(counterp1 = 0; counterp1 < nr_symmetry_matrices; ++counterp1) { 
		    symmetry_matrix[counterp1][0][0] = sym_matrixp1[counterp1][0][0];
		    symmetry_matrix[counterp1][0][1] = sym_matrixp1[counterp1][0][1];
		    symmetry_matrix[counterp1][1][0] = sym_matrixp1[counterp1][1][0];
		    symmetry_matrix[counterp1][1][1] = sym_matrixp1[counterp1][1][1];
		};
	    };
	    break;
	
	case p6mm:
	    nr_symmetry_matrices = 12;
      
	    /* Symmetry matrices for p6mm (assuming hexagonal unit vectors */
	    double sym_matrixp6mm[12][2][2] ={{{1.0, 0.0}, /* 1 */
					       {0.0, 1.0}},
			       
					      {{0.0, 1.0}, /* 2 */
					       {1.0, 0.0}},
			       
					      {{1.0, 1.0}, /* 3 */
					       {-1.0, 0.0}},
			       
					      {{-1.0, 0.0}, /* 4 */
					       {1.0, 1.0}},
			       
					      {{0.0, 1.0},  /* 5 */
					       {-1.0, -1.0}}, 
			       
					      {{-1.0, -1.0}, /* 6 */
					       {0.0, 1.0}},
			       
					      {{-1.0, 0.0}, /* 7 */
					       {0.0, -1.0}},
			       
					      {{0.0, -1.0}, /* 8 */
					       {-1.0, 0.0}},
			       
					      {{-1.0, -1.0}, /* 9 */
					       {1.0, 0.0}},
			       
					      {{1.0, 0.0}, /* 10 */
					       {-1.0, -1.0}},
			       
					      {{0.0, -1.0}, /* 11 */
					       {1.0, 1.0}},
			       
					      {{1.0, 1.0}, /* 12 */
					       {0.0, -1.0}}
	    };

	    {	
		unsigned int counterp6mm;
      
		for(counterp6mm = 0; counterp6mm < nr_symmetry_matrices; ++counterp6mm) { 
		    symmetry_matrix[counterp6mm][0][0] = sym_matrixp6mm[counterp6mm][0][0];
		    symmetry_matrix[counterp6mm][0][1] = sym_matrixp6mm[counterp6mm][0][1];
		    symmetry_matrix[counterp6mm][1][0] = sym_matrixp6mm[counterp6mm][1][0];
		    symmetry_matrix[counterp6mm][1][1] = sym_matrixp6mm[counterp6mm][1][1];
		};
	    };
	    break;

	case p6:
	    nr_symmetry_matrices = 6;

	    /* Symmetry matrices for p6? (assuming hexagonal unit vectors) */
	    double sym_matrixp6[6][2][2] = {{{1.0, 0.0}, /* 1 */
					     {0.0, 1.0}},
			       
					    {{1.0, 1.0}, /* 2 */
					     {-1.0, 0.0}},
			       
					    {{0.0, 1.0},  /* 3 */
					     {-1.0, -1.0}}, 
			       
					    {{-1.0, 0.0}, /* 4 */
					     {0.0, -1.0}},
				    
					    {{-1.0, -1.0}, /* 5 */
					     {1.0, 0.0}},
			       
					    {{0.0, -1.0}, /* 6 */
					     {1.0, 1.0}},
			       
	    };

	    {
		unsigned int counterp6;
      
		for(counterp6 = 0; counterp6 < nr_symmetry_matrices; ++counterp6) { 
		    symmetry_matrix[counterp6][0][0] = sym_matrixp6[counterp6][0][0];
		    symmetry_matrix[counterp6][0][1] = sym_matrixp6[counterp6][0][1];
		    symmetry_matrix[counterp6][1][0] = sym_matrixp6[counterp6][1][0];
		    symmetry_matrix[counterp6][1][1] = sym_matrixp6[counterp6][1][1];
		};
	    };
	    break;
	};	
    return;
}

/****************************************************************************************************************
 * static find_nr_param(struct displacements displ, size_t nr_atoms, size_t * displtox, int *xtpdispl)*
 * 														*
 * Find the number of parameter for least-squares fit.								*
 * displ: Structure with displacement vectors.									*
 * fixed_parameters: Array with the numbers of the fixed parameters						*
 * nr_atoms: Number of atoms in unit cell.									*
 * *displ_x: Pointer to translation matrix from user supplied parameter numbers to fit vector x components.	*
 * *x_displ: Pointer to translation matrix from fit vector x components	to user supplied parameter numbers.	*
 ****************************************************************************************************************/

static size_t find_nr_param(struct displacements displ[][MAX_ATOMS], int fixed_parameters[], size_t nr_atoms[], size_t nr_domains, size_t * displtox, int *xtodispl)
{
    size_t param_nr = 0;		/* Number of parameters, first one is scale factor */
    size_t atom_counter;
    size_t domain_counter;
    int param_array[MAX_FREE_PARAMETERS];      /* Stores serial numbers of parameters */
    size_t param_array_counter;
    int param_first_time = TRUE;			/* Was parameter found for first time? */
	
    memset(&(param_array[0]), 0, MAX_FREE_PARAMETERS);      /* Initialize to 0 */
    
    *(displtox + 0) =  0;		/* Scale factor has number 0 */
    *(xtodispl + 0) = 0;
    
    /* 1st displacement for atom */
    for (domain_counter = 0; domain_counter < nr_domains; ++domain_counter) {
	for (atom_counter = 0; atom_counter < nr_atoms[domain_counter]; atom_counter++) {	/* Was parameter nr already found? */ 

	    param_first_time = TRUE;

	    if (displ[domain_counter][atom_counter].nr1 < 1) {	/* Parameter not defined */
		continue;
	    };
	    
	    if (fixed_parameters[displ[domain_counter][atom_counter].nr1] == TRUE) {     /* Was parameter fixed? */
		continue;
	    };

	    for (param_array_counter = 0; param_array_counter < param_nr; param_array_counter++) {
		if (displ[domain_counter][atom_counter].nr1 == param_array[param_array_counter]) {	/* Was parameter nr already found? */ 
		    param_first_time = FALSE;		
		    break;
		};
	    };
	    
		
	    if (param_first_time == TRUE) {			/* If parameter nr not found before, add to param_array */
		param_array[param_nr] = displ[domain_counter][atom_counter].nr1;
		param_nr += 1;
		*(xtodispl + param_nr) = displ[domain_counter][atom_counter].nr1;		/* Translation from vector x for fit to user supplied paramter numbers */
		*(displtox + displ[domain_counter][atom_counter].nr1) = param_nr;		/* Translation from user supplied paramter numbers to vector x for fit */
	    };
	};	
    };

    /* 2nd displacement for atom */
    for (domain_counter = 0; domain_counter < nr_domains; ++domain_counter) {
	for (atom_counter = 0; atom_counter < nr_atoms[domain_counter]; atom_counter++) {	/* Was parameter nr already found? */ 
	  
	    param_first_time = TRUE;
  
	    if (displ[domain_counter][atom_counter].nr2 < 1) {	/* Parameter not defined */
		continue;
	    };

	    if (fixed_parameters[displ[domain_counter][atom_counter].nr2] == TRUE) {     /* Was parameter fixed? */
		continue;
	    };

	    for (param_array_counter = 0; param_array_counter < param_nr; param_array_counter++) {
		if(displ[domain_counter][atom_counter].nr2 == param_array[param_array_counter]) {	/* Was parameter nr already found? */ 
		    param_first_time = FALSE;
		    break;
		};
	    };
	
	    if (param_first_time == TRUE) {			/* If parameter nr not found before, add to param_array */
		param_array[param_nr] = displ[domain_counter][atom_counter].nr2;
		param_nr += 1;
		*(xtodispl + param_nr) = displ[domain_counter][atom_counter].nr2;		/* Translation from vector x for fit to user supplied paramter numbers */
		*(displtox + displ[domain_counter][atom_counter].nr2) = param_nr;		/* Translation from user supplied paramter numbers to vector x for fit */
	    };
	};	
    };

    /* 3rd displacement for atom */
    for (domain_counter = 0; domain_counter < nr_domains; ++domain_counter) {
	for (atom_counter = 0; atom_counter < nr_atoms[domain_counter]; atom_counter++) {	/* Was parameter nr already found? */ 
	    
	    param_first_time = TRUE;

	    if (displ[domain_counter][atom_counter].nr3 < 1) {	/* Parameter not defined */
		continue;
	    };

	    if (fixed_parameters[displ[domain_counter][atom_counter].nr3] == TRUE) {     /* Was parameter fixed? */
		continue;
	    };

	    for (param_array_counter = 0; param_array_counter < param_nr; param_array_counter++) {
		if(displ[domain_counter][atom_counter].nr3 == param_array[param_array_counter]) {		/* Was parameter nr already found? */ 
		    param_first_time = FALSE;
		    break;
		};	
	    };
		
	    if (param_first_time == TRUE) {			/* If parameter nr not found before, add to param_array */
		param_array[param_nr] = displ[domain_counter][atom_counter].nr3;
		param_nr += 1;
		*(xtodispl + param_nr) = displ[domain_counter][atom_counter].nr3;		/* Translation from vector x for fit to user supplied paramter numbers */
		*(displtox + displ[domain_counter][atom_counter].nr3) = param_nr;		/* Translation from user supplied paramter numbers to vector x for fit */
	    };
	};
    };	

    /* Debye-Waller factor for atom */
    for (domain_counter = 0; domain_counter < nr_domains; ++domain_counter) {
	for (atom_counter = 0; atom_counter < nr_atoms[domain_counter]; atom_counter++) {	/* Was parameter nr already found? */ 
	  
	    param_first_time = TRUE;
  
	    if (displ[domain_counter][atom_counter].nr_dw < 1) {	/* Parameter not defined */
		continue;
	    };

	    if (fixed_parameters[displ[domain_counter][atom_counter].nr_dw] == TRUE) {     /* Was parameter fixed? */
		continue;
	    };

	    for (param_array_counter = 0; param_array_counter < param_nr; param_array_counter++) {
		if(displ[domain_counter][atom_counter].nr_dw == param_array[param_array_counter]) {	/* Was parameter nr already found? */ 
		    param_first_time = FALSE;
		    break;
		};
	    };
	
	    if (param_first_time == TRUE) {			/* If parameter nr not found before, add to param_array */
		param_array[param_nr] = displ[domain_counter][atom_counter].nr_dw;
		param_nr += 1;
		*(xtodispl + param_nr) = displ[domain_counter][atom_counter].nr_dw;		/* Translation from vector x for fit to user supplied paramter numbers */
		*(displtox + displ[domain_counter][atom_counter].nr_dw) = param_nr;		/* Translation from user supplied paramter numbers to vector x for fit */
	    };
	};	
    };

    /* Components of anisotropic Debye-Waller factor for atom */
    for (domain_counter = 0; domain_counter < nr_domains; ++domain_counter) {
	for (atom_counter = 0; atom_counter < nr_atoms[domain_counter]; atom_counter++) {	/* Was parameter nr already found? */ 
  	    int ai;
	    param_first_time = TRUE;
  
	    for (ai=0; ai<6; ai++) {
	      if (displ[domain_counter][atom_counter].nr_aniso_dw[ai] < 1) {	/* Parameter not defined */
		continue;
	      };

	      if (fixed_parameters[displ[domain_counter][atom_counter].nr_aniso_dw[ai]] == TRUE) {     /* Was parameter fixed? */
		continue;
	      };

	      for (param_array_counter = 0; param_array_counter < param_nr; param_array_counter++) {
		if(displ[domain_counter][atom_counter].nr_aniso_dw[ai] == param_array[param_array_counter]) {	/* Was parameter nr already found? */ 
		  param_first_time = FALSE;
		  break;
		};
	      };
	      
	      if (param_first_time == TRUE) {			/* If parameter nr not found before, add to param_array */
		param_array[param_nr] = displ[domain_counter][atom_counter].nr_aniso_dw[ai];
		param_nr += 1;
		*(xtodispl + param_nr) = displ[domain_counter][atom_counter].nr_aniso_dw[ai];		/* Translation from vector x for fit to user supplied paramter numbers */
		*(displtox + displ[domain_counter][atom_counter].nr_aniso_dw[ai]) = param_nr;		/* Translation from user supplied paramter numbers to vector x for fit */
	      };
	    };
	};	
    };

    /* Occupancy for atom */
    for (domain_counter = 0; domain_counter < nr_domains; ++domain_counter) {
	for (atom_counter = 0; atom_counter < nr_atoms[domain_counter]; atom_counter++) {	/* Was parameter nr already found? */ 
	  
	    param_first_time = TRUE;
  
	    if (displ[domain_counter][atom_counter].nr_occ < 1) {	/* Parameter not defined */
		continue;
	    };

	    if (fixed_parameters[displ[domain_counter][atom_counter].nr_occ] == TRUE) {     /* Was parameter fixed? */
		continue;
	    };

	    for (param_array_counter = 0; param_array_counter < param_nr; param_array_counter++) {
		if(displ[domain_counter][atom_counter].nr_occ == param_array[param_array_counter]) {	/* Was parameter nr already found? */ 
		    param_first_time = FALSE;
		    break;
		};
	    };
	
	    if (param_first_time == TRUE) {			/* If parameter nr not found before, add to param_array */
		param_array[param_nr] = displ[domain_counter][atom_counter].nr_occ;
		param_nr += 1;
		*(xtodispl + param_nr) = displ[domain_counter][atom_counter].nr_occ;		/* Translation from vector x for fit to user supplied paramter numbers */
		*(displtox + displ[domain_counter][atom_counter].nr_occ) = param_nr;		/* Translation from user supplied paramter numbers to vector x for fit */
	    };
	};	
    };

    return param_nr;
}

	
/****************************************************************************************************************
 * int ls_fit_f(const gsl_vector * x, void *params, gsl_vector * f)						*
 * 														*
 * !!! ONLY FOR HEXAGONAL WITH ANGLE 60 DEGREE BETWEEN H AND K !!!						*
 * 														*
 *														*
 * Calculates |F_calc^2-F_exp^2|/sigma i.e f (F: structure factors) for use with GSL Levenberg-Marquardt		*
 * least-squares fit functions											*
 *														*
 ****************************************************************************************************************/
int ls_fit_f(const gsl_vector * x, void *params, gsl_vector * f)
{

    double scattering_length_sq;			/* The length of the scattering vector (used for
							   calculating the atomic scattering factor) */
    size_t f_counter;					/* Counter over structure factors */
    size_t domain_counter;					/* Counter over surface domains */
    double struc_f_real = 0.0;   			/* Real part of the structure factor */
    double struc_f_im = 0.0;   				/* Complex part of the structure factor */
    double struc_f_square = 0.0;			/* Square of the absolute structure factor */ 
    double residual = 0.0;				/* (F_calc^2-F_exp^2)/(2*sigma*F) */
    double h;					/* index h of current structure factor */
    double k; 					/* index k of current structure factor */
    double l;					/* index l of current structure factor */
    double (*domain_occ_ptr);				/* Pointer to array with domain occupancies */

    

#define PAR ((struct fit_all *) params)
#define CURR_DOMAIN_OCC *(domain_occ_ptr + domain_counter)
#define CURR_F ((struct structure_factor *)((((struct fit_all *) params)->struc_f_ptr) + f_counter))
#define CURR_BULK_F_REAL *((double *)((((struct fit_all *) params)->bulk_structure_factor_real_ptr) + f_counter))
#define CURR_BULK_F_IM *((double *)((((struct fit_all *) params)->bulk_structure_factor_im_ptr) + f_counter))

#define SCALE gsl_vector_get(x, (size_t) 0)

    domain_occ_ptr = ((struct fit_all *) params)->domain_occ_ptr; 
		
    /* Loop over all structure factors */
    for (f_counter = 0; f_counter < PAR->nr_struc_f; f_counter++) {

	/* Initialize struc_fs */
	struc_f_real = 0.0;
	struc_f_im = 0.0;
	struc_f_square = 0.0;


	/* Get h, k, l */
	h = ((struct structure_factor *)((((struct fit_all *) params)->struc_f_ptr) + f_counter))->h;
	k = ((struct structure_factor *)((((struct fit_all *) params)->struc_f_ptr) + f_counter))->k;
	l = ((struct structure_factor *)((((struct fit_all *) params)->struc_f_ptr) + f_counter))->l;

	/* Square of the scattering length */
/*	scattering_length_sq = (4.0/3.0/(PAR->lattice_par.a * PAR->lattice_par.a) * (h*h + k*k + h*k) + l*l / (PAR->lattice_par.c * PAR->lattice_par.c)) / 4.0; 
 */
	scattering_length_sq = calc_scattering_length_sq(PAR->lattice_par, h, k, l);


	/* Calculate the structure factor by summing the terms. */
	for (domain_counter = 0; domain_counter < PAR->nr_domains; ++domain_counter) {

	    struc_f_real = 0.0;
	    struc_f_im = 0.0;


	    calc_real_im_struc_f(x, params, h, k, l, scattering_length_sq, domain_counter, &struc_f_real, &struc_f_im);

	    /* Calculate the square of the absolute value of the structure factor  */
	    struc_f_real = struc_f_real + CURR_BULK_F_REAL;
	    struc_f_im = struc_f_im + CURR_BULK_F_IM;
	    struc_f_square +=  CURR_DOMAIN_OCC * SCALE * (struc_f_real*struc_f_real + struc_f_im*struc_f_im);
	};  /* End domain loop */

	residual = (struc_f_square - CURR_F->value*CURR_F->value)/(2.0 * CURR_F->sigma * CURR_F->value); 
		
	/* Set value of f_counter */
	gsl_vector_set (f, f_counter, residual);
    };
	
    return GSL_SUCCESS;

}

/****************************************************************************************************************
 * int ls_fit_df(const gsl_vector * x, void *params, gsl_matrix * J)  						*
 *														*
 * Calculates d(|F_calc^2-F_exp^2|) / dx_i i.e J (F: structure factors, x_i: free parameter) for use 			*
 * with GSL Levenberg-Marquardt least-squares fit functions							*
 *														*
 ****************************************************************************************************************/

int ls_fit_df(const gsl_vector * x, void *params,  gsl_matrix * J) 
{
    size_t par_counter;				/* Counter over parameters */
    double scattering_length_sq;		/* The length of the scattering vector (used for
						   calculating the atomic scattering factor) */
    size_t f_counter;				/* Counter over structure factors */
    size_t domain_counter;				/* Counter over surface domains */
    size_t atom_counter;  				/* Atom number used in loop */
    int aniso_counter;
    double dw_factor;   			/* debye-waller factor: exp(-1*dw_par*Q^2) */
    double aniso_dw_factor;                     /* Anisotropic Debye-Waller factor */
    double atom_scattering_f = 0.0;	    	/* Atomic scattering factor */
    double struc_f_real = 0.0;   		/* Real part of the structure factor */
    double struc_f_im = 0.0;   			/* Complex part of the structure factor */
    double struc_f_square = 0.0;			/* squared absolute value of the structure factor */
    double h;					/* index h of current structure factor */
    double k; 					/* index k of current structure factor */
    double l;					/* index l of current structure factor */
    double r1;					/* r1 coordinate of current atom */
    double r2;					/* r2 coordinate of current atom */
    double r3;					/* r3 coordinate of current atom */
    double dw;					/* Debye-Waller parameter of the current atom */
    double aniso_dw[6];                         /* Components of anisotropic Debye-Waller factor (beta11,beta22,beta33,beta12,beta13,beta23) */
    double occ;					/* Occupancy for the current atom */
    double df_term1 = 0.0;			/* First term for df */
    double df_term2 = 0.0;			/* Second term for df */
    double curr_df =0.0;			/* For summation of df */
    double curr_j = 0.0;			/* element of J */
    double curr_j0 = 0.0;			/* element of J for scale factor */

    struct atom_structure (*atom_ptr)[MAX_ATOMS];
    struct displacements (*displ_ptr)[MAX_ATOMS];

    int (*fixed_par_ptr);				/* Pointer to array with flag for fixed parameters */
    double (*domain_occ_ptr);			

    atom_ptr = ((struct fit_all *) params)->atoms_ptr;		/* Pointer to array with atom data */    

    displ_ptr = ((struct fit_all *) params)->displ_ptr;		/* Pointer to array with displacements */
    fixed_par_ptr = ((struct fit_all *) params)->fixed_par_ptr;
    domain_occ_ptr = ((struct fit_all *) params)->domain_occ_ptr; /* Pointer to array with domain occupancies */

#define PAR ((struct fit_all *) params)
#define CURR_ATOM  ((*(atom_ptr + domain_counter))[atom_counter]) 
#define CURR_DISPL  ((*(displ_ptr + domain_counter))[atom_counter]) 

#define CURR_F ((struct structure_factor *)((((struct fit_all *) params)->struc_f_ptr) + f_counter))
#define CURR_DOMAIN_OCC *(domain_occ_ptr + domain_counter)
#define CURR_BULK_F_REAL *((double *)((((struct fit_all *) params)->bulk_structure_factor_real_ptr) + f_counter))
#define CURR_BULK_F_IM *((double *)((((struct fit_all *) params)->bulk_structure_factor_im_ptr) + f_counter))
#define SCALE gsl_vector_get(x, (size_t) 0)		

		
    /* Loop over all structure factors */
    for (f_counter = 0; f_counter < PAR->nr_struc_f; f_counter++) { 		

	/* Initialize struct_fs */
	struc_f_real = 0.0;
	struc_f_im = 0.0;
	struc_f_square = 0.0;
	atom_scattering_f =0.0;

	/* Get h, k, l */
	h = ((struct structure_factor *)((((struct fit_all *) params)->struc_f_ptr) + f_counter))->h;
	k = ((struct structure_factor *)((((struct fit_all *) params)->struc_f_ptr) + f_counter))->k;
	l = ((struct structure_factor *)((((struct fit_all *) params)->struc_f_ptr) + f_counter))->l;

		
/*	fprintf (stdout, "h, k, l in ls_fit_df: %f %f %f\n", h, k ,l);	*/
		
	/* Square of the scattering length */
	/* This is (sin(theta)/lambda)^2 = q^2 / (16*pi^2) */ 
/*	scattering_length_sq = (4.0/3.0/(PAR->lattice_par.a * PAR->lattice_par.a) * (h*h + k*k + h*k) + l*l / (PAR->lattice_par.c * PAR->lattice_par.c)) / 4.0; */

	scattering_length_sq = calc_scattering_length_sq(PAR->lattice_par, h, k, l);


	/* Calculate the structure factor by summing the terms. */ 
	struc_f_square = 0.0;
	for (domain_counter = 0; domain_counter < PAR->nr_domains; ++domain_counter) {	/* loop over all domains */
	    struc_f_real = 0.0;
	    struc_f_im = 0.0;

	    calc_real_im_struc_f(x, params, h, k, l, scattering_length_sq, domain_counter, &struc_f_real, &struc_f_im);

	    struc_f_real = struc_f_real + CURR_BULK_F_REAL;
	    struc_f_im = struc_f_im + CURR_BULK_F_IM ;
	    struc_f_square += CURR_DOMAIN_OCC * struc_f_real*struc_f_real + struc_f_im*struc_f_im;	/* Add incoherently over domains */
	};	/* End domain loop 2*/


    /***** Loop over all parameters ******/
	for (par_counter = 0; par_counter < PAR->nr_parameters; par_counter++) {

	    curr_df = 0.0;
	    
	    /* Scale factor */
	    if (par_counter == 0) { 

		curr_j0 = struc_f_square/(2.0 * CURR_F->sigma * CURR_F->value);
		
		gsl_matrix_set(J, f_counter, (size_t) 0, curr_j0);	/* Set Jacobian matrix element */
		continue; 
	    };
		
	    /***** Calculate dI_hkl / d(x_(par_counter)) *******/
	    /* Find all displacements with number par_counter */
	    
	    for (domain_counter = 0; domain_counter < PAR->nr_domains; ++domain_counter) {	
		struc_f_real = 0.0;
		struc_f_im = 0.0;

/* Calculate the structure factor of the current domain. */ 
		calc_real_im_struc_f(x, params, h, k, l, scattering_length_sq, domain_counter, &struc_f_real, &struc_f_im);


		/* Loop over all atoms to sum the contribution to df from this domain */
		for (atom_counter = 0; atom_counter < *((PAR->nr_atoms_ptr)+domain_counter); atom_counter++) {

		    atom_scattering_f = atom_scatt_f(CURR_ATOM.el_nr, scattering_length_sq); /* Atomic scattering factor */

                /* Calculate current position of atom */
		    get_pos_dw_occ (x, params, &r1, &r2, &r3, &dw, &aniso_dw[0], &occ, atom_counter, domain_counter);

		    dw_factor = exp(-1.0 * dw * scattering_length_sq);	/* Debye-Waller factor  */
		    aniso_dw_factor = exp(-1.0 * (aniso_dw[0]*h*h + aniso_dw[1]*k*k + aniso_dw[2]*l*l + 2.0*aniso_dw[3]*h*k + 2.0*aniso_dw[4]*h*l + 2.0*aniso_dw[5]*k*l));  /* Anisotropic Debye-Waller factor */

		    if (CURR_DISPL.nr1 > 0 && (*(fixed_par_ptr + CURR_DISPL.nr1) == FALSE)) {		/* Is displacement paramtere nr1 defined and not fixed for current atom? */
			/* Is current parameter same as displacement parameter nr1 of current atom ? */
			if (*((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr1)) == par_counter) {
			    df_term1 = -1.0 * (occ * atom_scattering_f * dw_factor * aniso_dw_factor* sin(2.0*PI*(h*r1+k*r2+l*r3))) * (h * CURR_DISPL.dx1 + k * CURR_DISPL.dy1 + l * CURR_DISPL.dz1) * (struc_f_real + CURR_BULK_F_REAL);
			    df_term2 = (occ * atom_scattering_f * dw_factor * aniso_dw_factor* cos(2.0*PI*(h*r1+k*r2+l*r3))) * (h * CURR_DISPL.dx1 + k * CURR_DISPL.dy1 + l * CURR_DISPL.dz1) * (struc_f_im + CURR_BULK_F_IM);
			    curr_df += 2.0 * CURR_DOMAIN_OCC * SCALE * 2.0 * PI * (df_term1 + df_term2);
			};
		    };
				
		    if (CURR_DISPL.nr2 > 0 && (*(fixed_par_ptr + CURR_DISPL.nr2) == FALSE)) {		/* Is displacement parameter nr2 defined and not fixed for current atom? */
			/* Is current parameter same as displacement parameter nr2 of current atom ? */
			if (*((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr2)) == par_counter) {
			    df_term1 = -1.0 * (occ * atom_scattering_f * dw_factor * aniso_dw_factor* sin(2.0*PI*(h*r1+k*r2+l*r3))) * (h * CURR_DISPL.dx2 + k * CURR_DISPL.dy2 + l * CURR_DISPL.dz2) * (struc_f_real + CURR_BULK_F_REAL);
			    df_term2 = (occ * atom_scattering_f * dw_factor * aniso_dw_factor* cos(2.0*PI*(h*r1+k*r2+l*r3))) * (h * CURR_DISPL.dx2 + k * CURR_DISPL.dy2 + l * CURR_DISPL.dz2) * (struc_f_im + CURR_BULK_F_IM);
			    curr_df += 2.0 * CURR_DOMAIN_OCC * SCALE * 2.0 * PI * (df_term1 + df_term2);
			};
		    };

		    if (CURR_DISPL.nr3 > 0 && (*(fixed_par_ptr + CURR_DISPL.nr3) == FALSE)) {		/* Is displacement paramtere nr3 defined and not fixed for current atom? */
			/* Is current parameter same as displacement parameter nr3 of current atom ? */
			if (*((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr3)) == par_counter) {
			    df_term1 = -1.0 * (occ * atom_scattering_f * dw_factor * aniso_dw_factor * sin(2.0*PI*(h*r1+k*r2+l*r3))) * (h * CURR_DISPL.dx3 + k * CURR_DISPL.dy3 + l * CURR_DISPL.dz3) * (struc_f_real + CURR_BULK_F_REAL);
			    df_term2 =  (occ * atom_scattering_f * dw_factor * aniso_dw_factor * cos(2.0*PI*(h*r1+k*r2+l*r3))) * (h * CURR_DISPL.dx3 + k * CURR_DISPL.dy3 + l * CURR_DISPL.dz3) * (struc_f_im+ CURR_BULK_F_IM);
			    curr_df += 2.0 * CURR_DOMAIN_OCC * SCALE * 2.0 * PI * (df_term1 + df_term2);
			};
		    };

		    if (CURR_DISPL.nr_dw > 0 && (*(fixed_par_ptr + CURR_DISPL.nr_dw) == FALSE)) {		/* Is Debye-Waller fit parameter defined and not fixed for current atom? */
			/* Is current parameter same as Debye-Waller fit parameter of current atom ? */	
			if (*((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr_dw)) == par_counter) {
			    df_term1 = occ * atom_scattering_f * dw_factor * aniso_dw_factor *  cos(2.0*PI*(h*r1+k*r2+l*r3)) * (struc_f_real + CURR_BULK_F_REAL);
			    df_term2 = occ * atom_scattering_f * dw_factor * aniso_dw_factor * sin(2.0*PI*(h*r1+k*r2+l*r3)) * (struc_f_im + CURR_BULK_F_IM);

/* !!!!!! Commented line gives correct errors for Debye-Weller factor ! But Debye-Waller factor is included too much in fit , so fit does not converge well with regard to positions! 
With used line, the component of the gradient for Debye-Waller factors is a factor of dw_scale too small. This seems to improve convergence. 
DOES IT INTRODUCE SIGNIFICANT ERRORS ???????
If scale factor for Debye-Waller factor is 1, result should be same. Best to use this once good fit has been found. 
!!!!!!!!!!!!!!!!!! */
/*			    curr_df +=  -2.0 * CURR_DISPL.dw_scale * CURR_DOMAIN_OCC * SCALE * scattering_length_sq *  (df_term1 + df_term2); */
			    curr_df +=  -2.0 * CURR_DOMAIN_OCC * SCALE * scattering_length_sq *  (df_term1 + df_term2);
			};
		    };

		    for (aniso_counter = 0;aniso_counter <6; aniso_counter++) {
		      if (CURR_DISPL.nr_aniso_dw[aniso_counter] > 0 && (*(fixed_par_ptr + CURR_DISPL.nr_aniso_dw[aniso_counter]) == FALSE)) {		/* Is anisotropic Debye-Waller fit parameter defined and not fixed for current atom? */
			/* Is current parameter same as anisotropic Debye-Waller fit parameter of current atom ? */	
			if (*((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr_aniso_dw[aniso_counter])) == par_counter) {
			  double aniso_sum;

			  df_term1 =  occ * atom_scattering_f * dw_factor * aniso_dw_factor * cos(2.0*PI*(h*r1+k*r2+l*r3)) * (struc_f_real + CURR_BULK_F_REAL);
			  df_term2 = occ * atom_scattering_f * dw_factor * aniso_dw_factor *  sin(2.0*PI*(h*r1+k*r2+l*r3)) * (struc_f_im + CURR_BULK_F_IM);
			  aniso_sum =  CURR_DISPL.aniso_dw_scale[aniso_counter][0]*h*h + CURR_DISPL.aniso_dw_scale[aniso_counter][1]*k*k + CURR_DISPL.aniso_dw_scale[aniso_counter][2]*l*l + 2.0*CURR_DISPL.aniso_dw_scale[aniso_counter][3]*h*k  + 2.0*CURR_DISPL.aniso_dw_scale[aniso_counter][4]*h*l + 2.0*CURR_DISPL.aniso_dw_scale[aniso_counter][5]*k*l;
			  curr_df +=  -2.0 * CURR_DOMAIN_OCC * SCALE * aniso_sum * scattering_length_sq *  (df_term1 + df_term2);
			};
		      };
		    };


		    if (CURR_DISPL.nr_occ > 0 && (*(fixed_par_ptr + CURR_DISPL.nr_occ) == FALSE)) {	/* Is occupancy fit parameter defined and not fixed for current atom? */
			/* Is current parameter same as occupancy fit parameter of current atom ? */	
			if (*((size_t * )((((struct fit_all *) params)->displtox) + CURR_DISPL.nr_occ)) == par_counter) {
			    df_term1 = atom_scattering_f * dw_factor * aniso_dw_factor *  cos(2.0*PI*(h*r1+k*r2+l*r3)) * (struc_f_real + CURR_BULK_F_REAL);
			    df_term2 = atom_scattering_f * dw_factor * aniso_dw_factor *  sin(2.0*PI*(h*r1+k*r2+l*r3)) * (struc_f_im + CURR_BULK_F_IM);
			    curr_df += 2.0 * CURR_DOMAIN_OCC * SCALE * (df_term1 + df_term2);
			};
		    };



		};	/* End atom loop */
		
	    };	/* End domain loop */
	    curr_j =  curr_df / (2.0 * CURR_F->sigma * CURR_F->value);
	  

	    /* Set Jacobian matrix element */
	    if (par_counter != 0) {
		gsl_matrix_set(J, f_counter, par_counter , curr_j);	
		    
	    };	
	};    /* End parameter loop */ 					
    };	/* End structure factor loop */
	
    return GSL_SUCCESS;
}


/****************************************************************************************************************
 * int ls_fit_fdf(const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J) 				*
 *														*
 * Calculates |F_calc^2-F_exp^2|) i.e f and d(|F_calc^2-F_exp^2|) / dx_i i.e. J 				*
 * (F: structure factors, x_i: free parameter) for use with GSL Levenberg-Marquardt least-squares fit functions	*
 *														*
 ****************************************************************************************************************/
int ls_fit_fdf(const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J) 
{
    ls_fit_f(x, params, f);
    ls_fit_df(x, params, J);
		
    return GSL_SUCCESS;
}

/********************************************************************************************************
 * void print_fit_result(int end_status, int end_iteration, int max_iteration,  			*
 *gsl_multifit_fdfsolver * s), size_t nr_paramter, size_t nr_struc_f, int xtodispl[], gsl_matrix *covar)     	*
 *													*
 * Print the results of the fit, i.e. number of iterations, reason for stopping the fit, 		*
 * atomic positions and their changes, covariant matrix.						*
 * 													*
 * end_status: Status code returned by GSL fitting functions.						*
 * end_iteration: Number of last iteration.								*
 * max_iterations: Maximal number of iterations before fit is stopped.					*
 * *s: Pointer to solver										* 
 * start_atoms: Input atomic positions.									*
 * end_atoms: Final atomic positions.									*
 * *covar: Matrix with covariances.									*
 * nr_parameter: Number of free parameters in fit.							*
 * nr_struc_f: the number of structure factors								*
 * xtodispl: Transformation array from vector x for fit to user supplied parameter numbers		*
 ********************************************************************************************************/
void print_fit_result(int end_status, unsigned long int end_iteration, unsigned long int max_iteration, gsl_multifit_fdfsolver * s, size_t nr_parameter, const struct structure_factor exp_struc_f[], size_t nr_struc_f,  int xtodispl[], gsl_matrix *covar, struct atom_structure end_positions[][MAX_ATOMS],struct atom_structure atoms_bulk[MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth, const size_t nr_atoms[], const double domain_occ[], const size_t nr_domains, size_t nr_atoms_bulk, enum output print_intermediate)
{
	
    double chi1;	/* chi = sqrt(sum f_i) (for intensities, from fit ls_fit_f) */
    double chi2_square;  /*  chi^2  = sum ((f_exp - |f_calc|)^2 / sigma_exp^2) (for structure factors */
    double r_factor;	/* R-factor R = (sum |f_exp - |f_calc||) / (sum f_exp) */
    double error;	/* Error computed from the square root of the diagonal elements of the covariance matrix */
    size_t counter;				
    size_t counter2;
    int par_nr;		/* number of the current parameter (as defined by the user) */
    double par_value;	/* Value of the current parameter */

	
    fprintf(stdout, "\nFit results:\n");

    if (end_status == GSL_SUCCESS && end_iteration < max_iteration) {	/* Fit converged ? */
	fprintf(stdout, "Fit converged after %ld iterations.\n", end_iteration);
    };

    if (end_iteration >= max_iteration) {					/* Fit not converged ? */
	fprintf(stdout, "Fit not converged after %ld iterations.\n", end_iteration);
	if (print_intermediate != short_output) {
	  fprintf(stdout, "Consider increasing the maximum number of iterations or find better starting values.\n");
	};
    };
	
    /*!!!!!!!!!!!!!!!!! IS CHI^2 /DOF RIGHT? !!!!!!!!!!!!!!!!!!!*/ 
    /* chi^2 (for intensities, from fit ls_fit_f) */
    chi1 = gsl_blas_dnrm2(s->f);   /* Norm of f => chi */
    fprintf(stdout, "chi^2 = %f, chi^2 / (degree of freedom) = %f (Intensities)\n", chi1*chi1, chi1*chi1/(nr_struc_f-nr_parameter));

    /* chi^2 (for structure factors ) */
    chi2_square = calc_chi_square(exp_struc_f, nr_struc_f, end_positions, atoms_bulk,lattice_par, penetration_depth, nr_atoms, domain_occ, nr_domains, nr_atoms_bulk);
    fprintf(stdout, "chi^2 = %f, chi^2 / (degree of freedom) = %f (Structure factors)\n", chi2_square , chi2_square/(nr_struc_f-nr_parameter));

    /* R-factor (for structure factors )*/
    r_factor = calc_r_factor(exp_struc_f, nr_struc_f, end_positions, atoms_bulk, lattice_par, penetration_depth, nr_atoms, domain_occ, nr_domains, nr_atoms_bulk);
    fprintf(stdout, "R = %f\n\n", r_factor); 

    /* Print scale factor and error */
    error = sqrt(gsl_matrix_get(covar, (size_t) 0, (size_t) 0));  /* Calculate errors from covariance matrix */
    par_value = gsl_vector_get(s->x, (size_t) 0);		/* Value of the parameter */
    fprintf (stdout, "Scale factor:   %.14f +/- %f\n",  par_value, error);   

/* Print parameters and errors */
    for (counter = 1; counter < nr_parameter; counter++) {
	par_nr = xtodispl[counter];				/* get user-defined number of current parameter */
	error = sqrt(gsl_matrix_get(covar, counter, counter));  /* Calculate errors from covariance matrix (no coupling among parameters) */
	par_value = gsl_vector_get(s->x, counter);		/* Value of the parameter */
	fprintf (stdout, "Parameter Nr. %d:   %f +/- %f\n", par_nr, par_value, error);
    };	
	
    /* Print covariance matrix */
    if (print_intermediate != short_output) {
      fprintf (stdout, "\nCovariance matrix:\n");
      fprintf (stdout, "    ");
      for (counter = 0; counter < nr_parameter; counter++) {	/* Print top line with parameter numbers */
	fprintf (stdout, "      %d      ", xtodispl[counter]);
      };
      fprintf (stdout, "\n");
      for (counter = 0; counter < nr_parameter; counter++) {	
	fprintf (stdout, " %d ", (int) counter);		/* Print parameter number */
	for (counter2 = 0; counter2 <nr_parameter; counter2++)	{
	  fprintf (stdout, " %.10f", gsl_matrix_get(covar, counter, counter2));   /* Print values of covariant matrix */
	};
	fprintf (stdout, "\n");
      };
    };

    fflush(stdout);  	/* Flush stdout at end of fit */

    return;
}



/****************************************************************************************************************
 * void fitted_coordinates(struct fit_all all, gsl_multifit_fdfsolver *solver_ptr, struct atom_struc *end_pos)	*
 *		 size_t * displ_x) 									*
 *														*
 * Calculate positions at end of fit and write them to 	end_pos.						*
 * 														*
 * all: Structure with data for structure factor calculation. 							*
 * *fit_par: Pointer to parameter vector x from fit
 * *solver_ptr: Pointer to GSL least-squares fit solver.							*
 * *end_pos: Write calculated positions here.									*
 * *displ_x: Pointer to array for translation from user supplied parameter number to vector x for fit. 		*
 ****************************************************************************************************************/

void fitted_coordinates(struct fit_all all, gsl_vector *fit_par, struct atom_structure (*end_pos_ptr)[MAX_ATOMS]) 
{
    size_t atom_counter;
    size_t domain_counter;
/*    double par; */
    double r1;					/* r1 coordinate of current atom */
    double r2;					/* r2 coordinate of current atom */
    double r3;					/* r3 coordinate of current atom */
    double dw;					/* Debye-Waller parameter of current atom */
    double aniso_dw[6];                         /* Components of anisotropic Debye-Waller factor (beta11,beta22,beta33,beta12,beta13,beta23) */
    double occ;					/* Occupancy of current atom */

    struct atom_structure (*atom_ptr)[MAX_ATOMS];
    struct displacements (*displ_ptr)[MAX_ATOMS];

    atom_ptr = all.atoms_ptr;			/* Pointer to array with atom data */
    displ_ptr = all.displ_ptr;			/* Pointer to array with displacements */

#define CURR_END_POS ((*(end_pos_ptr + domain_counter))[atom_counter])	
#define CURR_START_POS  ((*(atom_ptr + domain_counter))[atom_counter]) 
#define CURR_DISPL3  ((*(displ_ptr + domain_counter))[atom_counter]) 

    for (domain_counter = 0; domain_counter < all.nr_domains; domain_counter++) {
	for (atom_counter = 0; atom_counter < *((all.nr_atoms_ptr)+domain_counter); atom_counter++)	{

	    get_pos_dw_occ(fit_par, &all, &r1, &r2, &r3, &dw, &aniso_dw[0], &occ, atom_counter, domain_counter);

	    CURR_END_POS.x = r1;			/* Positions */
	    CURR_END_POS.y = r2;
	    CURR_END_POS.z = r3;
	    CURR_END_POS.dw_par = dw;	/* Debye-Waller parameter */
	    {
	      int ai;
	      for (ai=0; ai<6; ai++) {
		CURR_END_POS.aniso_dw[ai] = aniso_dw[ai];
	      };
	    };
	    CURR_END_POS.occ = occ;		/* Occupancy */
	    CURR_END_POS.el_nr =  CURR_START_POS.el_nr;		/* Atom element number */
	};
    };

    return;
}




/****************************************************************************************************************
 * void do_ls_fit(struct atom_structure atoms[][MAX_ATOMS], const struct atom_structure atoms_bulk[MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth, struct displacements displ[][MAX_ATOMS],            *
 * int fixed_par[], double start_par[], double domain_occ[], size_t nr_atoms[], size_t nr_domains,                *
 * struct structure_factor exp_struc_f[], size_t nr_struc_f, unsigned int max_iteration, double delta_abs,      *
 * double delta_rel, struct atom_structure end_positions[][MAX_ATOMS],  unsigned long int rng_seed,             *
 * enum output print_intermediate, int n_fit, double par_var_width)                                             * 
 * 														*
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
 * start_par[nr_param]: Starting values for the fit parameters. 						*
 * nr_atoms: Number of atoms in the unit cell.									*
 * nr_atoms_bulk: Number of atoms in the bulk structure                                                 *
 * exp_struc_f[f nr]: Structure factor h,k,l of structure factor (f nr).					*
 * nr_struc_f: Number of structure factors.									*
 * max_iteration: Maximum number of iterations.									*
 * delta_abs: Absolute convergence criterion 									*
 * delta_rel: Relative convergence criterion (not very useful)							*
 * *end_positions: Pointer to array where final positions shoutld be stored.					*
 ****************************************************************************************************************/

void do_ls_fit(struct atom_structure atoms[][MAX_ATOMS], struct atom_structure atoms_bulk[], const struct latt_par_structure lattice_par, double penetration_depth, struct displacements displ[][MAX_ATOMS], int fixed_par[], double start_par[], double domain_occ[], size_t nr_atoms[], size_t  nr_domains, size_t nr_atoms_bulk, struct structure_factor exp_struc_f[], size_t nr_struc_f, unsigned long int max_iteration, double delta_abs, double delta_rel, struct atom_structure end_positions[][MAX_ATOMS],  unsigned long int rng_seed, enum output print_intermediate, long unsigned int n_fit, double par_var_width)
{
    gsl_multifit_function_fdf fdf;	/* Structure for function for calculating F and derivatives + other data. */
    /* int (* f): Pointer to function calculating |F_calc^2-F_exo^2|.
     * int (* df): Pointer to function calculating derivatives.X
     * int (* fdf): Pointer to function calculating |F_calc^2-F_exo^2| and derivatives.
     * size_t n: Number of structure factors.
     * size_t p: Number of free parameters for fit.
     * void * params: a pointer to the parameters of the function. */
					
    double bulk_structure_factor_real[MAX_STRUCTURE_FACTORS];   /* array with real part of bulk structurefactors */
    double bulk_structure_factor_im[MAX_STRUCTURE_FACTORS];   /* array with imaginary part of bulk structurefactors */

    size_t nr_param;			/* Number of free parameters */
    size_t par_counter;	
	
    size_t displtox[MAX_FREE_PARAMETERS];		/* Translation from user supplied parameter nrs fit vector x */
    int xtodispl[MAX_FREE_PARAMETERS];		/* Translation from fit vector x to user supplied parameter nrs */

    const gsl_rng_type * rng_type;		/* Random number generator type */
    gsl_rng * rng_ptr;				/* Random number generator pointer */	

    unsigned int iteration_counter = 0;		
    unsigned int f_counter;			/* Counter over structure factors */
    int status = GSL_CONTINUE;			/* Status of the fit */

    static struct fit_all all;	    /* Structure with atoms, displacements, structure factors */

    unsigned int fit_counter;
	
    gsl_multifit_fdfsolver * solver_ptr;   	/* Pointer to allocated solver */
    const gsl_multifit_fdfsolver_type * lmder_solver_ptr;  	/* Pointer to Levenberg-Marquardt solver function */
	
    gsl_vector *displ_init;			/* Initial displacements (set to 0) */
    gsl_matrix *covar;				/* Matrix for covariances */
    gsl_matrix *jacobian;                /* Matrix for jacobian */


    /* Setup random number generator */
    gsl_rng_env_setup();		
    rng_type  = gsl_rng_default;
    rng_ptr  = gsl_rng_alloc(rng_type);
    gsl_rng_set(rng_ptr, rng_seed);

    /* Find the number of parameters in which the F should be minimized and assign vector x for fit */
    nr_param = find_nr_param(displ, fixed_par, nr_atoms, nr_domains, displtox, xtodispl) + 1;   

    displ_init = gsl_vector_calloc(nr_param);		/* Vector with fit parameters */
    covar = gsl_matrix_alloc (nr_param, nr_param);	/* Matrix for covariances */
    jacobian = gsl_matrix_alloc (nr_struc_f, nr_param);	/* Matrix for Jacobian */

    /* Calculate bulk structure factors */
    {
      
      for (f_counter = 0; f_counter < nr_struc_f; f_counter++)
	{ 	
	calc_structure_factor_bulk(atoms_bulk, nr_atoms_bulk, lattice_par, penetration_depth, exp_struc_f[f_counter].h, exp_struc_f[f_counter].k, exp_struc_f[f_counter].l, &(bulk_structure_factor_real[f_counter]), &(bulk_structure_factor_im[f_counter]));
	};
    };


    /* Loop over the number of fits */
    for (fit_counter = 0; fit_counter < n_fit; fit_counter++) 
	{
	    iteration_counter = 0;
	    status = GSL_CONTINUE;

	    /* Starting values for parameters */
	    {
		double start_par_save;
		
		fprintf (stdout, "\n-----------------------------\n Fit Nr. %d", fit_counter + 1);
		if (print_intermediate != short_output) {
		  fprintf	(stdout, "\nStarting parameters:\n");
		  /* Scale parameter */
		  fprintf (stdout, "Parameter Nr. 0: %f\n", start_par[xtodispl[0]]);
		};
		gsl_vector_set(displ_init, (size_t) 0, start_par[xtodispl[0]]);
		for (par_counter = 1; par_counter < nr_param; par_counter++) {
		    /* Change starting parameters according to start_par_width) */
		    start_par_save = start_par[xtodispl[par_counter]] + gsl_rng_uniform(rng_ptr) * 2.0 * par_var_width - par_var_width;
		    /* Print parameter */
		    if (print_intermediate != short_output) {
 
		      fprintf (stdout, "Parameter Nr. %d: %f\n", xtodispl[par_counter], start_par_save);
		    };
		    gsl_vector_set(displ_init, par_counter,  start_par_save);
		};
		fprintf (stdout, "\n");
	    };
	    
	    /* Starting values for parameters */

	    /*	    for (par_counter = 0; par_counter < nr_param; par_counter++) {
		gsl_vector_set(displ_init, par_counter, start_par[xtodispl[par_counter]]);
		}; */

	    /* Print fixed parameters */
	    if (print_intermediate != short_output) {  
	      fprintf (stdout, "Fixed parameters:\n");
	      for (par_counter = 1; par_counter < MAX_FREE_PARAMETERS; par_counter++) {
		if (fixed_par[par_counter] == TRUE) {
		    
		  fprintf (stdout, "Parameter Nr. %d fixed to: %f\n", (int) par_counter, start_par[par_counter]);
		};
	      };
	      fprintf (stdout, "\n");
	    };

	    /* Assign structure with atomic data, structure factors, ... */
	    all.atoms_ptr = &atoms[0];
	    all.atoms_bulk_ptr = &atoms_bulk[0];
	    all.displ_ptr = &displ[0];
	    all.struc_f_ptr =  exp_struc_f;
	    all.domain_occ_ptr = domain_occ;
	    all.lattice_par = lattice_par;
	    all.bulk_structure_factor_real_ptr = &bulk_structure_factor_real[0];
	    all.bulk_structure_factor_im_ptr = &bulk_structure_factor_im[0];
	    all.penetration_depth = penetration_depth;
	    all.start_par_ptr = &start_par[0];
	    all.fixed_par_ptr = &fixed_par[0];
	    all.displtox = displtox;
	    all.nr_atoms_ptr = &nr_atoms[0];
	    all.nr_atoms_bulk = nr_atoms_bulk;
	    all.nr_domains = nr_domains;
	    all.nr_struc_f = nr_struc_f;
	    all.nr_parameters = nr_param;


	    /* Setup structure with functions to minimize */
	    fdf.f = &ls_fit_f;
	    fdf.df = &ls_fit_df;
	    fdf.fdf = &ls_fit_fdf;
	    fdf.n = nr_struc_f;
	    fdf.p = nr_param;
	    fdf.params = &all;
		
	    /* Initialize the solver for the functions defined in fdf to Levenberg-Marquardt with 
	     * nr_struc_f data points, nr_param free parameters and initial values for the free parameters 
	     * given by displ_init. */
	    lmder_solver_ptr = gsl_multifit_fdfsolver_lmder;
	    solver_ptr  = gsl_multifit_fdfsolver_alloc (lmder_solver_ptr, nr_struc_f, nr_param);
	    gsl_multifit_fdfsolver_set (solver_ptr, &fdf, displ_init);
    
	    /* Print current parameters */
	    if (print_intermediate == intermediate_status) {
		size_t counter;
		int curr_par_nr;
		double curr_par;
		double curr_dx;
		fprintf (stdout, "\nStart, status %d, f: %f\n", status, gsl_blas_dnrm2(solver_ptr->f));
		fprintf (stdout, "Current parameters:\n");
		
		curr_par = gsl_vector_get (gsl_multifit_fdfsolver_position(solver_ptr), (size_t) 0);	/* Get scale */
		curr_dx = gsl_vector_get (solver_ptr->dx, (size_t) 0);	
		fprintf (stdout, "Scale: %.14f, change: %.14f\n", curr_par, curr_dx);
	
		for (counter = 1; counter < nr_param; counter++) {
		    curr_par_nr = xtodispl[counter];	   /* get user-defined number of current parameter */
		    curr_par = gsl_vector_get (gsl_multifit_fdfsolver_position(solver_ptr), counter);
		    curr_dx = gsl_vector_get (solver_ptr->dx, counter);
		    fprintf(stdout, "Parameter nr. %d: %f, change %f\n", curr_par_nr, curr_par, curr_dx);
		};
	    };	

	    /* Do iterations */
	    while ((status == GSL_CONTINUE) && (iteration_counter < max_iteration))	{	
		/* Iterate solver one time */
		status = gsl_multifit_fdfsolver_iterate (solver_ptr);	
		
		if (status != GSL_SUCCESS && status != GSL_CONTINUE) {			/* Problem ? */
		    fprintf (stderr, "Problem: Iteration %d returned %d \n", iteration_counter+1, status);
		    fprintf (stdout, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		    fprintf (stdout, "!!! Problem: Iteration returned %d !!!\n", status);
		    fprintf (stdout, "!!! %s !!!\n", gsl_strerror (status));
		    fprintf (stdout, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		};	


		/* Print current parameters */
		if (print_intermediate == intermediate_status) {
		    size_t counter;
		    int curr_par_nr;
		    double curr_par;
		    double curr_dx;
		    fprintf (stdout, "\nIteration %d, status %d. f: %f\n", iteration_counter+1, status, gsl_blas_dnrm2(solver_ptr->f));
		    fprintf (stdout, "Current parameters:\n");
		    
		    curr_par = gsl_vector_get (gsl_multifit_fdfsolver_position(solver_ptr), (size_t) 0);	/* Get scale */
		    curr_dx = gsl_vector_get (solver_ptr->dx, (size_t) 0);	
		    fprintf (stdout, "Scale: %.14f, change: %.14f\n", curr_par, curr_dx);
	
		    for (counter = 1; counter < nr_param; counter++) {
			curr_par_nr = xtodispl[counter];	   /* get user-defined number of current parameter */
			curr_par = gsl_vector_get (gsl_multifit_fdfsolver_position(solver_ptr), counter);
			curr_dx = gsl_vector_get (solver_ptr->dx, counter);
			fprintf(stdout, "Parameter nr. %d: %f, change: %f\n", curr_par_nr, curr_par, curr_dx);
		    };
		};
	
	
	    
		/* Check for convergence 
		 * Stop if the change dx_i in each free parameter x_i is smaller than 
		 * delta_abs + delta_rel * |x_i| */
		status = gsl_multifit_test_delta (solver_ptr->dx, solver_ptr->x, delta_abs, delta_rel);
		if (print_intermediate == intermediate_status) {
		    fprintf (stdout, "Convergence status: %d\n",  status);
		};
		iteration_counter++;
	    };			
	
	    /* Calculate final coordinates */
	    fitted_coordinates(all, solver_ptr->x, end_positions);
	
	    /* Calculate covariance matrix */
         gsl_multifit_fdfsolver_jac (solver_ptr, jacobian);
	    gsl_multifit_covar (jacobian, 0.0, covar);
	
	    /* Print state at end of fit */
	    print_fit_result(status, iteration_counter, max_iteration,  solver_ptr, nr_param, exp_struc_f, nr_struc_f, xtodispl, covar, end_positions,  atoms_bulk, lattice_par, penetration_depth,  nr_atoms, domain_occ, nr_domains,nr_atoms_bulk, print_intermediate );
	    
	    /* Free memory of solver */
	    gsl_multifit_fdfsolver_free (solver_ptr);
	};

    /* Free memory of covariance matrix */
    gsl_matrix_free(covar);
    /* Free memory of initial displacements */
    gsl_vector_free(displ_init);
	    
    return;
}

/****************************************************************************************************************
 * double sa_chi_square(void *params) 										*
 * 														*
 * Calculate  the energy of a configuration i.e. chi^2 for parameter values.					*
 * 														*
 * params: 
 ****************************************************************************************************************/

double sa_chi_square(void *params) 
{
    double chi_square = 0.0;		/* chi^2 = f^2, i.e. the value to be minimized by the fit */
    gsl_vector (*fit_par_ptr);	/* Pointer to fit parameters */

    double scattering_length_sq;			/* The length of the scattering vector (used for
							   calculating the atomic scattering factor) */
    size_t f_counter;					/* Counter over structure factors */
    size_t domain_counter;					/* Counter over surface domains */
    static double struc_f_square_calc[MAX_STRUCTURE_FACTORS]; /* calculated structure factors */
    double struc_f_real;	/* Real part of the structure factor */
    double struc_f_im;	/* Imaginary part of the structure factor */
    double struc_f_bulk_real = 0.0;   			/* Real part of the bulk structure factor */
    double struc_f_bulk_im = 0.0;   				/* Complex part of the bulk structure factor */    double residual = 0.0;				/* (F_calc^2-F_exp^2)/(2*sigma*F) */
    double sum1 = 0.0;			/* Summing for scale */
    double sum2 = 0.0;			
    double scale;				/* Scale factor between intensities!*/
    double h;					/* index h of current structure factor */
    double k; 					/* index k of current structure factor */
    double l;					/* index l of current structure factor */
    double (*domain_occ_ptr);				/* Pointer to array with domain occupancies */
    struct atom_structure (*atom_bulk_ptr);	/* Pointer to array with bulk atomic data */

#define PAR ((struct fit_all *) params)
#define CURR_F ((struct structure_factor *)((((struct fit_all *) params)->struc_f_ptr) + f_counter))
#define CURR_DOMAIN_OCC *(domain_occ_ptr + domain_counter)

    fit_par_ptr = ((struct fit_all *) params)->fit_par_ptr;
    domain_occ_ptr = ((struct fit_all *) params)->domain_occ_ptr; 
    atom_bulk_ptr = ((struct fit_all *) params)->atoms_bulk_ptr;

    for (f_counter = 0; f_counter <  PAR->nr_struc_f; f_counter++) {	/* Loop over all structure factors */

	/* Initialize struc_fs */
       	struc_f_square_calc[f_counter] = 0.0; 

	/* Get h, k, l */
	h = ((struct structure_factor *)((((struct fit_all *) params)->struc_f_ptr) + f_counter))->h;
	k = ((struct structure_factor *)((((struct fit_all *) params)->struc_f_ptr) + f_counter))->k;
	l = ((struct structure_factor *)((((struct fit_all *) params)->struc_f_ptr) + f_counter))->l;

	/* Square of the scattering length  */
/*	scattering_length_sq = (4.0/3.0/(PAR->lattice_par.a * PAR->lattice_par.a) * (h*h + k*k + h*k) + l*l / (PAR->lattice_par.c * PAR->lattice_par.c)) / 4.0; */
   	scattering_length_sq = calc_scattering_length_sq(PAR->lattice_par, h, k, l);

	/* Calculate bulk amplitude */
	calc_structure_factor_bulk(atom_bulk_ptr, PAR->nr_atoms_bulk, PAR->lattice_par, PAR->penetration_depth, h, k, l, &struc_f_bulk_real, &struc_f_bulk_im);

        /* Loop over domains and add the structure factor incoherently (i.e. absolute values)*/
	for(domain_counter = 0; domain_counter < PAR->nr_domains; domain_counter++)   {

	    struc_f_real = 0.0;
	    struc_f_im = 0.0;

	    calc_real_im_struc_f(fit_par_ptr, params, h, k, l, scattering_length_sq,  domain_counter, &struc_f_real, &struc_f_im);
	    struc_f_real = struc_f_real + struc_f_bulk_real;
	    struc_f_im = struc_f_im + struc_f_bulk_im;
	    struc_f_square_calc[f_counter] += CURR_DOMAIN_OCC * (struc_f_real*struc_f_real + struc_f_im*struc_f_im);
	    
	};  /* End domain loop */
	
/* This is for scale determination */
	sum1 +=  CURR_F->value* CURR_F->value * struc_f_square_calc[f_counter] / (2.0 * CURR_F->sigma * CURR_F->value * 2.0 * CURR_F->sigma * CURR_F->value);
	sum2 +=  struc_f_square_calc[f_counter]*struc_f_square_calc[f_counter] / (2.0 * CURR_F->sigma * CURR_F->value * 2.0 * CURR_F->sigma * CURR_F->value);

    }; /* End structure factor loop */
 
    /* Average scale (weighted with error) for intensities!*/
    scale = sum1 / sum2;

    /* Loop over structure factors 2 to calculate chi_square */
    for (f_counter = 0; f_counter < PAR->nr_struc_f; f_counter++) {	/* Loop over all structure factors */
	
    residual = (scale * struc_f_square_calc[f_counter] - CURR_F->value*CURR_F->value) / (2.0 * CURR_F->sigma * CURR_F->value); 
    chi_square += residual * residual;

    };	/* End structure factor loop 2 */

/*    fprintf (stdout, "Chi^2: %f, scale: %f\n", chi_square, scale);  */

    return (chi_square);
}

/****************************************************************************************************************
 * void sa_step (const gsl_rng *rng_ptr, void *params, double step_size)					*
 * 														*
 * Change the current values of the fit parameter for a new step. Random step  with random number generator	*
 * rng_ptr and maximum step size step_size									*
 * 														*
 * !!!! At the moment Debye-Waller parameter and positions have same step width !!!!				*
 *      Better to change that ?											*
 * !!!! Maximum value of parameter limited to 1 								*
 *      Therefore also Debye-Waller parameter limited to 1 -> Change scale for DW parameter in fit file !	*
 ****************************************************************************************************************/


void sa_step (const gsl_rng *rng_ptr, void *params, double step_size)
{
    gsl_vector (*fit_par_ptr);	/* Pointer to fit parameters */
    double old_par;		/* Old value of the current fit parameter */
    double new_par;		/* New value of the current fit parameter */
    size_t par_counter;		/* Count over fit parameters */
    double random_step;		/* Random step for the current parameter */



    fit_par_ptr = ((struct fit_all *) params)->fit_par_ptr;

    /* Loop over all fit parameters and change their value */
    for (par_counter = 1; par_counter <  ((struct fit_all *) params)->nr_parameters; par_counter++) 
	{
/*	    fprintf(stderr, "sa_step, %d\n", par_counter); */
	    new_par = 0.0;
	    old_par = gsl_vector_get(fit_par_ptr, par_counter);		/* Current value of fit parameter */

	    do {		/* Find new value for fit parameter, but must be <=1.0 */
		random_step = gsl_rng_uniform(rng_ptr) * 2.0 * step_size - step_size ;	/* Random step for fit parameter */
		new_par = old_par + random_step;				/* New value of fit parameter */
	    } while (new_par > 1.0 || new_par <-1.0);

	    gsl_vector_set(fit_par_ptr, par_counter, new_par);
	};

    return;
}


/****************************************************************************************************************
 * double sa_distance (void *params1, void *params2)								*
 * 														*
 * Calculate the "distance" between params1->fit_par and params2->fit_par.					*
 * 														*
 * !!!! At the moment just take maximum absolute difference between the vector elements of fit_par  !!!!	*
 * 	Better distance ?											*
 ****************************************************************************************************************/

double sa_distance (void *params1, void *params2)
{
    double distance;		/* "Distance" between params1->fit_par and params2->fit_par. */ 
    double min_diff;		/* Minimum element of params1->fit_par and params2->fit_par. */
    double max_diff;		/* Maximum element of params1->fit_par - params2->fit_par. */
    gsl_vector * difference; 	/* params1->fit_par - params2->fit_par. */
    gsl_vector (*fit_par_ptr1);	/* Pointer to fit parameters 1 */
    gsl_vector (*fit_par_ptr2);	/* Pointer to fit parameters 2 */

/* fprintf(stderr, "sa_distance\n"); */


    fit_par_ptr1 = ((struct fit_all *) params1)->fit_par_ptr;
    fit_par_ptr2 = ((struct fit_all *) params2)->fit_par_ptr;

    /* Allocate difference */
    difference = gsl_vector_alloc(((struct fit_all *) params1)->nr_parameters);

    /* Subtract  params1->fit_par from params2->fit_par. */
    gsl_vector_memcpy(difference, fit_par_ptr1);
    gsl_vector_sub(difference, fit_par_ptr2);
	
    /* Find the maximum absolute difference */
    gsl_vector_minmax(difference, &min_diff, &max_diff);
    if (fabs(min_diff) > max_diff) {
	distance = min_diff;
    }
    else {
	distance = max_diff;
    };
    

   return distance;
}

/****************************************************************************************************************
 * void print_sa_step (void *params);										*
 * 														*
 * Print current fit parameters for step 									*
 ****************************************************************************************************************/

void sa_print_step (void *params)
{
    gsl_vector (*fit_par_ptr);	/* Pointer to fit parameters */

    fit_par_ptr = ((struct fit_all *) params)->fit_par_ptr;
    gsl_vector_fprintf(stdout, fit_par_ptr, "%f");
    fprintf (stdout, "\n");
    return;
}


/****************************************************************************************************************
 * void sa_copy(void *params_source, void *params_dest)								*
 * 														*
 * Copy the parameter structure params from *params_source to *params_dest for use in GSL simulated 		*
 * annealing function							 					*
 ****************************************************************************************************************/

void sa_copy(void *params_source, void *params_dest)
{
    gsl_vector *fit_par_source;			/* Values of the parameters for the fit (source) */
    gsl_vector *fit_par_dest;			/* Values of the parameters for the fit (destination); */

/* fprintf(stderr, "sa_copy\n"); */

    /* Get the fit parameter vector in destination, otherwise it is overwritten */
    fit_par_dest = ((struct fit_all *) params_dest)->fit_par_ptr;

    /* Copy params */
    memcpy(params_dest, params_source, sizeof(struct fit_all));

    /* Copy the fit parameters	*/
    fit_par_source = ((struct fit_all *) params_source)->fit_par_ptr;
    gsl_vector_memcpy(fit_par_dest, fit_par_source);
    ((struct fit_all *) params_dest)->fit_par_ptr = fit_par_dest;

    return;
}

/****************************************************************************************************************
 * void sa_copy_construct(void *params_source)									*
 *														*
 * Create a new copy of the parameter structure *params_source for use in GSL simulated 			*
 * annealing function												*
 ****************************************************************************************************************/

void * sa_copy_construct(void *params_source)
{
    gsl_vector *fit_par_source;			/* Values of the parameters for the fit (source) */
    gsl_vector *fit_par_dest;			/* Values of the parameters for the fit (destination); */
    struct fit_all *params_dest;		/* New copy of all parameters */

/* fprintf(stderr, "sa_copy_cosntruct\n");    */
    
    /* Copy params */
    params_dest = malloc(sizeof(struct fit_all));
    memcpy(params_dest, params_source, sizeof(struct fit_all));

    
   /* Copy the fit parameters	*/
    fit_par_dest = gsl_vector_calloc(((struct fit_all *) params_source)->nr_parameters);
    fit_par_source = ((struct fit_all *) params_source)->fit_par_ptr;
    gsl_vector_memcpy(fit_par_dest, fit_par_source);
    ((struct fit_all *) params_dest)->fit_par_ptr = fit_par_dest;

    return params_dest;
}

/************************************************************************************************
 * void sa_destroy(void *params) 								*
 *												*
 * Destroys the parameter structure *params for use in GSL simulated 				*
 * annealing function.										*
 ************************************************************************************************/
void sa_destroy(void *params) 
{
    gsl_vector *fit_par;			/* Values of the parameters for the fit */

/* fprintf(stderr, "sa_destroy\n"); */

    /* Free vector with parameters */
    fit_par = ((struct fit_all *) params)->fit_par_ptr;
    gsl_vector_free(fit_par);

    /* Free params */
    free(params);

    return;
}



/****************************************************************************************************************
 * void print_sa_result (void *params, int xtodispl[])								*
 * 														*
 * Print results of the simulated annealing fit.								*
 ****************************************************************************************************************/

void print_sa_result (void *params, int xtodispl[], struct atom_structure end_positions[][MAX_ATOMS])
{
    double chi1_square;	/* chi^2 (for intensities) */
    double chi2_square;	/* chi^2  = sum ((f_exp - |f_calc|)^2 / sigma_exp^2) (for structure factors) */
    double r_factor;	/* R-factor R = (sum |f_exp - |f_calc||) / (sum f_exp) */
   
    size_t counter;				
    int par_nr;		/* number of the current parameter (as defined by the user) */
    double par_value;	/* Value of the current parameter */
    double (*domain_occ_ptr);				/* Pointer to array with domain occupancies */
     struct structure_factor (*exp_struc_f_ptr);	/* Pointer to experimental structure factors */

    domain_occ_ptr = ((struct fit_all *) params)->domain_occ_ptr; 
    exp_struc_f_ptr = ((struct fit_all *) params)->struc_f_ptr;

#define PAR ((struct fit_all *) params)

    fprintf(stdout, "\n Results of simulated annealing run:\n");

/* chi^2 (for intensities) */
    chi1_square = sa_chi_square(params);
    fprintf(stdout, "chi^2 = %f, chi^2 / (degree of freedom) = %f (Intensities) \n", chi1_square, chi1_square/(PAR->nr_struc_f-PAR->nr_parameters));

    /* chi^2 (for structure factors ) */
    chi2_square = calc_chi_square(exp_struc_f_ptr, PAR->nr_struc_f, end_positions, PAR->atoms_bulk_ptr, PAR->lattice_par, PAR->penetration_depth, PAR->nr_atoms_ptr, domain_occ_ptr, PAR->nr_domains, PAR->nr_atoms_bulk); 
    fprintf(stdout, "chi^2 = %f, chi^2 / (degree of freedom) = %f (Structure factors)\n", chi2_square , chi2_square/(PAR->nr_struc_f-PAR->nr_parameters));

    /* R-factor (for structure factors )*/
    r_factor = calc_r_factor(exp_struc_f_ptr, PAR->nr_struc_f, end_positions, PAR->atoms_bulk_ptr, PAR->lattice_par, PAR->penetration_depth, PAR->nr_atoms_ptr, domain_occ_ptr, PAR->nr_domains, PAR->nr_atoms_bulk); 
    fprintf(stdout, "R = %f\n\n", r_factor); 

    /* Print parameters */
    for (counter = 1; counter < PAR->nr_parameters; counter++)	{
	par_nr = xtodispl[counter];				/* get user-defined number of current parameter */
	par_value = gsl_vector_get(PAR->fit_par_ptr, counter);		/* Value of the parameter */
	fprintf (stdout, "Parameter Nr. %d:   %f \n", par_nr, par_value);
    };	

    fflush(stdout);  	/* Flush stdout at end of fit */

    return;
}


/****************************************************************************************************************
 * void do_sa_fit(struct atom_structure atoms[][],  struct atom_structure atoms_bulk[MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth, struct displacements displ[][], double start_par[],		*
 * double domain_occ[], size_t nr_atoms[], size_t nr_domains, struct structure_factor exp_struc_f[],	*
 * size_t nr_struc_f, gsl_siman_params_t sa_par, struct atom_structure (*end_positions)[MAX_ATOMS],	*
 *  unsigned long int rng_seed, enum output  print_intermediate)						*
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
 * penetration_depth: Penetration depth of X-rays into crystal (reciprocal units of lattice parameter)  *
 * displ[domain nr][atom nr]: Displacements for atom (atom nr) in domain (domain nr).				*
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
 * print_intermediate: Should intermediate information be printed during fit?
 ****************************************************************************************************************/

void do_sa_fit(struct atom_structure atoms[][MAX_ATOMS], struct atom_structure atoms_bulk[MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth,  struct displacements displ[][MAX_ATOMS], int fixed_par[], double start_par[], double domain_occ[], size_t nr_atoms[], size_t nr_domains, size_t nr_atoms_bulk, struct structure_factor exp_struc_f[], size_t nr_struc_f, gsl_siman_params_t sa_par, struct atom_structure (*end_positions)[MAX_ATOMS], unsigned long int rng_seed, enum output  print_intermediate, unsigned long int n_fit, double par_var_width)
{

    size_t nr_param;			/* Number of free parameters */
    size_t par_counter;	
	
    size_t displtox[MAX_FREE_PARAMETERS];		/* Translation from user supplied parameter nrs fit vector x */
    int xtodispl[MAX_FREE_PARAMETERS];		/* Translation from fit vector x to user supplied parameter nrs */
    static struct fit_all all;	    /* Structure with atoms, displacements, structure factors */

    gsl_vector *fit_par;			/* Values of the parameters for the fit; */
    const gsl_rng_type * rng_type;		/* Random number generator type */
    gsl_rng * rng_ptr;				/* Random number generator pointer */
    
    unsigned int fit_counter ;

    /* Setup random number generator */
    gsl_rng_env_setup();		
    rng_type  = gsl_rng_default;
    rng_ptr  = gsl_rng_alloc(rng_type);
    gsl_rng_set(rng_ptr, rng_seed);

    /* Find the number of parameters in which the F should be minimalized and assign vector x for fit */
    nr_param = find_nr_param(displ, fixed_par, nr_atoms, nr_domains, displtox, xtodispl) + 1;   

    fit_par = gsl_vector_calloc(nr_param);	/* Vector with fit parameters */

    /* Loop over the number of fits */
    for (fit_counter = 0; fit_counter < n_fit; fit_counter++) 
	{
  
	    /* Starting values for parameters */
	    {
		double start_par_save;
		
		fprintf (stdout, "\n-----------------------------\n Fit Nr. %d", fit_counter + 1);
		fprintf	(stdout, "\nStarting parameters:\n");
		for (par_counter = 0; par_counter < nr_param; par_counter++) {
		    /* Change starting parameters according to start_par_width) */
		    start_par_save = start_par[xtodispl[par_counter]] + gsl_rng_uniform(rng_ptr) * 2.0 * par_var_width - par_var_width;
		    /* Print parameter */
		    fprintf (stdout, "Parameter Nr. %d: %f\n", xtodispl[par_counter], start_par_save);
		    gsl_vector_set(fit_par, par_counter,  start_par_save);
		};
		fprintf (stdout, "\n");
	    };

	    /* Print fixed parameters */
	    fprintf (stdout, "Fixed parameters:\n");
	    for (par_counter = 1; par_counter < MAX_FREE_PARAMETERS; par_counter++) {
		if (fixed_par[par_counter] == TRUE) {
		    
		  fprintf (stdout, "Parameter Nr. %d fixed to: %f\n", (int) par_counter, start_par[par_counter]);
		};
	    };
	    fprintf (stdout, "\n");


	    /* Assign structure with atomic data, structure factors, ... */
	    all.atoms_ptr = &atoms[0];
	    all.atoms_bulk_ptr = &atoms_bulk[0];
	    all.displ_ptr = &displ[0];
	    all.struc_f_ptr =  exp_struc_f;
	    all.domain_occ_ptr = domain_occ;
	    all.lattice_par = lattice_par;
	    all.penetration_depth = penetration_depth;
	    all.start_par_ptr = &start_par[0];
	    all.fixed_par_ptr = &fixed_par[0];
	    all.displtox = displtox;
	    all.nr_atoms_ptr = nr_atoms;
	    all.nr_atoms_bulk = nr_atoms_bulk;
	    all.nr_domains = nr_domains;
	    all.nr_struc_f = nr_struc_f;
	    all.nr_parameters = nr_param;
	    all.fit_par_ptr = fit_par;

	    /* Do simulated annealing run */
	    switch (print_intermediate) {    /* Should intermediate information at each temperature step be printed? */
	    case no_intermediate_status:
	      gsl_siman_solve(rng_ptr, &all, sa_chi_square, sa_step, sa_distance, NULL, sa_copy, sa_copy_construct, sa_destroy, (size_t) 0, sa_par); 	
		break;
	    case short_output:
	      gsl_siman_solve(rng_ptr, &all, sa_chi_square, sa_step, sa_distance, NULL, sa_copy, sa_copy_construct, sa_destroy, (size_t) 0, sa_par); 	
		break;
	    case intermediate_status:
	      gsl_siman_solve(rng_ptr, &all, sa_chi_square, sa_step, sa_distance, sa_print_step, sa_copy, sa_copy_construct, sa_destroy, (size_t) 0, sa_par); 
		break;
	    default:	/* Should never get here */
		fprintf(stderr, "Error: print_intermediate has strange value!");
		exit(8);
		break;
	    };
	    
    

	    /* Calculate final coordinates */
	    fitted_coordinates(all, fit_par, end_positions);


	    /* Print results */
	    print_sa_result (&all, xtodispl, end_positions);



	};	

    gsl_vector_free(fit_par);		/* Free vector for fit parameters */
    gsl_rng_free(rng_ptr);		/* Free random number generator */
    return;
}


/********************************************************************************************************
 * double calc_electron_diff_map(double r1, double r2, struct structure_factor struc_f[], size_T nr_f, 	*
 *                             const struct atom_structure  atoms[MAX_ATOMS],  size_t nr_atoms[], double scale)*
 *                             										*
 * Calculates and returns the value of the electron difference map at r1, r2. 				*
 *													*
 *
 * !!!!!!!!!!!! CHECK  SCALE => REALLY IN UNITS OF e????????????????????
 *													
 * !!!!!!!!!!!! ONLY RELIABLE FOR ONE DOMAIN !!!!!!!!!!!!!!!!!!!!!!!!!!
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
 double calc_electron_diff_map(double r1, double r2, struct structure_factor exp_struc_f[], size_t nr_f,  struct atom_structure  atoms[][MAX_ATOMS],const struct latt_par_structure lattice_par, double penetration_depth,  size_t nr_atoms[],  size_t nr_domains, double domain_occ[], double scale) 
/* double calc_electron_diff_map(double r1, double r2, struct structure_factor exp_struc_f[], size_t nr_f, const struct atom_structure  atoms[MAX_ATOMS], size_t nr_atoms,  double scale) */
{
    double calc_struc_f_real = 0.0;   			/* Real part of the structure factor for the model */
    double calc_struc_f_im = 0.0;   			/* Complex part of the structure factor for the model */
    double calc_struc_f_abs[MAX_DOMAINS];		/* Absolute value of the structure factor */ 
    double calc_struc_f_phase[MAX_DOMAINS];	/* Phase of the structure factor */ 
    double calc_struc_f_square_sum = 0.0;		/* Sum of the intensity over all domains  */ 
    double domain_frac;					/* Approximate fraction of experimental structure factor corresponding to current domain */
    size_t atom_counter;  			/* Atom number used in loop */
    size_t domain_counter;  			/* Domain number used in loop */
    
    double dw_factor;   			        /* debye-waller factor: exp(-1*dw_par*Q^2) */
    double atom_scattering_f = 0.0;			/* Atomic scattering factor */
    double scattering_length_sq;			/* The length of the scattering vector (used for
							   calculating the atomic scattering factor) */

    size_t matrix_counter;			/* Counter used for the symmetry matrices */
    double electron_diff_value = 0.0;			/* Value of the electron difference map at r1,r2 */
    double unit_cell_area;				/* Area of the surface unit cell (Ang^2) */
    size_t f_counter;				/* Counter used to sum all structure factors */
    double h;                           		/* Current h */
    double k; 	        		                /* Current k */			
    double l;   	                        	/* Current l */			
    double value = 0.0;                    		/* Current structure factor value */


/* test for creating symmetry related reflections 
FILE *testfile;
testfile = fopen("test.dat","w");		*/


    /* Loop over all reflections */
    for (f_counter = 0; f_counter < nr_f; ++f_counter) {



	/* generate symmetrical equivalent reflections */
	for (matrix_counter = 0; matrix_counter < nr_symmetry_matrices; ++matrix_counter) {
	    h = exp_struc_f[f_counter].h * symmetry_matrix[matrix_counter][0][0] + exp_struc_f[f_counter].k * symmetry_matrix[matrix_counter][0][1];  
	    k = exp_struc_f[f_counter].h * symmetry_matrix[matrix_counter][1][0] + exp_struc_f[f_counter].k * symmetry_matrix[matrix_counter][1][1];
	    l = exp_struc_f[f_counter].l;
	    value = exp_struc_f[f_counter].value;

/* test for creating symmetry related reflections 
   fprintf (testfile, "%f %f %f %f 0.1 \n", h, k, l, value); */

	    
	    /******* Calculate structure factor from model   ********/
		    
	    /* Square of the scattering length  */
/*	    scattering_length_sq = (4.0/3.0/(lattice_par.a*lattice_par.a) * (h*h + k*k + h*k) + l*l / (lattice_par.c*lattice_par.c)) / 4.0; */
	    scattering_length_sq = calc_scattering_length_sq(lattice_par, h, k, l);

	    calc_struc_f_square_sum = 0.0;

	    /* Calculate the structure factor by summing the terms. */
 	    for (domain_counter = 0; domain_counter < nr_domains; domain_counter++) {	/* Loop over domains */
		calc_struc_f_real = 0.0;
		calc_struc_f_im = 0.0;
		for (atom_counter = 0; atom_counter < nr_atoms[domain_counter]; ++atom_counter) {
		    /* Add the contribution of atom j to the structure factor. */

/*		    dw_factor = exp(-1.0 * atoms[atom_counter].dw_par * scattering_length_sq);	*/ /* Debye-Waller factor  */
/*		    calc_struc_f_real += (atoms[atom_counter].occ * atom_scattering_f * dw_factor * cos(2.0*PI*(h*atoms[atom_counter].x+k*atoms[atom_counter].y+l*atoms[atom_counter].z))) ;
		    calc_struc_f_im +=  (atoms[atom_counter].occ * atom_scattering_f * dw_factor * sin(2.0*PI*(h*atoms[atom_counter].x+k*atoms[atom_counter].y+l*atoms[atom_counter].z))); */

		    atom_scattering_f = atom_scatt_f(atoms[domain_counter][atom_counter].el_nr, scattering_length_sq); /* Atomic scattering factor */
		    dw_factor = exp(-1.0 * atoms[domain_counter][atom_counter].dw_par * scattering_length_sq);	/* Debye-Waller factor  */
		    calc_struc_f_real += (atoms[domain_counter][atom_counter].occ * atom_scattering_f * dw_factor * cos(2.0*PI*(h*atoms[domain_counter][atom_counter].x+k*atoms[domain_counter][atom_counter].y+l*atoms[domain_counter][atom_counter].z))) ;
		    calc_struc_f_im +=  (atoms[domain_counter][atom_counter].occ * atom_scattering_f * dw_factor * sin(2.0*PI*(h*atoms[domain_counter][atom_counter].x+k*atoms[domain_counter][atom_counter].y+l*atoms[domain_counter][atom_counter].z)));
		};
		/* Calculate the absolute value of the structure factor*/
		/*calc_struc_f_abs[domain_counter] = sqrt(calc_struc_f_real*calc_struc_f_real + calc_struc_f_im*calc_struc_f_im);	    */
		calc_struc_f_abs[domain_counter] = sqrt(calc_struc_f_real*calc_struc_f_real + calc_struc_f_im*calc_struc_f_im);
		calc_struc_f_square_sum += domain_occ[domain_counter] * calc_struc_f_abs[domain_counter] * calc_struc_f_abs[domain_counter];

	    

		/* Calculate the phase factor */
/*	    calc_struc_f_phase = atan2(calc_struc_f_im, calc_struc_f_real); */
		calc_struc_f_phase[domain_counter] = atan2(calc_struc_f_im, calc_struc_f_real);

	    }; 
/* if (r1 ==0.2 && r2 ==0.1) {
 fprintf (stdout, "%f %f %f %f %f %f %f %f\n", h,k,l,value, calc_struc_f_abs, calc_struc_f_phase*360.0/2.0/PI,  360.0*(h*r1+k*r2),value-calc_struc_f_abs);
 }; */

	    /* Sum all reflections for value of electron difference map */
/* 	    electron_diff_value += (value - scale * calc_struc_f_abs) * cos((2.0*PI * (h * r1 + k * r2)) - calc_struc_f_phase); */

 	    for (domain_counter = 0; domain_counter < nr_domains; domain_counter++) {	/* Loop over domains */
		domain_frac = calc_struc_f_abs[domain_counter] / sqrt(calc_struc_f_square_sum); 	/* Approximate fraction of experimental structure factor corresponding to current domain */
		electron_diff_value += (value/scale * domain_frac -  calc_struc_f_abs[domain_counter]) * cos((2.0*PI * (h * r1 + k * r2)) - calc_struc_f_phase[domain_counter]);
	    };
	};
    };
		
    /* Normalize with unit cell area */
    unit_cell_area = lattice_par.a * lattice_par.b * sin (lattice_par.gamma); /* Area of the surface unit cell (Ang^2) */
    electron_diff_value = electron_diff_value / unit_cell_area;	
	
/* test for creating symmetry related reflections 
   fclose (testfile); */

/* if (r1 ==0.2 && r2 ==0.1) {
fprintf (stdout, "\n\n");	
}; */

    return electron_diff_value;
}	  


/************************************************************************************************
 * double find_scale(const struct structure_factor exp_struc_f[], const size_t nr_f,  struct atom_structure  atoms[MAX_DOMAINS][MAX_ATOMS], const struct atom_structure atoms_bulk[MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth, const size_t nr_atoms[], const double domain_occ[], const size_t nr_domains, size_t nr_atoms_bulk) 
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
 * domain_occ: Occupation of domain                                                             *
 * nr_domains: Number of domains.								*
 ************************************************************************************************/
double find_scale(const struct structure_factor exp_struc_f[], const size_t nr_f,  struct atom_structure  atoms[MAX_DOMAINS][MAX_ATOMS], const struct atom_structure atoms_bulk[MAX_ATOMS], const struct latt_par_structure lattice_par, double penetration_depth, const size_t nr_atoms[], const double domain_occ[], const size_t nr_domains, size_t nr_atoms_bulk) 
{
    size_t f_counter;		
    double struc_factor;		/* Absolute value of structure factor */ 
    double sum1 = 0.0;			/* Summing for scale */
    double sum2 = 0.0;			
    double scale;			/* Scale factor */

    for (f_counter = 0; f_counter < nr_f; f_counter++) {	/* Loop over all structure factors */

        struc_factor = calc_structure_factor(atoms, atoms_bulk, lattice_par, penetration_depth, exp_struc_f[f_counter].h, exp_struc_f[f_counter].k, exp_struc_f[f_counter].l, nr_atoms, domain_occ, nr_domains, nr_atoms_bulk);

	
	/* This is for scale determination */
	sum1 +=  exp_struc_f[f_counter].value * struc_factor / (exp_struc_f[f_counter].sigma*exp_struc_f[f_counter].sigma);
	sum2 +=  struc_factor*struc_factor / (exp_struc_f[f_counter].sigma*exp_struc_f[f_counter].sigma);

    };

    /* Average scale (weighted with error) */
    scale = sum1 / sum2;

    return scale;
}

/************************************************************************
 * double atom_scatt_f(int element_nr, scattering_length_sq)		*
 * 									*
 * Calculate the atomic scattering factor with the Cromer-Mann formula  *
 * 									*
 * element_nr: Element (atomic)number of atom				*
 * scattering_length_sq: Square of the scattering vector		*
 ************************************************************************/

double atom_scatt_f(int element_nr, double scattering_length_sq)
{
    unsigned int cm_counter;	/* counter used to sum the terms of the Cromer-Mann formula */
    double atom_scattering_f = 0.0;		/* Atomic scattering factor */

    /* Check if element data is available (doesn't check if nr > 14 !) */
    if (cm_par[element_nr][0].a == 0.0) {
	fprintf (stderr, "ERROR: Atomic scattering data is not available for element with atomic number %d!\n", element_nr);
	exit(8);
    };    

     /* Calculate the atomic scattering factor with the Cromer-Mann formula */
    for (cm_counter = 0; cm_counter < 4; ++cm_counter) {			/* Sum the first 4 terms of the Cromer-Mann formula  */
	atom_scattering_f = atom_scattering_f + cm_par[element_nr][cm_counter].a * exp(-cm_par[element_nr][cm_counter].b * scattering_length_sq);
    };
    atom_scattering_f = atom_scattering_f + cm_par[element_nr][0].c;          /* Constant term of the Cromer-Mann formula */
    return atom_scattering_f;

}


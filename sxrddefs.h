/* Defintions for sxrdcalc */

/* Already included? */
#ifndef SXRDDEFS_INCLUDED
#define SXRDDEFS_INCLUDED

/* For debugging */
/* #define DEBUG1 */
/* #include <mcheck.h> */

/* Include files   */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>


/* Define version and date */
#define VERSION "3.3.3"
#define VERSION_DATE "August 2019"

/* Define some maximal numbers */
#define MAX_STRUCTURE_FACTORS 100000	/* Maximum number of structure factors */
#define MAX_ATOMS 1000			/* Maximum number of atoms in the unit cell */
#define FILENAME_LENGTH 255		/* Lenght of file names */
#define LINE_LENGTH 255			/* Maximal length of a line in the input file */
#define MAX_FREE_PARAMETERS (size_t) 100	/* The number of free parameters in the model */
#define MAX_INPUT_FILE_LINES 2000	/* Maximal length of the input file */
#define MAX_DOMAINS 10			/* Maximal number of different domains */

/* Put following into input file */
#define STEP_PATTERSON_R1 100		/* Number of steps in the Patterson map in r1 direction */
#define STEP_PATTERSON_R2 100		/* Number of steps in the Patterson map in r2 direction */

/* not used */
#define MAX_PATTERSON_R1 1		/* Maximum r1 in the Patterson map (in fractional coordinates) */
#define MAX_PATTERSON_R2 1		/* Maximum r2 in the Patterson map (in fractional coordinates) */

/* false and true */
#define FALSE 0
#define TRUE 1



/************************************************
 *  Global Structures				*
 ************************************************/

struct latt_par_structure {			/* Lattice parameters a,b,c and angles alpha, beat, gamma */
    double a;
    double b;
    double c;
    double alpha;
    double beta;
    double gamma;  
} ;


struct structure_factor{ 	/* Structure containing the structure factor */
    double h;			/* h index of structure factor */
    double k;			/* k index of structure factor */
    double l;			/* l index of structure factor */
    double value;		/* value of structure factor */
    double sigma;		/* error of structure factor */
};


struct atom_structure {			/* Structure containing the data of atoms in unit cell*/
    int el_nr; 			/* Number of element */
    double x;			/* atom position in unit cell x (in fractional coordinates)*/
    double y;			/* atom position in unit cell y (in fractional coordinates)*/
    double z;			/* atom position in unit cell z (in fractional coordinates)*/
    double scatt_f;		/* atomic scattering factor */
    double dw_par;		/* Debye-Waller factor of the atoms  */
    double occ;			/* Occupany of atom */
    double aniso_dw[6];          /* anisotropic Debye-Waller factor (b11,b22,b33,b12,b13,23) */
} ;

struct displacements {			/* Structure for atomic displacements used in fit */
    int nr1;			/* Number of parameter for 1st displacement vector */
    double dx1;			/* 1st displacement vector x (fractional coordinates) */
    double dy1;			/* 1st displacement vector y (fractional coordinates) */
    double dz1;			/* 1st displacement vector z (fractional coordinates) */
    int nr2;			/* Number of parameter for 2st displacement vector */
    double dx2;			/* 2st displacement vector x (fractional coordinates) */
    double dy2;			/* 2st displacement vector y (fractional coordinates) */
    double dz2;			/* 2st displacement vector z (fractional coordinates) */
    int nr3;			/* Number of parameter for 3st displacement vector */
    double dx3;			/* 3st displacement vector x (fractional coordinates) */
    double dy3;			/* 3st displacement vector y (fractional coordinates) */
    double dz3;			/* 3st displacement vector z (fractional coordinates) */
    int nr_dw;			/* Number of parameter for Debye-Waller parameter */
    double dw_scale;		/* Scale for Debye-Waller parameter */
    int nr_occ;			/* Number of parameter for occupancy */
    int nr_aniso_dw[6];          /* Numbers of parameters for anisotropic Debye-Waller factor */
    double aniso_dw_scale[6][6];          /* anisotropic Debye-Waller factor (6 vectors with the six elements (b11,b22,b33,b12,b13,23)) */
};

struct cromer_mann_par {					/* Structure for Cromer-Mann parameters used to calculate atomic scattering factors */
    double a;							/* Parameter a */
    double b;							/* Parameter b */
    double c;							/* Parameter c */
};


/********************************************************************
 * Global variable declarations  (real definition in sxrdcalc.c)    *
 ********************************************************************/

#ifdef UNDEF
extern char input_struc_file[FILENAME_LENGTH] ;   /* Input structure file name */
extern char calc_struc_f_file[FILENAME_LENGTH];       /* Filename to which structure factors should be written */
extern char exp_struc_f_file[FILENAME_LENGTH];		/* Experimental structure factors file */
extern char patt_file_name[FILENAME_LENGTH]; 	/* Output file for Patterson map */
#endif

/* Some structures */ 
/* extern struct latt_par_structure lattice_par; */ 	/* Lattice parameters a,b,c and angles alpha, beat, gamma */
extern struct structure_factor struc_f_calc[MAX_STRUCTURE_FACTORS];		/* calculated structure factor */
extern struct structure_factor struc_f_exp[MAX_STRUCTURE_FACTORS];		/* experimental structure factor */
/*extern struct atom_structure  atoms[MAX_ATOMS]; */

/* symmetry matrix */
extern double symmetry_matrix [12][2][2];
extern unsigned int nr_symmetry_matrices;


/********************************************************
 * Constants 						*
 ********************************************************/

extern const double PI;		/*  Pi from Windows calculator */


/********************************************************
 * Enumerations enum					*
 ********************************************************/
						  
/* Define the symmetries for input to set_symmetry_matrix */
enum SURFACE_SYMMETRY {p1, p6, p6mm};	

enum output {no_intermediate_status, intermediate_status, short_output};

#endif /* SXRDDEFS_INCLUDED */

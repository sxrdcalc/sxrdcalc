/****************************************************************************************
 *											*
 *    Program name: sxrdcalc 								* 
 *    Version: 3.32 - January 2011							*
 *    Written by W. Voegeli 								*
 * 											*
 *    Does calculations for surface x-ray diffraction.					* 
 *											*
 *    Input variables are read from first command line argument.			*
 * 											*
 *    Can calculate at present:								*
 *     - Patterson map from experimental data						*
 *     - Electron difference map (only for one domain!)					*
 *     - Structure factors 								*
 *     - Levenberg-Marquardt least-squares fit of model to experimental intensities	*
 *     - Simulated annealing fit							*
 *											*
 *    NOTE: Patterson map and electron difference map only for hexagonal system with gamma=120deg and p6 or p6mm symmetry.	*
 *	    For other geometries, add symmtry matrices.    *
 ****************************************************************************************/

/****************************************************************************************
 * Uses functions in following files:							*
 * - sxrdcalc.c	(main program, input file prasing)					*
 * - sxrdio.c   (some input and output functions)					*
 * - calcs.c    (functions that do the actual calculations)				*
 * - standard c library functions							*
 * - GNU Scientific Library								*
 ****************************************************************************************/
 

/*******************************************
 * Include files			   *
 *******************************************/

#include "sxrddefs.h"     /* Global definitions */
#include "sxrdio.h"	  /* I/O */
#include "calcs.h"	  /* Some calculations */
#include <gsl/gsl_siman.h>

/*******************************************
 * Global Variable declarations            *
 *******************************************/

/* symmetry matrix */
enum SURFACE_SYMMETRY symmetry;
double symmetry_matrix [12][2][2];
unsigned int nr_symmetry_matrices;



/* Global structures */
struct structure_factor struc_f_calc[MAX_STRUCTURE_FACTORS];		/* calculated structure factor */
struct structure_factor struc_f_exp[MAX_STRUCTURE_FACTORS];		/* experimental structure factor */


const double PI=3.1415926535897932384626433832795;		/*  Pi from Windows calculator */

/********************************
 * Variables local to this file *
 ********************************/
/* Filenames */
static char *program_name;                      /* Name of the program (from the command line) */
static char *input_file_name;			/* Name of the input file (from the command line) */
static char input_struc_file[MAX_DOMAINS][FILENAME_LENGTH];	/* Input structure file names */
static char fit_struc_file[MAX_DOMAINS][FILENAME_LENGTH];	/* Input file names for structure refinement */
static char bulk_struc_file[FILENAME_LENGTH] = "";               /* Input file name for bulk structure */
static char calc_struc_f_file[FILENAME_LENGTH];   		/* Filename for output structure factors  */
static char exp_struc_f_file[FILENAME_LENGTH];			/* Experimental structure factors file */
static char patt_file_name[FILENAME_LENGTH];			/* Output file for Patterson map */
static char fit_coord_file[MAX_DOMAINS][FILENAME_LENGTH];	/* Output file names for fitted coordinates */
static char el_diff_file_name[FILENAME_LENGTH];			/* Output file for electron difference map */

static unsigned int nr_fit_coord_files_read;     	/* Number of output files for fitted coordinates read */

static enum {patterson, theo_range, theo_data, ls_fit, sa_fit, electron_difference_map, electron_map} todo;

static enum output print_intermediate = no_intermediate_status; 	/* Print intermediate information in simulated annealing? */

static unsigned long int n_fit = 1;				/* Number of fits to start with different starting parameters */
static double par_var_width = 0.0;			/*Width of the gaussian function which is used to generate the starting parameters for the different fits */
static char comment[LINE_LENGTH];		/* Comment for output file */ 

static size_t nr_atoms[MAX_DOMAINS];   		/* Number of atoms in input unit cell */
static size_t nr_atoms_bulk = 0;   		/* Number of atoms in bulk unit cell */

static size_t nr_struc_f_exp;		/* Total number of experimental structure factors */
static size_t nr_struc_f_calc;		/* Total number of calculated structure factors */

static size_t nr_domains = 1;		/* number of domains */
static double domain_occ[MAX_DOMAINS];	/* Occupancy of each domain */

static double penetration_depth = 1.0E50; /* penetration depth of X-rays into crystal [reciprocal units of lattice vector] */

static double scale;			/* Scale factor between experimental and calculated structure factors */
static double theo_scale = 1.0;		/* Scale factor for theoretically calculated structure factors */

static unsigned  long int rng_seed = 0; 			/* Seed for random number generator for simulated annealing */

/* Structures */
static struct atom_structure  atoms_unit_cell[MAX_DOMAINS][MAX_ATOMS];     	/* Structure containing the data of atoms in unit cell*/
static struct atom_structure  atoms_bulk[MAX_ATOMS];     	/* Structure containing the data of atoms in bulk unit cell*/
static struct atom_structure  end_positions[MAX_DOMAINS][MAX_ATOMS];     	/* Structure containing the fitted coordinates of atoms */
static struct latt_par_structure lattice_par; 	/* Lattice parameters a,b,c and angles alpha, beta, gamma */


static struct displacements atom_displacements[MAX_DOMAINS][MAX_ATOMS];	/* Atomi displacement vectors used in fit */

static gsl_siman_params_t sa_par;					/* Parameters for simulated annealing */

static double start_par[MAX_FREE_PARAMETERS];					/* Starting values for the parameters */

static int fixed_parameters[MAX_FREE_PARAMETERS];	/* Fixed parameters (TRUE=fixed) */

/* Boolean */
static int do_read = 0;			/* Was the do identifier read? */

/* Step sizes for Patterson map */
static unsigned long int step_patterson_r1 = STEP_PATTERSON_R1;		 
static unsigned long int step_patterson_r2 = STEP_PATTERSON_R2; 

/* Step sizes for electron difference map */
static unsigned long int step_el_diff_r1 = STEP_PATTERSON_R1;		 
static unsigned long int step_el_diff_r2 = STEP_PATTERSON_R1; 

/*  Ranges and step sizes of h, k, l */
static double min_h = 0.0;
static double max_h = 10.0;
static double step_h = 1.0;
static double min_k = 0.0;
static double max_k = 10.0;
static double step_k = 1.0;
static double min_l = 0.25;
static double max_l = 0.25;
static double step_l = 1.0;

/* Lowest index bulk reflections */
static double bulk_h = 12345678.12345678;
static double bulk_k = 12345678.12345678;
static double bulk_l = 12345678.12345678;

/* For fit */
static unsigned long int max_iteration = 100;		/*  Maximum number of iterations */ 
static double  delta_abs = 1.0;		/* Absolute convergence criterion */
static double  delta_rel = 1.0;		/* Relative convergence criterion (not much sense)*/

/**************************************
 * Functions                          *
 **************************************/



/********************************************************
 * parse_line(char identifier[], char value[])		*
 *							*
 * Parse line of input file	      			*
 ********************************************************/

void parse_line(const char identifier[LINE_LENGTH], const char value[LINE_LENGTH])
{	
    char **endptr = NULL;  	/* A pointer for strtod, strtol */
    int identifier_parsed = 0;	/* Boolean: Was the identifier correctly parsed? */
 
    /* 'do' : what to calculate */
    if (strcmp(identifier, "do") == 0)      { 
	identifier_parsed = 1;		

	if (do_read) {					/* do already read ? */
	    fprintf (stderr, "Several \"do\" identifiers in input file!\n");
	    exit(8);
	};
	do_read = 1;
    
	if (strcmp(value, "patterson") == 0)	  {     /* Calculate Patterson map */ 
	    todo = patterson;
      
	};
	
	if (strcmp(value, "electron_diff_map") == 0)	  	  {	/* Calculate electron difference map */
	    todo = electron_difference_map;
	};

	if (strcmp(value, "electron_map") == 0)	  	  {	/* Calculate electron map */
	    todo = electron_map;
	};
    
	if (strcmp(value, "theo_range") == 0)	  {	/* Calculate structure factors from model structure */
	    todo = theo_range;				/* for range of h, k, l */
	};

	if (strcmp(value, "theo_data") == 0)	  {	/* Calculate structure factors from model structure */
	    todo = theo_data;					/* (same as in structure factor input file) */
	};

	if (strcmp(value, "ls_fit") == 0)	  	  {	/* Structure refinement by least-squares fit */
	    todo = ls_fit;
	};

	if (strcmp(value, "sa_fit") == 0)	  	  {	/* Structure refinement by simulated annealing fit */
	    todo = sa_fit;
	};

    
    };   /* end 'todo' */
  
    /* Filenames */
    if (strcmp(identifier, "patterson_file") == 0)      {   /* Filename for Patterson map */
	strcpy (patt_file_name, value);
	identifier_parsed = 1;		
    };
  
    if (strcmp(identifier, "electron_diff_file") == 0)      {   /* Filename for electron difference map */
	strcpy (el_diff_file_name, value);
	identifier_parsed = 1;		
    };

    if (strcmp(identifier, "f_in_file") == 0)      {   /* Filename for input structure factors */
	strcpy (exp_struc_f_file, value);
	identifier_parsed = 1;		
    };

    if (strcmp(identifier, "f_out_file") == 0)      {   /* Filename for output structure factors */
	strcpy (calc_struc_f_file, value);
	identifier_parsed = 1;		
    };

    if (strcmp(identifier, "bulk_struc_in_file") == 0)      {   /* Filename for bulk structure */
	strcpy (&(bulk_struc_file[0]), value);
	identifier_parsed = 1;		
    };

    if (strcmp(identifier, "struc_in_file1") == 0)      {   /* Filename for model structure file (domain 1)*/
	strcpy (&(input_struc_file[0][0]), value);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "struc_in_file2") == 0)      {   /* Filename for model structure file (domain 2)*/
	strcpy (&(input_struc_file[1][0]), value);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "struc_in_file3") == 0)      {   /* Filename for model structure file (domain 3)*/
	strcpy (&(input_struc_file[2][0]), value);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "struc_in_file4") == 0)      {   /* Filename for model structure file (domain 4)*/
	strcpy (&(input_struc_file[3][0]), value);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "struc_in_file5") == 0)      {   /* Filename for model structure file (domain 5)*/
	strcpy (&(input_struc_file[4][0]), value);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "struc_in_file6") == 0)      {   /* Filename for model structure file (domain 6)*/
	strcpy (&(input_struc_file[5][0]), value);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "struc_in_file7") == 0)      {   /* Filename for model structure file (domain 7)*/
	strcpy (&(input_struc_file[6][0]), value);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "struc_in_file8") == 0)      {   /* Filename for model structure file (domain 8)*/
	strcpy (&(input_struc_file[7][0]), value);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "struc_in_file9") == 0)      {   /* Filename for model structure file (domain 9)*/
	strcpy (&(input_struc_file[8][0]), value);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "struc_in_file10") == 0)      {   /* Filename for model structure file (domain 10)*/
	strcpy (&(input_struc_file[9][0]), value);
	identifier_parsed = 1;		
    };

    if (strcmp(identifier, "fit_struc_in_file1") == 0)      {   /* Filename for input file for structure refinement (domain 1)*/
	strcpy (&(fit_struc_file[0][0]), value);
	identifier_parsed = 1;		
    };
  
    if (strcmp(identifier, "fit_struc_in_file2") == 0)      {   /* Filename for input file for structure refinement (domain 2)*/
	strcpy (&(fit_struc_file[1][0]), value);
	identifier_parsed = 1;		
    };
  
    if (strcmp(identifier, "fit_struc_in_file3") == 0)      {   /* Filename for input file for structure refinement (domain 3)*/
	strcpy (&(fit_struc_file[2][0]), value);
	identifier_parsed = 1;		
    };
  
    if (strcmp(identifier, "fit_struc_in_file4") == 0)      {   /* Filename for input file for structure refinement (domain 4)*/
	strcpy (&(fit_struc_file[3][0]), value);
	identifier_parsed = 1;		
    };
  
    if (strcmp(identifier, "fit_struc_in_file5") == 0)      {   /* Filename for input file for structure refinement (domain 5)*/
	strcpy (&(fit_struc_file[4][0]), value);
	identifier_parsed = 1;		
    };
  
    if (strcmp(identifier, "fit_struc_in_file6") == 0)      {   /* Filename for input file for structure refinement (domain 6)*/
	strcpy (&(fit_struc_file[5][0]), value);
	identifier_parsed = 1;		
    };
  
    if (strcmp(identifier, "fit_struc_in_file7") == 0)      {   /* Filename for input file for structure refinement (domain 7)*/
	strcpy (&(fit_struc_file[6][0]), value);
	identifier_parsed = 1;		
    };
  
    if (strcmp(identifier, "fit_struc_in_file8") == 0)      {   /* Filename for input file for structure refinement (domain 8)*/
	strcpy (&(fit_struc_file[7][0]), value);
	identifier_parsed = 1;		
    };
  
    if (strcmp(identifier, "fit_struc_in_file9") == 0)      {   /* Filename for input file for structure refinement (domain 9)*/
	strcpy (&(fit_struc_file[8][0]), value);
	identifier_parsed = 1;		
    };
  
    if (strcmp(identifier, "fit_struc_in_file10") == 0)      {   /* Filename for input file for structure refinement (domain 10)*/
	strcpy (&(fit_struc_file[9][0]), value);
	identifier_parsed = 1;		
    };
  
    if (strcmp(identifier, "fit_coord_out_file1") == 0)      {   /* Filename for output of fitted coordinates (domain 1)*/
	strcpy (&(fit_coord_file[0][0]), value);
	identifier_parsed = 1;		
	++nr_fit_coord_files_read; 
    };
  
    if (strcmp(identifier, "fit_coord_out_file2") == 0)      {   /* Filename for output of fitted coordinates (domain 2)*/
	strcpy (&(fit_coord_file[1][0]), value);
	identifier_parsed = 1;		
	++nr_fit_coord_files_read; 
    };
  
    if (strcmp(identifier, "fit_coord_out_file3") == 0)      {   /* Filename for output of fitted coordinates (domain 3)*/
	strcpy (&(fit_coord_file[2][0]), value);
	identifier_parsed = 1;		
	++nr_fit_coord_files_read; 
    };
  
    if (strcmp(identifier, "fit_coord_out_file4") == 0)      {   /* Filename for output of fitted coordinates (domain 4)*/
	strcpy (&(fit_coord_file[3][0]), value);
	identifier_parsed = 1;		
	++nr_fit_coord_files_read; 
    };
  
    if (strcmp(identifier, "fit_coord_out_file5") == 0)      {   /* Filename for output of fitted coordinates (domain 5)*/
	strcpy (&(fit_coord_file[4][0]), value);
	identifier_parsed = 1;		
	++nr_fit_coord_files_read; 
    };
  
    if (strcmp(identifier, "fit_coord_out_file6") == 0)      {   /* Filename for output of fitted coordinates (domain 6)*/
	strcpy (&(fit_coord_file[5][0]), value);
	identifier_parsed = 1;		
	++nr_fit_coord_files_read; 
    };
  
    if (strcmp(identifier, "fit_coord_out_file7") == 0)      {   /* Filename for output of fitted coordinates (domain 7)*/
	strcpy (&(fit_coord_file[6][0]), value);
	identifier_parsed = 1;		
	++nr_fit_coord_files_read; 
    };
  
    if (strcmp(identifier, "fit_coord_out_file8") == 0)      {   /* Filename for output of fitted coordinates (domain 8)*/
	strcpy (&(fit_coord_file[7][0]), value);
	identifier_parsed = 1;		
	++nr_fit_coord_files_read; 
    };
  
    if (strcmp(identifier, "fit_coord_out_file9") == 0)      {   /* Filename for output of fitted coordinates (domain 9)*/
	strcpy (&(fit_coord_file[8][0]), value);
	identifier_parsed = 1;		
	++nr_fit_coord_files_read; 
    };
  
    if (strcmp(identifier, "fit_coord_out_file10") == 0)      {   /* Filename for output of fitted coordinates (domain 10)*/
	strcpy (&(fit_coord_file[9][0]), value);
	identifier_parsed = 1;		
	++nr_fit_coord_files_read; 
    };
  
    /* Number of domains */
    if (strcmp(identifier, "nr_domains") == 0) {
	nr_domains = strtoul(value, endptr, 10);
	identifier_parsed = 1;		
    };

    /* Domain occupancies */
    if (strcmp(identifier, "domain_occ1") == 0) {		/* Domain occupancy 1 */
	domain_occ[0] = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "domain_occ2") == 0) {		/* Domain occupancy 2 */
	domain_occ[1] = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "domain_occ3") == 0) {		/* Domain occupancy 3 */
	domain_occ[2] = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "domain_occ4") == 0) {		/* Domain occupancy 4 */
	domain_occ[3] = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "domain_occ5") == 0) {		/* Domain occupancy 5 */
	domain_occ[4] = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "domain_occ6") == 0) {		/* Domain occupancy 6 */
	domain_occ[5] = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "domain_occ7") == 0) {		/* Domain occupancy 7 */
	domain_occ[6] = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "domain_occ8") == 0) {		/* Domain occupancy 8 */
	domain_occ[7] = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "domain_occ9") == 0) {		/* Domain occupancy 9 */
	domain_occ[8] = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "domain_occ10") == 0) {     	/* Domain occupancy 10 */
	domain_occ[9] = strtod(value, endptr);
	identifier_parsed = 1;		
    };

    /* Penetration depth of X-rays into crystal (reciprocal units of lattice paramters) */
    if (strcmp(identifier, "penetration_depth")  == 0) {
	penetration_depth = strtod(value,  endptr);
	identifier_parsed = 1;		
    };

    /* Symmetry */
    if (strcmp(identifier, "symmetry") == 0)      {   /* Symmetry */ 
	identifier_parsed = 1;		
	if (strcmp(value, "p6") == 0) {	/* p6 symmetry */
	    symmetry = p6;
	};

	if (strcmp(value, "p6mm") == 0) {	/* p6mm symmetry */
	    symmetry = p6mm;
	};
    };

    /* Number of steps for Patterson map */
    if (strcmp(identifier, "step_patterson_r1") == 0) {
	step_patterson_r1 = strtoul(value, endptr, 10);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "step_patterson_r2") == 0) {
	step_patterson_r2 = strtoul(value, endptr, 10);
	identifier_parsed = 1;		
    };
    
    /* Number of steps for electron difference map */
    if (strcmp(identifier, "step_el_diff_r1") == 0) {
	step_el_diff_r1 = strtoul(value, endptr, 10);
	identifier_parsed = 1;		
    };
 
    if (strcmp(identifier, "step_el_diff_r2") == 0) {
	step_el_diff_r2 = strtoul(value, endptr, 10);
	identifier_parsed = 1;		
    };

    /* Scale factor for theoretically calculated structure factors */
    if (strcmp(identifier, "theo_scale")  == 0) {
	theo_scale = strtod(value,  endptr);
	identifier_parsed = 1;		
    };
 
    /* Minimum values , maximum values and step sizes  for h, k, l */
    if (strcmp(identifier, "minh") == 0) {
	min_h = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "maxh") == 0) {
	max_h = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "steph") == 0) {
	step_h = strtod(value, endptr);
	identifier_parsed = 1;		
    };

    if (strcmp(identifier, "mink") == 0) {
	min_k = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "maxk") == 0) {
	max_k = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "stepk") == 0) {
	step_k = strtod(value, endptr);
	identifier_parsed = 1;		
    };

    if (strcmp(identifier, "minl") == 0) {
	min_l = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "maxl") == 0) {
	max_l = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "stepl") == 0) {
	step_l = strtod(value, endptr);
	identifier_parsed = 1;		
    };	

    /* Get index of bulk reflections */
    if (strcmp(identifier, "bulkh") == 0) {
	bulk_h = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "bulkk") == 0) {
	bulk_k = strtod(value, endptr);
	identifier_parsed = 1;		
    };
    if (strcmp(identifier, "bulkl") == 0) {
	bulk_l = strtod(value, endptr);
	identifier_parsed = 1;		
    };

    /********** Fit general ***************/

    if (strcmp(identifier, "n_fit") == 0) {			/* Number of fits to start with different starting parameters */
	n_fit =  strtoul(value, endptr,10);
	identifier_parsed = 1;		
    };

    if (strcmp(identifier, "par_var_width") == 0) {			/* Width of the gaussian function which is used to generate the starting parameters for the different simulations */
	par_var_width =  strtod(value, endptr);
	identifier_parsed = 1;		
    };

    /********* For least-squares fit ***************/



    if  (strcmp(identifier, "max_iteration") == 0) {	/* Maximum number of iterations */
	max_iteration = strtoul(value, endptr,10);
	identifier_parsed = 1;		
    };

    if  (strcmp(identifier, "delta_abs") == 0) {		/* Absolute convergence criterion */				
	delta_abs = strtod(value, endptr);
	identifier_parsed = 1;		
    };

    if  (strcmp(identifier, "delta_rel") == 0) {		/* Relative convergence criterion (not much sense) */		
	delta_rel = strtod(value, endptr);
	identifier_parsed = 1;		
    };

    /********* For simulated annealing ***********/
    if  (strcmp(identifier, "print_intermediate") == 0) {/* Print intermediate information in simulated annealing? */
	if (strcmp(value, "yes") == 0) {
	    if (print_intermediate != short_output) {
	        print_intermediate = intermediate_status;
	    };
	    identifier_parsed = 1;		
	};
	if (strcmp(value, "no") == 0) {		
	    if (print_intermediate != short_output) {
	        print_intermediate = no_intermediate_status;
	    };
	    identifier_parsed = 1;		
	};
    };

    if  (strcmp(identifier, "short_output") == 0) {/* Amount of output? */
	if (strcmp(value, "yes") == 0) {
	    print_intermediate = short_output;
	    identifier_parsed = 1;		
	};
	if (strcmp(value, "no") == 0) {		
	    identifier_parsed = 1;		
	};
    };

    if  (strcmp(identifier, "n_tries") == 0) {			/* 	The number of points to try for each step */
	sa_par.n_tries =  (int) strtol(value, endptr,10);
	identifier_parsed = 1;		
    };

    if  (strcmp(identifier, "iters_fixed_T") == 0) {		/* The number of iterations at each temperature */
	sa_par.iters_fixed_T =  (int) strtol(value, endptr,10);
	identifier_parsed = 1;		
    };

    if  (strcmp(identifier, "step_size") == 0) {		/* The maximum step size in the random walk */
	sa_par.step_size =  strtod(value, endptr);
	identifier_parsed = 1;		
    };

    if  (strcmp(identifier, "sa_k") == 0) {			/* Boltzmann K */ 
	sa_par.k             =  strtod(value, endptr);
	identifier_parsed = 1;		
    };

    if  (strcmp(identifier, "t_initial") == 0) {		/* Initial temperatur */
	sa_par.t_initial =  strtod(value, endptr);
	identifier_parsed = 1;		
    };

    if  (strcmp(identifier, "mu_t") == 0) {		/* damping factor for temperature */
	sa_par.mu_t =  strtod(value, endptr);
	identifier_parsed = 1;		
    };

    if  (strcmp(identifier, "t_min") == 0) {		/* Minimum temperature */
	sa_par.t_min =  strtod(value, endptr);
	identifier_parsed = 1;		
    };

    if  (strcmp(identifier, "sa_seed") == 0) {		/* Seed for random number generator */
	rng_seed =  strtoul(value, endptr, 10);
	identifier_parsed = 1;		
    };

    if  (strcmp(identifier, "rng_seed") == 0) {		/* Seed for random number generator */
	rng_seed =  strtoul(value, endptr, 10);		/* SAME AS ABOVE ! */
	identifier_parsed = 1;		
    };


    /* No match found? -> Print error message. */
    if (!identifier_parsed) {	
	fprintf (stderr, "Error: %s %s not understood in input file!", identifier, value);
	exit(8);
    };	    

    return;
}

/************************************************************************
 * void print_input(void)						*
 *									*
 * Prints the input read from the input file.				*
 ************************************************************************/

void print_input(void)
{      
    /* Program and version */
    fprintf(stdout, "\n---------------------------------------\n");
    fprintf(stdout, "Program %s for surface x-ray diffraction calculations.\n", program_name);
    fprintf(stdout, "Version %s - %s\n\n", VERSION, VERSION_DATE);

    /* Inputfile */
    fprintf(stdout, "\n Inputfile: %s\n", input_file_name);

    /* Patterson map? */
    if (todo == patterson) {
	fprintf (stdout, "Calculating Patterson map.\n");
	fprintf (stdout, "   %s\n" , comment);	
	fprintf (stdout, "Input file for structure factors: %s\n", exp_struc_f_file);
	fprintf (stdout, "Output file for Patterson map: %s\n\n", patt_file_name);
	switch (symmetry)
	    {
	    case p6:
		fprintf (stdout, "Surface symmetry used: p6.\n");    
		break;
	    case p6mm:
		fprintf (stdout, "Surface symmetry used: p6mm.\n");    
		break;
	    case p1:
		fprintf (stdout, "Surface symmetry used: p1.\n");    
		break;
	    default:
		fprintf (stderr, "Don't know which symmetry!\n");    
	    };
    };

    
    /* Electron difference map? */
    if (todo == electron_difference_map) {
	fprintf (stdout, "Calculating electron difference map.\n");
	fprintf (stdout, "   %s\n" , comment);	
	fprintf (stdout, "Input file for structure factors: %s\n", exp_struc_f_file);
	{
	    unsigned int domain_counter;
	    for (domain_counter = 0; domain_counter < nr_domains; ++domain_counter) {
		fprintf (stdout, "Model input file for domain %d: %s\n", domain_counter + 1, &(input_struc_file[domain_counter][0]));
		fprintf (stdout, "Occupancy of domain %d: %f\n", domain_counter + 1, domain_occ[domain_counter]);
	    }; 
	};
	fprintf (stdout, "Output file for electron difference map: %s\n\n", el_diff_file_name);
	switch (symmetry)
	    {
	    case p6:
		fprintf (stdout, "Surface symmetry used: p6.\n");    
		break;
	    case p6mm:
		fprintf (stdout, "Surface symmetry used: p6mm.\n");    
		break;
	    case p1:
		fprintf (stdout, "Surface symmetry used: p1.\n");    
		break;
	    default:
		fprintf (stderr, "Don't know which symmetry!\n");    
	    };
    };    

  /* Theoretical electron  map? */
    if (todo == electron_map) {
	fprintf (stdout, "Calculating electron map from model structure.\n");
	fprintf (stdout, "   %s\n" , comment);	
	{
	    unsigned int domain_counter;
	    for (domain_counter = 0; domain_counter < nr_domains; ++domain_counter) {
		fprintf (stdout, "Model input file for domain %d: %s\n", domain_counter + 1, &(input_struc_file[domain_counter][0]));
		fprintf (stdout, "Occupancy of domain %d: %f\n", domain_counter + 1, domain_occ[domain_counter]);
	    }; 
	};
	fprintf (stdout,"h from %.2f to %.2f, step %.2f\n",  min_h, max_h, step_h);
	fprintf (stdout, "k from %.2f to %.2f, step %.2f\n",  min_k, max_k, step_k);
	fprintf (stdout, "l from %.2f to %.2f, step %.2f\n",  min_l, max_l, step_l);
	
	switch (symmetry)
	    {
	    case p6:
		fprintf (stdout, "Surface symmetry used: p6.\n");    
		break;
	    case p6mm:
		fprintf (stdout, "Surface symmetry used: p6mm.\n");    
		break;
	    case p1:
		fprintf (stdout, "Surface symmetry used: p1.\n");    
		break;
	    default:
		fprintf (stderr, "Don't know which symmetry!\n");    
	    };
	fprintf (stdout, "Output file for electron map: %s\n\n", el_diff_file_name);
    };    

    /* Structure factors for range of h, k, l? */
    if (todo == theo_range) {
	fprintf (stdout, "Calculating theoretical structure factors from model.\n");    
	fprintf (stdout, "   %s\n" , comment);	
	fprintf (stdout, "h from %.2f to %.2f, step %.2f\n",  min_h, max_h, step_h);
	fprintf (stdout, "k from %.2f to %.2f, step %.2f\n",  min_k, max_k, step_k);
	fprintf (stdout, "l from %.2f to %.2f, step %.2f\n",  min_l, max_l, step_l);
	fprintf (stdout, "Scale factor: %.8f\n", theo_scale);
	if (strlen(bulk_struc_file) > 0) {
  	    fprintf (stdout, "Penetration depth: %.2f\n\n", penetration_depth);
	};
	fprintf (stdout, "Output file for structure factors: %s\n\n", calc_struc_f_file);
	if (strlen(bulk_struc_file) > 0) {
	    fprintf (stdout, "Input file for bulk structure: %s\n", bulk_struc_file);
	};
	fprintf (stdout, "Number of domains: %d \n", (int) nr_domains);
	{
	    unsigned int domain_counter;
	    for (domain_counter = 0; domain_counter < nr_domains; ++domain_counter) {
		fprintf (stdout, "Model input file for domain %d: %s\n", domain_counter + 1, &(input_struc_file[domain_counter][0]));
		fprintf (stdout, "Occupancy of domain %d: %f\n", domain_counter + 1, domain_occ[domain_counter]);
	    };
	};
    };

    /* Structure factors for same range as in experimental data? */
    if (todo == theo_data) {
	fprintf (stdout, "Calculating theoretical structure factors from model.\n");    
	fprintf (stdout, "   %s\n" , comment);	
    
	fprintf (stdout, "Input file for experimental structure factors: %s\n", exp_struc_f_file);
	fprintf (stdout, "Output file for structure factors: %s\n\n", calc_struc_f_file);
	fprintf (stdout, "Scale factor: %.8f\n", theo_scale);
	if (strlen(bulk_struc_file) > 0) {
	    fprintf (stdout, "Penetration depth: %.2f\n\n", penetration_depth);
	};
	if (strlen(bulk_struc_file) > 0) {
	    fprintf (stdout, "Input file for bulk structure: %s\n", bulk_struc_file);
	};
	fprintf (stdout, "Number of domains: %d \n", (int) nr_domains);
	{
	    unsigned int domain_counter;
	    for (domain_counter = 0; domain_counter < nr_domains; ++domain_counter) {
		fprintf (stdout, "Model input file for domain %d: %s\n", domain_counter + 1, &(input_struc_file[domain_counter][0]));
		fprintf (stdout, "Occupancy of domain %d: %f\n", domain_counter + 1, domain_occ[domain_counter]);
	    };
	};
    };
 
    /* Least-squares fit ? */
    if (todo == ls_fit) {
	fprintf (stdout, "Least-squares fit of model to experimental structure factors.\n");
	fprintf (stdout, "   %s\n" , comment);

	fprintf (stdout, "Number of domains: %d \n\n", (int) nr_domains);
	fprintf (stdout, "Input file for experimental structure factors: %s\n", exp_struc_f_file);
	if (strlen(bulk_struc_file) > 0) {
	    fprintf (stdout, "Input file for bulk structure: %s\n", bulk_struc_file);
	};
	if (strlen(bulk_struc_file) > 0) {
	    fprintf (stdout, "Penetration depth: %.2f\n\n", penetration_depth);
	};
	{
	    unsigned int domain_counter;
	    for (domain_counter = 0; domain_counter < nr_domains; ++domain_counter) {
		fprintf (stdout, "Fit Input file for domain %d: %s\n", domain_counter + 1, &(fit_struc_file[domain_counter][0]));
		fprintf (stdout, "Output file for domain %d: %s\n", domain_counter + 1, &(fit_coord_file[domain_counter][0]));
		fprintf (stdout, "Occupancy of domain %d: %f\n", domain_counter + 1, domain_occ[domain_counter]);
	    };
	};
	fprintf (stdout, "\n");
	
	fprintf (stdout, " Fit parameters:\n");
	fprintf (stdout, "  Number of fit runs: %ld\n", n_fit);
	fprintf (stdout, "  Variation of start parameters for each run: %f\n", par_var_width);
	fprintf (stdout, "  Maximum number of iterations: %ld\n", max_iteration);
	fprintf (stdout, "  Absolute stopping criteria (delta_abs): %.12f\n", delta_abs);
	fprintf (stdout, "  Relative stopping criteria (delta_rel): %.12f\n", delta_rel);
	fprintf (stdout, "  Seed for random number generator: %ld\n",rng_seed);
    };
    
    /* Simulated annealing */
    if (todo == sa_fit) {
	fprintf (stdout, "Simulated annealing fit of model to experimental structure factors.\n");
	fprintf (stdout, "   %s\n" , comment);

	fprintf (stdout, "Number of domains: %d \n\n", (int) nr_domains);
	fprintf (stdout, "Input file for experimental structure factors: %s\n", exp_struc_f_file);
	if (strlen(bulk_struc_file) > 0) {
	    fprintf (stdout, "Input file for bulk structure: %s\n", bulk_struc_file);
	};
	{
	    unsigned int domain_counter;
	    for (domain_counter = 0; domain_counter < nr_domains; ++domain_counter) {
		fprintf (stdout, "Fit Input file for domain %d: %s\n", domain_counter + 1, &(fit_struc_file[domain_counter][0]));
		fprintf (stdout, "Output file for domain %d: %s\n", domain_counter + 1, &(fit_coord_file[domain_counter][0]));
		fprintf (stdout, "Occupancy of domain %d: %f\n", domain_counter + 1, domain_occ[domain_counter]);
	    };
	};
	if (strlen(bulk_struc_file) > 0) {
	    fprintf (stdout, "Penetration depth: %.2f\n\n", penetration_depth);
	};
	fprintf (stdout, "\n");

	fprintf (stdout, " Simulated annealing parameters:\n");
	fprintf (stdout, "  Number of fit runs: %ld\n", n_fit);
	fprintf (stdout, "  Variation of start parameters for each run: %f\n", par_var_width);
	fprintf (stdout, "  Number of points to try for each step: %d\n", sa_par.n_tries);
	fprintf (stdout, "  Number of iterations at each temperature: %d\n", sa_par.iters_fixed_T);
	fprintf (stdout, "  Maximum step size: %f\n", sa_par.step_size);
	fprintf (stdout, "  Boltzmann K: %f\n", sa_par.k );
	fprintf (stdout, "  Initial temperatur: %f\n", sa_par.t_initial);
	fprintf (stdout, "  Damping factor for temperature mu_T: %f\n", sa_par.mu_t);
	fprintf (stdout, "  Minimum temperature: %f\n", sa_par.t_min);
	fprintf (stdout, "  Seed for random number generator: %ld\n",rng_seed);
    };
    

    fprintf(stdout, "---------------------------------------\n\n");
		        
    return;
}

/************************************************************************
 * void   read_input_file(char filename[])				*
 * 									*
 * Reads the input file and parses the input identifiers in it		*
 *									*
 * filename: Name of the input file					*
 ************************************************************************/

void read_input_file(char *filename)
{ 
    FILE *infile;				/* Input file */
    char line[LINE_LENGTH];	  	/* Line currently read */
    char *line_identifier_ptr;		/* Pointer to identifier in line */
    char *line_value_ptr;			/* Pointer to value in line */ 
    char *line_ptr;			/* Pointer into line */ 
    size_t identifier_length;		/* Length of the identifier in line */
    size_t value_length;			/* Length of the value in line */
    unsigned int line_counter;			/* Counts the lines read */
    char curr_identifier[LINE_LENGTH];	/* Current identifier read (first string in line) */
    char curr_value[LINE_LENGTH];		/* Current value read (second string in line) */


 

    infile = fopen(filename, "r");	/* Open input file */

    if (infile == NULL) {                     /* Check for error */
	fprintf(stderr, "Error: Could not open input file %s \n", filename);
	exit(8);
    };	

    /**********************************************
     * Loop for read and parsing the input file 	*
     **********************************************/

    for (line_counter = 0; line_counter < MAX_INPUT_FILE_LINES; ++line_counter)  {
	identifier_length = 0;
	value_length = 0;


	/*****  Read a line *****/
	fgets(line, (int) sizeof(line), infile);  
    
	if (feof(infile) != 0) {           /* End of file? */
	    break;
	};


	/***** Get the identifier and value *****/ 
	line_ptr = line;		/* Set line_ptr to start of line */
    
	if (*line_ptr == '\n') {		    /* Empty line? */ 
	    continue;
	};
    
	if (strlen(line) == 0) {		/* Should not happen !*/
	    fprintf (stderr, " Read empty line! Internal bug? ");
	    exit(8);
	};

	while(isspace(*line_ptr)) {		/* remove leading spaces and tabs */
	    ++line_ptr;
	};	
    
	line_identifier_ptr = line_ptr;     /* Set line_identifier_ptr to first non-whitespace character */
    
	if (*line_ptr == '#') {		    /* Comment line? */
	    continue;
	};

	while ((isspace(*line_ptr) == 0) && (*line_ptr != '=')) { 	    /* Find length of identifier */
	    ++line_ptr;
	    ++identifier_length;
	};
    
	if (*line_ptr == '\0') {			/* End of line before reading value? */
	    fprintf (stderr, "Error: Invalid line in input file: line %d\n", (line_counter+1));
	    fprintf (stderr, "%s %s %s\n.", line, curr_identifier, curr_value);
	    exit(8);
	};
	   
	/* remove whitespace and equal sign */
	while ((isspace(*line_ptr)) || (*line_ptr == '='))  {
	    ++line_ptr;
	};	   

	/* Set start of value */
	line_value_ptr = line_ptr;

	/* Find length of value */
	while ((isspace(*line_ptr) == 0) && (*line_ptr != '\0')) {
	    ++line_ptr;
	    ++value_length;
	};

	/* Check length of identifier and value */ 
	if ((identifier_length == 0) || (value_length == 0)) {
	    fprintf (stderr, "Error: Invalid line in input file: line %d\n", (line_counter+1));
	    fprintf (stderr, " %s %s %s\n.", line, curr_identifier, curr_value);
	    exit(8);
	};      

	/* Get identifier and value */		
	strncpy(curr_identifier, line_identifier_ptr, identifier_length);
	*(curr_identifier + identifier_length) = '\0';   

	strncpy(curr_value, line_value_ptr,  value_length);      
	*(curr_value + value_length) = '\0';   


	/** For comment identifier: read whole line */
	if (strcmp(curr_identifier, "comment") == 0) {
	    strcpy (comment, line_value_ptr);
	    continue;
	};	

	/***** Parse line *****/
	parse_line(curr_identifier, curr_value);    


    };   /* End loop for reading input file */ 

    fclose(infile);
  
    /**************************
     * Checks for consistency *
     **************************/ 

    /* was do found? */
    if (!do_read) {
	fprintf (stderr, "Please tell me what to do! (with \"do\" identifier in input file)");
	exit(8);
    };


    /* Is the number of files read same as number of domains */
    if ((todo == ls_fit) && (nr_fit_coord_files_read != nr_domains)) {
	fprintf (stderr, "Error: Number of output files for fitted coordinates different from number of domains!");
	exit(8);
    };

    /* Is the number of domains larger than MAX_DOMAINS */
    if (nr_domains > MAX_DOMAINS) {
        fprintf (stderr, "Error: Maximum number of domains is %d!\n", MAX_DOMAINS);
	fprintf (stderr, "Change MAX_DOMAINS and reading of structure files in source code and recompile.\n");
    };
/* Print input data */
    print_input();
  
    return;
}


/********************************************************************************
 * void structure_factor_range(void)						*
 * 										*
 * Calculates the structure factors for a range of h, k, l.			*
 ********************************************************************************/

void structure_factor_range(void)
{
    unsigned int h_counter;			/* Counter for calculating all structure factors */
    unsigned int k_counter;			/* Counter for calculating all structure factors */
    unsigned int l_counter;			/* Counter for calculating all structure factors */
    unsigned int nr_steps_h;			/* number of steps in h for calculating structure factors */
    unsigned int nr_steps_k;			/* number of steps in h for calculating structure factors */
    unsigned int nr_steps_l;			/* number of steps in h for calculating structure factors */
    double curr_h;			/* h currently used */
    double curr_k;			/* k currently used */
    double curr_l;			/* l currently used */
    double struc_factor;		/* Absolute value of structure factor */ 
    double test_h;			/* For testing whether bulk rod */
    double test_k;			/* For testing whether bulk rod */

    nr_struc_f_calc = 0;          /* First structure factor */

    /* Determine number of steps needed for h,k,l */
    nr_steps_h =  abs((int) floor((max_h - min_h) / step_h) + 1);
    nr_steps_k =  abs((int) floor((max_k - min_k) / step_k) + 1);
    nr_steps_l =  abs((int) floor((max_l - min_l) / step_l) + 1);
      
    for (h_counter = 0; h_counter < nr_steps_h; ++h_counter) {     /* Loop over h */
	curr_h = min_h + ((double) h_counter * step_h);  	/* Calculate current h */
	test_h =  curr_h / bulk_h;			       

	for (k_counter = 0; k_counter < nr_steps_k; ++k_counter) {      /* Loop over k */
	    curr_k = min_k + ((double) k_counter * step_k);  	/* Calculate current k */
	    test_k =  curr_k / bulk_k;

	    /* If designated as bulk rod -> do not calculate */
	    if (bulk_h != 12345678.12345678 && bulk_k != 12345678.12345678) {
	      if (((double) (((int)test_h) * bulk_h) == curr_h) && ((double) (((int)test_k) * bulk_k) == curr_k)) {
		continue;
	      };
	    };
		
	    for (l_counter = 0; l_counter < nr_steps_l; ++l_counter) {      /* Loop over k */

		curr_l = min_l + ((double) l_counter * step_l);  	/* Calculate current l */

		struc_factor = calc_structure_factor(atoms_unit_cell, atoms_bulk, lattice_par, penetration_depth, curr_h, curr_k, curr_l, nr_atoms, domain_occ, nr_domains, nr_atoms_bulk);

		struc_f_calc[nr_struc_f_calc].value = struc_factor * theo_scale;
		struc_f_calc[nr_struc_f_calc].h = curr_h;
		struc_f_calc[nr_struc_f_calc].k = curr_k;
		struc_f_calc[nr_struc_f_calc].l = curr_l;
		++nr_struc_f_calc;
	    };	 
	};
    }; 

    return;
}     


/************************************************************************************************
 * void structure_factor_data(struct structure_factor exp_struc_f[], unsigned int nr_struc_f))	*
 * 												*
 * Calculates the structure factors for the same h, k, l as in the experimental data.		*
 ************************************************************************************************/

void structure_factor_data(struct structure_factor exp_struc_f[], size_t nr_struc_f)
{
    unsigned int f_counter;			/* Counter over structure factors */
    double struc_factor;		/* Absolute value of structure factor */ 


	
    for (f_counter = 0; f_counter < nr_struc_f; f_counter++)
	{ 	

        struc_factor = calc_structure_factor(atoms_unit_cell, atoms_bulk,lattice_par, penetration_depth,  exp_struc_f[f_counter].h, exp_struc_f[f_counter].k, exp_struc_f[f_counter].l, nr_atoms, domain_occ, nr_domains, nr_atoms_bulk);
			
	    struc_f_calc[f_counter].h = exp_struc_f[f_counter].h;
	    struc_f_calc[f_counter].k = exp_struc_f[f_counter].k;
	    struc_f_calc[f_counter].l = exp_struc_f[f_counter].l;
	    struc_f_calc[f_counter].value = struc_factor * theo_scale;

	};	 
    return;
}    



/**************************************
 * Main                               *
 **************************************/

int main(int argc, char *argv[])
{
/* mtrace to catch memeory leaks */
/*     mtrace(); */


    /* Get the program name */
    program_name = argv[0];
  
  
    /******* Check for command line syntax *******/
    if (argc != 2) {                     /* Right number of arguments? */
	fprintf (stderr, "Please give file containing input variables as argument!\n");
	usage(program_name);                /* print usage information */
	exit(8);                /* and exit */
    };

    /* Get inputfile */
    input_file_name = argv[1];

    /********* Read and parse input file ********/  
    read_input_file(input_file_name);
  
  
  
    switch (todo) {			/* Decide what to do */
    
  
	/*****************
	 * Patterson map *
	 *****************/
    case patterson:
    
	/***** Read experimental structure factor data **********/
	nr_struc_f_exp = read_exp_struc_f(exp_struc_f_file);
    
	/**** Set symmetry matrices ***************/
	set_symmetry_matrix(symmetry);

	/***** Calculate Patterson map and write it to output file *******/
	patterson_map(patt_file_name, step_patterson_r1, step_patterson_r2, nr_struc_f_exp, comment); 
    
	break;

	
	/********************************
	 * Electron difference map	*
	 ********************************/
    case electron_difference_map:
    {
      /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
      /* !!!!!!!!!!! Bulk not correct (included in scale, but nowhere else !!!!!!!!!!!!! */
      /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

	unsigned int r1_counter;	/* Counts steps of the electron difference map in the r1 direction */
	unsigned int r2_counter;	/* Counts steps of the electron difference map in the r2 direction */
/*	unsigned int domain_counter; */	/* Counts over surface domains */
	double curr_r1;         /* Current r1 value */
	double curr_r2;		/* Current r2 value */
	double curr_value;      /* Current value of the electron difference map at r1,r2 */
	int number;             /* Number of charcters written */
	FILE *el_diff_file;	/* Output file for the electron difference map */
		

	/* WARNING */
	fprintf (stdout, "!!!   WARNING: calculation of electron difference map probably not correct !!! ");
	
	/***** Read experimental structure factor data **********/
	nr_struc_f_exp = read_exp_struc_f(exp_struc_f_file);

	/***** Read input structure files ******/
	{	
	    unsigned int struc_file_counter;		/* Counter for the input file */
	    for (struc_file_counter = 0; struc_file_counter < nr_domains; ++struc_file_counter) {
		nr_atoms[struc_file_counter] = read_struc(&(atoms_unit_cell[struc_file_counter][0]), &lattice_par, &(input_struc_file[struc_file_counter][0]));   
	    };
	    if (strlen(bulk_struc_file) > 0) {
	        nr_atoms_bulk = read_struc(&(atoms_bulk[0]), &lattice_par, &(bulk_struc_file[0]));   
	    }
	    else {
	        nr_atoms_bulk = 0;
	    };
	};
		     

	/**** Set symmetry matrices ***************/
	set_symmetry_matrix(symmetry);

	/* Open file for electron difference map */
	el_diff_file = fopen(el_diff_file_name,"w");
	if (el_diff_file == NULL) { 			/* Check for error */
	    fprintf(stderr, "Error: Could not open file %s for output of electron difference map!\n", el_diff_file_name);
	    exit(8);
	};

	/* Find the scale factor */
	scale = find_scale(struc_f_exp, nr_struc_f_exp, atoms_unit_cell, atoms_bulk, lattice_par, penetration_depth, nr_atoms, domain_occ, nr_domains, nr_atoms_bulk);
fprintf(stdout, "scale: %f\n",scale);
	number = fprintf(el_diff_file, comment);					/* Write comment in first line */
	/* Loops for calculating the electron difference map and writing it to output file */
	for (r1_counter=0; r1_counter < step_el_diff_r1; ++r1_counter) {   /* Loop over the r1 direction */
	    curr_r1 = (double)r1_counter / (double)step_el_diff_r1;		   /* Current r1 value */
	    for (r2_counter=0; r2_counter < step_el_diff_r2; ++r2_counter) {  	   /* Loop over the r2 direction */
		curr_r2 =  (double)r2_counter / (double)step_el_diff_r2;		   /* Current r2 value */

		/* Calculate value of the electron difference map at (curr_r1. curr_r2)  */
		    curr_value =  calc_electron_diff_map(curr_r1, curr_r2, struc_f_exp, nr_struc_f_exp, atoms_unit_cell, lattice_par, penetration_depth, nr_atoms, nr_domains, domain_occ, scale); 

/*   		curr_value = 0.0; 
 	        for (domain_counter = 0; domain_counter < nr_domains; domain_counter++) { */   /*	 Loop over domains */ 
/*		    curr_value += domain_occ[domain_counter] * calc_electron_diff_map(curr_r1, curr_r2, struc_f_exp, nr_struc_f_exp, &(atoms_unit_cell[domain_counter][0]), nr_atoms, scale); 
		    }; */

		/* Write current value of electron difference map to file */
		number = fprintf(el_diff_file, "%f %f %f\n", curr_r1 , curr_r2, curr_value);
		if (number ==-1) {							/* Check for error */
		    fprintf (stderr, "Error: Problem writing to file for electron difference map: %s !", el_diff_file_name);
		    exit(8);
		};
	    };
	};
		
	/* Close file for electron difference map */
	fclose(el_diff_file);
    };
	    
    break;

	
  	/************************
	 * Electron  map	*
	 ************************/
    case electron_map:
    {
	unsigned int r1_counter;	/* Counts steps of the electron difference map in the r1 direction */
	unsigned int r2_counter;	/* Counts steps of the electron difference map in the r2 direction */
/*	unsigned int domain_counter; */	/* Counts over surface domains */
	double curr_r1;         /* Current r1 value */
	double curr_r2;		/* Current r2 value */
	double curr_value;      /* Current value of the electron difference map at r1,r2 */
	double curr_h;          /* To make dummy structure factors */
	double curr_k;          /* To make dummy structure factors */ 
	unsigned int h_counter;			/* Counter for calculating all structure factors */
	unsigned int k_counter;			/* Counter for calculating all structure factors */
	unsigned int nr_steps_h;			/* number of steps in h for calculating structure factors */
	unsigned int nr_steps_k;			/* number of steps in h for calculating structure factors */
	unsigned int nr_steps_l;			/* number of steps in h for calculating structure factors */
	int number;             /* Number of charcters written */
	FILE *el_diff_file;	/* Output file for the electron difference map */
	double test_h;			/* For testing whether bulk rod */
	double test_k;			/* For testing whether bulk rod */
	

	/***** Read input structure files ******/
	{	
	    unsigned int struc_file_counter;		/* Counter for the input file */
	    for (struc_file_counter = 0; struc_file_counter < nr_domains; ++struc_file_counter) {
		nr_atoms[struc_file_counter] = read_struc(&(atoms_unit_cell[struc_file_counter][0]), &lattice_par, &(input_struc_file[struc_file_counter][0]));   
	if (strlen(bulk_struc_file) > 0) {
	  fprintf (stderr, "\nCAREFUL: Electron map does not use bulk.\n\n");
	};

	    };
	};

	/**** Set symmetry matrices ***************/
	set_symmetry_matrix(symmetry);

	/* Open file for electron difference map */
	el_diff_file = fopen(el_diff_file_name,"w");
	if (el_diff_file == NULL) { 			/* Check for error */
	    fprintf(stderr, "Error: Could not open file %s for output of electron difference map!\n", el_diff_file_name);
	    exit(8);
	};

	/* Set the scale factor */
	scale = 1.0;
	fprintf(stdout, "scale: %f\n",scale);

	/***** Make dummy structure factors ******/
	nr_struc_f_exp = 0;

	/* Determine number of steps needed for h,k,l */
	nr_steps_h = abs((int) (floor((max_h - min_h) / step_h))) + 1;
	nr_steps_k = abs((int) (floor((max_k - min_k) / step_k))) + 1;
	nr_steps_l = abs((int) (floor((max_l - min_l) / step_l))) + 1;
      
	for (h_counter = 0; h_counter < nr_steps_h; ++h_counter) {     /* Loop over h */
	  curr_h = min_h + ((double) h_counter * step_h);  	/* Calculate current h */
	  test_h =  curr_h / bulk_h;			       

	  for (k_counter = 0; k_counter < nr_steps_k; ++k_counter) {      /* Loop over k */
	    curr_k = min_k + ((double) k_counter * step_k);  	/* Calculate current k */
	    test_k =  curr_k / bulk_k;

	    /* If bulk rod -> do not calculate (l not taken into account)*/
	    if (((double) (((int)test_h) * bulk_h) == curr_h) && ((double) (((int)test_k) * bulk_k) == curr_k)) {
	      continue;
	    };

	    struc_f_exp[nr_struc_f_exp].value = 0.0;
	    struc_f_exp[nr_struc_f_exp].h = curr_h;
	    struc_f_exp[nr_struc_f_exp].k = curr_k;
	    struc_f_exp[nr_struc_f_exp].l = 0.0;
	    ++nr_struc_f_exp;
	  };
	};

      

	/************ Calculate electron map ***********/
	number = fprintf(el_diff_file, comment);					/* Write comment in first line */
	/* Loops for calculating the electron map and writing it to output file */
	for (r1_counter=0; r1_counter < step_el_diff_r1; ++r1_counter) {   /* Loop over the r1 direction */
	    curr_r1 = (double)r1_counter / (double)step_el_diff_r1;		   /* Current r1 value */
	    for (r2_counter=0; r2_counter < step_el_diff_r2; ++r2_counter) {  	   /* Loop over the r2 direction */
		curr_r2 =  (double)r2_counter / (double)step_el_diff_r2;		   /* Current r2 value */

		/* Calculate value of the electron difference map at (curr_r1. curr_r2)  (for experimental structure factors all 0 -> theoretical electron map */
		    curr_value =  calc_electron_diff_map(curr_r1, curr_r2, struc_f_exp, nr_struc_f_exp, atoms_unit_cell, lattice_par, penetration_depth, nr_atoms, nr_domains, domain_occ, scale); 

/*   		curr_value = 0.0; 
 	        for (domain_counter = 0; domain_counter < nr_domains; domain_counter++) { */   /*	 Loop over domains */ 
/*		    curr_value += domain_occ[domain_counter] * calc_electron_diff_map(curr_r1, curr_r2, struc_f_exp, nr_struc_f_exp, &(atoms_unit_cell[domain_counter][0]), nr_atoms, scale); 
		    }; */

		/* Write current value of electron map to file */
		number = fprintf(el_diff_file, "%f %f %f\n", curr_r1 , curr_r2, curr_value);
		if (number ==-1) {							/* Check for error */
		    fprintf (stderr, "Error: Problem writing to file for electron map: %s !", el_diff_file_name);
		    exit(8);
		};
	    };
	};
		
	/* Close file for electron map */
	fclose(el_diff_file);
    };
	    
    break;

  
    /******************************************************************
     * : Calculate theoretical structure factor for range of h, k, l. *
     ******************************************************************/
    case theo_range:
    
	/***** Read input structure files ******/
    {	
	unsigned int struc_file_counter;		/* Counter for the input file */
	for (struc_file_counter = 0; struc_file_counter < nr_domains; ++struc_file_counter) {
	    nr_atoms[struc_file_counter] = read_struc(&(atoms_unit_cell[struc_file_counter][0]), &lattice_par, &(input_struc_file[struc_file_counter][0]));   
	};
	if (strlen(bulk_struc_file) > 0) {
	    nr_atoms_bulk = read_struc(&(atoms_bulk[0]), &lattice_par, &(bulk_struc_file[0]));   
	}
	else {
	    nr_atoms_bulk = 0;
	};

    };
    
    /***** Calculate inplane structure factors ******/
    structure_factor_range();

         
    /***** Write the calculated structure factors to output file ******/
    write_structure_factor(calc_struc_f_file, nr_struc_f_calc, comment); 
    
    break;

    /*************************************************************
     * Calculate structure factors for the same range as in data *
     *************************************************************/
    
    case theo_data:

	/***** Read input structure files ******/
    {	
	unsigned int struc_file_counter;		/* Counter for the input structure file */
	for (struc_file_counter = 0; struc_file_counter < nr_domains; ++struc_file_counter) {
	    nr_atoms[struc_file_counter] = read_struc(&(atoms_unit_cell[struc_file_counter][0]), &lattice_par, &(input_struc_file[struc_file_counter][0]));   
	};
	if (strlen(bulk_struc_file) > 0) {
	    nr_atoms_bulk = read_struc(&(atoms_bulk[0]), &lattice_par, &(bulk_struc_file[0]));   
	}
	else {
	    nr_atoms_bulk = 0;
	};
    };

    /***** Read structure factor file ****/
    nr_struc_f_exp = read_exp_struc_f(exp_struc_f_file);
      
    /***** Calculate inplane structure factors ******/
    structure_factor_data(struc_f_exp ,nr_struc_f_exp);

    /***** Write the calculated structure factors to output file ******/
    write_structure_factor(calc_struc_f_file, nr_struc_f_exp, comment); 
       
    break;

    /**************************************************************
     * Structure refinement by least-squares fit to F		*
     **************************************************************/
    case ls_fit:
    
	/* Initialize starting values for parameters */
    { 
	unsigned int par_counter;	
	for (par_counter = 0; par_counter < MAX_FREE_PARAMETERS; par_counter++) {
	    start_par[par_counter] = 0.0;
	};
	
	memset(&fixed_parameters, FALSE, MAX_FREE_PARAMETERS);  /* Initialize to false = 0 */
	
    };	
	
/***** Read fit input files *****/
    {			
	unsigned int fit_file_counter;		/* Counter for the fit input file */
	for (fit_file_counter = 0; fit_file_counter < nr_domains; ++fit_file_counter) {
	    nr_atoms[fit_file_counter] = read_fit_file(&(atoms_unit_cell[fit_file_counter][0]), &lattice_par, &(atom_displacements[fit_file_counter][0]), &fixed_parameters[0], &start_par[0], &(fit_struc_file[fit_file_counter][0])); 
	};
    };
 
    /*** Read bulk structure ***/
    if (strlen(bulk_struc_file) > 0) {
        nr_atoms_bulk = read_struc(&(atoms_bulk[0]), &lattice_par, &(bulk_struc_file[0]));   
    }
    else {
        nr_atoms_bulk = 0;
    };
   
    /***** Read structure factor file ****/
    nr_struc_f_exp = read_exp_struc_f(exp_struc_f_file);
   
    /***** Least squares fit *******/
    do_ls_fit(atoms_unit_cell, atoms_bulk, lattice_par, penetration_depth, atom_displacements, fixed_parameters, start_par, domain_occ,  nr_atoms, nr_domains, nr_atoms_bulk, struc_f_exp, nr_struc_f_exp, max_iteration, delta_abs, delta_rel, end_positions, rng_seed, print_intermediate, n_fit, par_var_width);	
		
    /***** Print the fitted coordinates to file ******/
    {
	unsigned int domain_nr;
	for (domain_nr = 0; domain_nr < nr_domains; domain_nr++) {
	    write_fitted_coordinates (&(end_positions[domain_nr][0]), nr_atoms[domain_nr], lattice_par, comment, &(fit_coord_file[domain_nr][0]));
	};
    };
    break;
   
    
     /**************************************************************
     * Structure refinement by simulated annealing		*
     **************************************************************/
    case sa_fit:
    
	/* Initialize starting values for parameters */
    { 
	unsigned int par_counter;	
	for (par_counter = 0; par_counter < MAX_FREE_PARAMETERS; par_counter++) {
	    start_par[par_counter] = 0.0;
	};
    };	
	
/***** Read fit input files *****/
    {			
	unsigned int fit_file_counter;		/* Counter for the fit input file */
	for (fit_file_counter = 0; fit_file_counter < nr_domains; ++fit_file_counter) {
	    nr_atoms[fit_file_counter] = read_fit_file(&(atoms_unit_cell[fit_file_counter][0]), &lattice_par, &(atom_displacements[fit_file_counter][0]), &fixed_parameters[0], &start_par[0], &(fit_struc_file[fit_file_counter][0])); 
	};
    };
    
    /*** Read bulk structure ***/
    if (strlen(bulk_struc_file) > 0) {
        nr_atoms_bulk = read_struc(&(atoms_bulk[0]), &lattice_par, &(bulk_struc_file[0]));   
    }
    else {
        nr_atoms_bulk = 0;
    };

    /***** Read structure factor file ****/
    nr_struc_f_exp = read_exp_struc_f(exp_struc_f_file);
   
    /***** Simulated annealing fit *******/
    do_sa_fit(atoms_unit_cell, atoms_bulk, lattice_par, penetration_depth, atom_displacements, fixed_parameters, start_par, domain_occ,  nr_atoms, nr_domains,nr_atoms_bulk,  struc_f_exp, nr_struc_f_exp, sa_par, end_positions, rng_seed, print_intermediate, n_fit, par_var_width);	
		
    /***** Print the fitted coordinates to file ******/
    {
	unsigned int domain_nr;
	for (domain_nr = 0; domain_nr < nr_domains; domain_nr++) {
	    write_fitted_coordinates (&(end_positions[domain_nr][0]), nr_atoms[domain_nr], lattice_par, comment, &(fit_coord_file[domain_nr][0]));
	};
    };
    break;

   
    default:		
	/* Should never get here! */
	fprintf(stderr, "Internal error! Reached default value in switch!\n");
	fprintf(stdout, "Internal error! Reached default value in switch!\n");
	exit(8);
	break;
    };
  
    
    return(0);
}


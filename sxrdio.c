#include "sxrddefs.h"
#include "sxrdio.h"

/* Function prototypes */
int find_el_nr(char atom_name[20]);
void find_el_name(int el_nr, char *atom_name);

/********************************************************************************************************
 * unsigned int read_struc( struct atom_structure (*read_atom)[MAX_ATOMS], struct latt_par_structure *lattice_par, char filename[])                      *
 *                                                                                                      *
 * Reads structure file into lattice_par and atoms                                              *
 * Returns number of atoms read.                                                                        *
 *                                                                                                      *
 * filename: Name of structure file                                                                     *
 *                                                                                                      *
 * Format of structure file:                                                                            *
 * 1st line: Any comment. is ignored                                                                    *
 * 2nd line: lattice paramters and angles: a b c alpha beta gamma                                       *
 * other lines: Atomic positions: Atom name x y z                                                       *
 * Atom name is two characters abbrevation, x,y,z in fractional coordinates.                            *
 ********************************************************************************************************/

unsigned int read_struc(struct atom_structure *read_atom, struct latt_par_structure *lattice_par, char filename[])
{
    int number;                                   /* Number of variables read */
    int i;                /* counter */
    unsigned int atom_j;                                           /* Number of atom currently read */
    int line_nr;             /* Line number */
    char atom_name[20];				  /* Name of atom */
    char comment_line[LINE_LENGTH];			/* line with comment (is ignored) */
    char line[LINE_LENGTH];             /* input line */ 
    FILE *struc_file;                                     /* File containing the input structure data */
    double dummy[20];

    struc_file = fopen(filename,"r");             /* Open file with input structure */

    if (struc_file == NULL) {                     /* Check for error */
	printf("Error: Could not open structure file %s \n", filename);
	fclose(struc_file);                   /* Close file with input structure */
        exit(8);
    };

    fgets(comment_line, (int) sizeof(comment_line), struc_file);           /* Discard first line used for comments */
    number = fscanf(struc_file, "%lf %lf %lf %lf %lf %lf %lf", &lattice_par->a, &lattice_par->b, &lattice_par->c, &lattice_par->alpha, &lattice_par->beta, &lattice_par->gamma, dummy) ;    /* Read lattice parameters and angles. */
    if (number != 6) {
	fprintf(stderr, "Error: Input model structure file format wrong (lattice parameters)!\n");      /* wrong number of variables read? */
	fprintf(stdout, "Error: Input model structure file format wrong (lattice parameters)!\n");      /* wrong number of variables read? */
	fclose(struc_file);                   /* Close file with input structure */
	exit(8);
    };
    if (feof(struc_file) != 0) {           /* End of file? */
        atom_j = 0;	/* Neccessary to get right number of atoms */
	fclose(struc_file);                   /* Close file with input structure */
	
	return atom_j ;                                /* Return number of atoms read */
    };				

    line_nr = 2;
    atom_j = 0;
    while(1) {
      line_nr++;

      fgets(line, (int) sizeof(line), struc_file);
      if (feof(struc_file)) {           /* End of file? */
	break;			
      };
      if (ferror(struc_file)) {
	fprintf(stderr, "Error while reading file %s, line %d", filename, line_nr);
      };
      
      if (strlen(line) == 0 || strcmp(line, "\n") == 0 || strncmp(line, "#", 1) == 0) {    /*   empty line or comment line */
	continue;
      };
      
      if (strncmp (line, "displ", 5) == 0 || strncmp (line, "dw_par", 5) == 0 || strncmp (line, "occ", 3) == 0 || strncmp (line, "aniso_dw_par ", 12) == 0 || strncmp (line, "start_par ", 10) == 0|| strncmp (line, "fix ", 4) == 0 ) {  /* Fit keywords - ignore */
	continue;
      }

      else if  (strncmp (line, "aniso_dw ", 9) == 0) {       /* Anisotropic DW parameter */
	number = sscanf(line, "aniso_dw %lf %lf %lf %lf %lf %lf", &(read_atom[atom_j-1].aniso_dw[0]),  &(read_atom[atom_j-1].aniso_dw[1]), &(read_atom[atom_j-1].aniso_dw[2]),  &(read_atom[atom_j-1].aniso_dw[3]),   &(read_atom[atom_j-1].aniso_dw[4]),  &(read_atom[atom_j-1].aniso_dw[5]));  /*  (b11,b22,b33,b12,b13,23) */
	if (number != 6) {
	  fprintf (stderr, "Error: Need six values for anisotropic DW factor in line %d, file %s", line_nr, filename);
	  fclose(struc_file);
	  exit(8);
	};
	continue;
      }

      else if (strncmp(line, "pos", 3) == 0 ) {     /* pos keyword - atomic positions */
	number = sscanf(line, "pos %s %lf %lf %lf %lf %lf", atom_name, &(read_atom[atom_j].x), &(read_atom[atom_j].y), &(read_atom[atom_j].z), &(read_atom[atom_j].dw_par), &(read_atom[atom_j].occ));
	if (number == 5) {	/* No Occupancy */
	    read_atom[atom_j].occ = 1.0;
	};		
	atom_j++;
      }

      else {         /* No keyword - atomic positions */
	number = sscanf(line, "%s %lf %lf %lf %lf %lf", atom_name, &(read_atom[atom_j].x), &(read_atom[atom_j].y), &(read_atom[atom_j].z), &(read_atom[atom_j].dw_par), &(read_atom[atom_j].occ));
	if (number == 5) {	/* No Occupancy */
	  read_atom[atom_j].occ = 1.0;
	};		
	atom_j++;
	/*	fprintf("AAAAAA  %s  %f  %f  %f  %f  %f", atom_name, (read_atom[atom_j].x), (read_atom[atom_j].y), (read_atom[atom_j].z), (read_atom[atom_j].dw_par), (read_atom[atom_j].occ));*/
	
      };

      if (number != 5 && number != 6 && atom_name[0] != '#') {
	fprintf(stderr, "Error: Input model structure file %s format wrong (atomic data) in line %d!\n%d", filename, line_nr, number);    /* wrong number of variables read? */
	fprintf(stdout, "Error: Input model structure file %s format wrong (atomic data) in line %d!\n", filename, line_nr);    /* wrong number of variables read? */
	fclose(struc_file);                   /* Close file with input structure */
	exit(8); 
      };
      
      for (i=0; i<6; i++) {     /* Set anisotropic DW to 0 */
	read_atom[atom_j-1].aniso_dw[i] = 0.0;
      };
      
      read_atom[atom_j-1].el_nr = find_el_nr(atom_name);  	/* Assign element number from atom name */
    
      if (atom_j == MAX_ATOMS) {      /* Maximum number of atoms reached */
	fprintf(stderr, "ERROR: Maximum number of atoms is %d!", MAX_ATOMS);
	fclose(struc_file);                   /* Close file with input structure */
	exit(8); 
      };

    };

    fclose(struc_file);                   /* Close file with input structure */

    return atom_j ;                                /* Return number of atoms read */
}

 
/********************************************************************************************************
 * unsigned int read_exp_struc_f(char filename[])								*
 * 													*
 * Reads file with experimental structure factors into struc_f_exp					*
 * Returns number of structure factors read								*
 *													*
 * filename: Name of file with experimental structure factors 						*
 * 													*
 * Format of file:									 		*
 * other lines: h   k   l   structure factor   sigma							*
 *													*
 * sigma is the error of the measured  structure factors						*
 ********************************************************************************************************/
unsigned int read_exp_struc_f(char filename[])
{
    int number; 					/* Number of variables read */
    unsigned int struc_f_j;					/* Number of structure factor currently read */
    FILE *exp_struc_f_file;				/* File containing the experimental structure factors data */
  
    exp_struc_f_file = fopen(filename,"r");		/* Open file with experimental structure factors */
    if (exp_struc_f_file == NULL) { 			/* Check for error */
	fprintf(stderr, "Error: Could not open file with experimental structure factors %s \n", filename);
	fprintf(stdout, "Error: Could not open file with experimental structure factors %s \n", filename);
	exit(8);
    };
    
    /*  number = fscanf(exp_struc_f_file, "\n");			*/	/* Discard first line used for comments */ 
  
    for (struc_f_j = 0; struc_f_j<MAX_STRUCTURE_FACTORS; struc_f_j++) {   	/* Read h,k,l, structure factor and sigma */
	number = fscanf(exp_struc_f_file, "%lf %lf %lf %lf %lf", &(struc_f_exp[struc_f_j].h), &(struc_f_exp[struc_f_j].k), &(struc_f_exp[struc_f_j].l), &(struc_f_exp[struc_f_j].value), &(struc_f_exp[struc_f_j].sigma));
    
	if (feof(exp_struc_f_file) != 0) {           /* End of file */
	    break;
	};
	if (number == -1) {     /* Error? */
      
	    fprintf(stderr, "Some problem with file with experimental structure factors %s !\n ", filename);
	    fprintf(stdout, "Some problem with file with experimental structure factors %s !\n ", filename);
	    fclose(exp_struc_f_file); 			/* Close file with experimental structure factors */
	    exit(8);
	};
	if (number != 5) {
	    fprintf(stderr, "input format of file with experimental structure factors is wrong! Continuing.\n");	/* wrong number of variables read? */
	    fprintf(stdout, "input format of file with experimental structure factors is wrong! Continuing.\n");	/* wrong number of variables read? */
	}; 
    };     
  
    fclose(exp_struc_f_file); 			/* Close file with experimental structure factors */
  
    return struc_f_j ;				/* Return the number of structure factors read */
}


/********************************************************************************************************
 * int read_fit_file(struct atom_structure (*read_atom)[MAX_ATOMS], struct latt_par_structure *lattice_par, struct displacements (*read_displ)[MAX_ATOMS], int (*fixed_parameters[], double *start_par, char filename[])                      *
 *                                                                                                      *
 * Reads input file for structure refinement into lattice_par, read_atom, read_displ                    *
 * Returns number of atoms read.                                                                        *
 *                                                                                                      *
 * filename: Name of fit input files                                                                    *
 *                                                                                                      *
 * Format of fit input file:                                                                            *
 * 1st line: Any comment. is ignored                                                                    *
 * 2nd line: lattice paramters and angles: a b c alpha beta gamma                                       *
 * other lines:												*
 * Atomic positions: 											*
 * 	Start with keyword "pos", then:   Atom name   x   y   z   debye waller parameter                *
 * 	Atom name is two characters abbrevation, x,y,z in fractional coordinates.                       *
 * Displacement vectors:										*
 * 	Start with keyword "displ1", "displ2" or "displ3", 						*
 * 	then: Nr of displacement parameter   dx   dy   dz		 				*
 * Debye-Waller factor:											*
 * 	Start with keyword "dw_par",									* 
 * 	then: Nr of parameter scale  	      		 						*
 * Fixed parameters: Parameters can be fixed by the keyword "fix" and the Nr of the parameter.		*
 * 													*
 ********************************************************************************************************/

int read_fit_file(struct atom_structure *read_atom, struct latt_par_structure *lattice_par,  struct displacements *read_displ, int *fixed_parameters, double *start_par, char filename[])
{
    int number;                                   /* Number of variables read */
    int line_counter;				/* Counter for lines in input file */
    int atom_number = -1;				/* Number of atom currently read */
    int atom_counter;
    int i;                         /* Counter */
    char comment_line[LINE_LENGTH];			/* line with comment (is ignored) */
    /* Variables for reading one line */
    char keyword[LINE_LENGTH];				/* pos: Atomic positions and DW-parameter, displ: Displacements, dw_par: Parameter for Debye-Waller, occ: Occupancy */ 
    char string[LINE_LENGTH];				/* either atom name or number of fit parameter.  */
    double number1;					/* x-position or x-displacement */
    double number2;					/* y-position or y-displacement */
    double number3;					/* z-position or z-displacement */
    double number4;					/* Debye-Waller parameter or empty */
    double number5;					/* Occupancy or empty */
    double number6;					/* empty or beta_23*/

    FILE *fit_file;                                     /* File containing the input data for fit */

    /* Initialize displacement structure */
    for(atom_counter = 0; atom_counter < MAX_ATOMS; atom_counter++) {
	read_displ[atom_counter].nr1 = -1;
	read_displ[atom_counter].nr2 = -1;
	read_displ[atom_counter].nr3 = -1;	
	read_displ[atom_counter].nr_dw = -1;
	read_displ[atom_counter].nr_occ = -1;
	for (i=0; i<6; i++) {
	  read_displ[atom_counter].nr_aniso_dw[i] = -1;
	};
    };
    
  
    fit_file = fopen(filename,"r");             /* Open file with fit data */

    if (fit_file == NULL) {                     /* Check for error */
	printf("Error: Could not open input file for fit %s \n", filename);
        exit(8);
    };

    fgets(comment_line, (int) sizeof(comment_line), fit_file);           /* Discard first line used for comments */
    number = fscanf(fit_file, "%lf %lf %lf %lf %lf %lf", &lattice_par->a, &lattice_par->b, &lattice_par->c, &lattice_par->alpha, &lattice_par->beta, &lattice_par->gamma) ;    /* Read lattice parameters and angles. */
    if (number != 6) {
	fprintf(stderr, "Error:  Format of input file for fit %s wrong (lattice parameters)!\n", filename);      /* wrong number of variables read? */
	fprintf(stdout, "Error:  Format of input file for fit %s wrong (lattice parameters)!\n", filename);      /* wrong number of variables read? */
	fclose(fit_file);                   /* Close file with fit data */
	exit(8);
    };

    /* Read atom data or displacements */
    for (line_counter = 0; line_counter < MAX_ATOMS; line_counter++) {               /* Read atom names and coordinates or displacements */
	strcpy(keyword, ""); 
	number = fscanf(fit_file, "%s %s %lf %lf %lf %lf %lf %lf", keyword, string, &number1, &number2, &number3, &number4, &number5, &number6);
	
	if (keyword[0] != '#') {	    /* Comment ? */
	           
	    if (strcmp(keyword, "pos") == 0) { 			/* Positions, DW parameter and occupancy? */
		if (number != 6 && number !=7) {
		    fprintf(stderr, "Error: Format of input file for fit %s wrong in line %d (atomic data)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fprintf(stdout, "Error: Format of input file for fit %s wrong in line %d (atomic data)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fclose(fit_file);                   /* Close file with fit data */
		    exit(8);
		};
	    
		atom_number++;
		(read_atom + atom_number)->el_nr = find_el_nr(string);  /* Assign element number from atom name */
		(read_atom + atom_number)->x = number1;
		(read_atom + atom_number)->y = number2;
		(read_atom + atom_number)->z = number3;
		(read_atom + atom_number)->dw_par = number4;
		if (number == 6) {		/* No occupancy */
		    (read_atom + atom_number)->occ = 1.0;
		};
		if (number == 7) {		/* With occupancy */
		    (read_atom + atom_number)->occ = number5;
		};
		/* set anisotropic Debyer-Waller factor to 0*/
		for (i=0; i<6; i++) { 
		  (read_atom + atom_number)->aniso_dw[i] = 0.0;
		};
	    };

	    if (strcmp(keyword, "aniso_dw") == 0) { 			/* Anisotropic Debye-Waller factor ? */
	      if (number != 7) {
		    fprintf(stderr, "Error: Format of  input file for fit %s wrong in line %d (anisotropic Debye-Waller factor)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fprintf(stdout, "Error: Format of  input file for fit %s wrong in line %d (Anisotropic Debye-Waller factor)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fclose(fit_file);                   /* Close file with fit data */
		    exit(8);
		};
	      (read_atom + atom_number)->aniso_dw[0] = strtod(string, NULL);
	      (read_atom + atom_number)->aniso_dw[1] = number1;
	      (read_atom + atom_number)->aniso_dw[2] = number2;
	      (read_atom + atom_number)->aniso_dw[3] = number3;
	      (read_atom + atom_number)->aniso_dw[4] = number4;
	      (read_atom + atom_number)->aniso_dw[5] = number5;
	    };
  

	    if (strcmp(keyword, "displ1") == 0) { 			/* Displacement parameters ? */
		if (number != 5) {
		    fprintf(stderr, "Error: Format of  input file for fit %s wrong in line %d (displacement parameters)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fprintf(stdout, "Error: Format of  input file for fit %s wrong in line %d (displacement parameters)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fclose(fit_file);                   /* Close file with fit data */
		    exit(8);
		};
		(read_displ + atom_number)->nr1 = atoi(string);
		(read_displ + atom_number)->dx1 = number1;
		(read_displ + atom_number)->dy1 = number2;
		(read_displ + atom_number)->dz1 = number3;
	    };

	    if (strcmp(keyword, "displ2") == 0) { 			/* Displacement parameters ? */
		if (number != 5) {
		    fprintf(stderr, "Error: Format wrong of input file for fit %s in line %d (displacement parameters)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fprintf(stdout, "Error: Format wrong of input file for fit %s in line %d (displacement parameters)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fclose(fit_file);                   /* Close file with fit data */
		    exit(8);
		};
		(read_displ + atom_number)->nr2 = atoi (string);
		(read_displ + atom_number)->dx2 = number1;
		(read_displ + atom_number)->dy2 = number2;
		(read_displ + atom_number)->dz2 = number3;
	    };
	    
	    if (strcmp(keyword, "displ3") == 0) { 			/* Displacement parameters ? */
		if (number != 5) {
		    fprintf(stderr, "Error: Format of input file for fit %s wrong in line %d (displacement parameters)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fprintf(stdout, "Error: Format of input file for fit %s wrong in line %d (displacement parameters)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fclose(fit_file);                   /* Close file with fit data */		
		    exit(8);
		};
		(read_displ + atom_number)->nr3 = atoi(string);
		(read_displ + atom_number)->dx3 = number1;
		(read_displ + atom_number)->dy3 = number2;
		(read_displ + atom_number)->dz3 = number3;
	    };

	    if (strcmp(keyword, "dw_par") == 0) { 			/* Parameter for Debye-Waller factor ? */
		if (number != 3) {
		    fprintf(stderr, "Error: Format of input file for fit %s wrong in line %d (Debye-Waller fit parameter)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fprintf(stdout, "Error: Format of input file for fit %s wrong in line %d (Debye-Waller fit parameter)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fclose(fit_file);                   /* Close file with fit data */
		    exit(8);
		};
		(read_displ + atom_number)->nr_dw = atoi(string);
		(read_displ + atom_number)->dw_scale = number1;
	    };

	    if (strcmp(keyword, "aniso_dw_par") == 0) { 			/* Parameter for anisotropic Debye-Waller factor ? */
		if (number != 8) {
		  fprintf(stderr, "Error: Format of input file for fit %s wrong in line %d (anisotropic Debye-Waller fit parameter)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		  fprintf(stdout, "Error: Format of input file for fit %s wrong in line %d (anisotropic Debye-Waller fit parameter)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		  fclose(fit_file);                   /* Close file with fit data */
		  exit(8);
		};
		for (i=0; i<6; i++) {     /* look for empty parameter and assign */
		  if 	((read_displ + atom_number)->nr_aniso_dw[i] == -1) {
		    (read_displ + atom_number)->nr_aniso_dw[i] = atoi(string);
		    (read_displ + atom_number)->aniso_dw_scale[i][0] = number1;
		    (read_displ + atom_number)->aniso_dw_scale[i][1] = number2;
		    (read_displ + atom_number)->aniso_dw_scale[i][2] = number3;
		    (read_displ + atom_number)->aniso_dw_scale[i][3] = number4;
		    (read_displ + atom_number)->aniso_dw_scale[i][4] = number5;
		    (read_displ + atom_number)->aniso_dw_scale[i][5] = number6;
		    break;
		  };
		};
	    };
	    
	    if (strcmp(keyword, "occ") == 0) { 			/* Parameter occupancy? */
		if (number != 2) {
		    fprintf(stderr, "Error: Format of input file for fit %s wrong in line %d (Occupancy fit parameter)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fprintf(stdout, "Error: Format of input file for fit %s wrong in line %d (Occupancy fit parameter)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fclose(fit_file);                   /* Close file with fit data */
		    exit(8);
		};
		(read_displ + atom_number)->nr_occ = atoi(string);
	    };


	    if (strcmp(keyword, "start_par") == 0) { 			/* Start values for parameters ? */
		if (number != 3) {
		    fprintf(stderr, "Error: Format of input file for fit %s wrong in line %d (starting value for parameters)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fprintf(stdout, "Error: Format of input file for fit %s wrong in line %d (starting value for parameters)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fclose(fit_file);                   /* Close file with fit data */
		    exit(8);
		};
		*(start_par + atoi(string))  = number1;
	    };

	    if (strcmp(keyword, "fix") == 0) {		/* Number of fixed parameter? */
		if (number != 2) {
		    fprintf(stderr, "Error: Format of input file for fit %s wrong in line %d (fixed parameters)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fprintf(stdout, "Error: Format of input file for fit %s wrong in line %d (fixed parameters)! Wrong number of columns!\n", filename, line_counter + 2);    /* wrong number of variables read? */
		    fclose(fit_file);                   /* Close file with fit data */
		    exit(8);
		};
		*(fixed_parameters + atoi(string)) = TRUE;
	    };

	};
	if (feof(fit_file) != 0) {           /* End of file? */
	    break;
	};
    
	if (strcmp(keyword, "displ1") != 0 && strcmp(keyword, "displ2") != 0 &&  strcmp(keyword, "displ3") != 0 && strcmp(keyword, "pos") != 0  && strcmp(keyword, "aniso_dw") != 0 && strcmp(keyword, "start_par") != 0 && strcmp(keyword, "dw_par") != 0 && strcmp(keyword, "aniso_dw_par") != 0 && strcmp(keyword, "occ") != 0 && strcmp(keyword, "fix") != 0 && keyword[0] != '#') { 	/* Don't understand keyword? */
	    fprintf(stderr, "Error: Format of input file %s for fit wrong in line %d (Line doesn't start with pos, displ1, displ2, displ3, dw_par, occ or start_par)\n", filename, line_counter + 2);
	    fprintf(stdout, "Error: Format of input file %s for fit wrong in line %d (Line doesn't start with pos, displ1, displ2, displ3, dw_par, occ or start_par)\n", filename, line_counter + 2);
	    fprintf(stderr, "%s %s %f %f %f %f\n", keyword, string, number1, number2, number3, number4);
	    fclose(fit_file);                   /* Close file with fit data */	
	    exit(8);
	};
  
    };

    fclose(fit_file);                   /* Close file with fit data */

    return atom_number + 1;
}


/************************************************************************************************
 * int write_structure_factor(char filename[], size_t nr_f, char write_comment)			*
 *												*
 * Write the calculated structure factors to the file given by filename				*
 * 												*
 * filename: Name of output file								*
 * nr_f: Number of structure factors								*
 *												*
 * Format of output file:									*
 * First line:											*
 * Other lines: h  k  l  structure factor.							*
 ************************************************************************************************/

void write_structure_factor(char filename[], size_t nr_f, char write_comment[])
{
    int number; 					/* Number of characters written */
    unsigned int f_nr;						/* Number of structure factor currently written */
    FILE *calc_f_file;					/* Output file for the calculated structure factors */
  
    /* Open file for output */
    calc_f_file = fopen(filename,"w");		
    if (calc_f_file == NULL) { 			/* Check for error */
	fprintf(stderr, "Error: Could not open file %s for output of calculated structure factors!\n", filename);
	fprintf(stdout, "Error: Could not open file %s for output of calculated structure factors!\n", filename);
	exit(8);
    };
  
    /* Write data to file */
    number = fprintf(calc_f_file, write_comment);					/* Write comment in first line */
    if (number ==-1) {							/* Check for error */
	fprintf (stderr, "Error: Problem writing to file for calculated structure factors: %s !", filename);
	fprintf (stdout, "Error: Problem writing to file for calculated structure factors: %s !", filename);
	fclose(calc_f_file);
	exit(8);
    };
  
    for (f_nr=0; f_nr<nr_f; ++f_nr) {
	number = fprintf(calc_f_file, "%f %f %f %f 0.1 \n", struc_f_calc[f_nr].h, struc_f_calc[f_nr].k, struc_f_calc[f_nr].l, struc_f_calc[f_nr].value);
	if (number ==-1) {							/* Check for error */
	    fprintf (stderr, "Error: Problem writing to file for calculated structure factors: %s !", filename);
	    fprintf (stdout, "Error: Problem writing to file for calculated structure factors: %s !", filename);
	    fclose(calc_f_file);
	    exit(8);
	};
    };   
    /* Close output file */ 
    fclose(calc_f_file);
    return;
}


/****************************************************************************************************************
 * write_fitted_coordinates (struct atom_structure positions, size_t nr_atoms, 					*
 * struct latt_par_structure lattice_par. char *filename[])							*
 *	 													*
 * Write the fitted coordinates to file given by filename							*
 *														*
 * positions: fitted positions											*
 * nr_atoms: Number of atoms											*
 * filename: Filenames for different domains									*
 ****************************************************************************************************************/

void write_fitted_coordinates (struct atom_structure *positions,  size_t nr_atoms, struct latt_par_structure lattice_par, char write_comment[], char *filename)
{
    char atom_name[20];			/* Name of atom */
    int number; 					/* Number of characters written */
    unsigned int atom_nr;						/* Number of atom currently written */
    int aniso_counter;
    int write_aniso;        /* Boolean whether anisotropic DW factor should be written */
    FILE *fit_coord_file;					/* Output file for the calculated structure factors */
  
    /* Open file for output */
    fit_coord_file = fopen(filename,"w");		
    if (fit_coord_file == NULL) { 			/* Check for error */
	fprintf(stderr, "Error: Could not open file %s for output of fitted coordinates!\n", filename);
	fprintf(stdout, "Error: Could not open file %s for output of fitted coordinates!\n", filename);
	exit(8);
    };
  
    /* Write data to file */
    number = fprintf(fit_coord_file, write_comment);					/* Write comment in first line */
    if (number ==-1) {							/* Check for error */
	fprintf (stderr, "Error: Problem writing to file for output of fitted coordinates %s !", filename);
	fprintf (stdout, "Error: Problem writing to file for output of fitted coordinates %s !", filename);
	fclose(fit_coord_file);
	exit(8);
    };
    number = fprintf(fit_coord_file, "%f %f %f %f %f %f\n", lattice_par.a, lattice_par.b, lattice_par.c, lattice_par.alpha, lattice_par.beta, lattice_par.gamma ); 	/* Write lattice parameter and angles */
     if (number ==-1) {							/* Check for error */
	fprintf (stderr, "Error: Problem writing to file for output of fitted coordinates %s !", filename);
	fprintf (stdout, "Error: Problem writing to file for output of fitted coordinates %s !", filename);
	fclose(fit_coord_file);
	exit(8);
    };

    for (atom_nr = 0; atom_nr < nr_atoms; ++atom_nr) {
	find_el_name((positions + atom_nr)->el_nr, atom_name);     /* Assign atom name from element number */
	number = fprintf(fit_coord_file, "%s %f %f %f %f %f \n", atom_name, (positions + atom_nr)->x, (positions + atom_nr)->y, (positions + atom_nr)->z, (positions + atom_nr)->dw_par, (positions + atom_nr)->occ);
	if (number == -1) {							/* Check for error */
	    fprintf (stderr, "Error: Problem writing to file for fitted coordinates: %s !", filename);
	    fprintf (stdout, "Error: Problem writing to file for fitted coordinates: %s !", filename);
	    fclose(fit_coord_file);
	    exit(8);
	};
	/* Does atom have anisotropic Debye-Waller factor? => Write out in extra line */
	write_aniso = FALSE;
	for (aniso_counter=0; aniso_counter<6; aniso_counter++) {
	  if (((positions + atom_nr)->aniso_dw)[aniso_counter] != 0.0) {
	    write_aniso = TRUE;
	  };
	};
	if (write_aniso == TRUE) {
	  fprintf(fit_coord_file, "aniso_dw %f %f %f %f %f %f\n", ((positions + atom_nr)->aniso_dw)[0],  ((positions + atom_nr)->aniso_dw)[1],  ((positions + atom_nr)->aniso_dw)[2],  ((positions + atom_nr)->aniso_dw)[3],  ((positions + atom_nr)->aniso_dw)[4],  ((positions + atom_nr)->aniso_dw)[5]);
	};
    };   

    /* Close output file */ 
    fclose(fit_coord_file);

    return;
}



	


/*************************************************************************************
 * void usage(char name)                                                             * 
 *                                                                                   *
 * Print information how to use the program and exit                                 *
 *************************************************************************************/

void usage(char name[])
{ 
    fprintf(stderr, "  %s  -  Ver. %s,  %s\n", name, VERSION, VERSION_DATE);
    fprintf(stderr, "    Calculates structure factors from model structure, \n");
    fprintf(stderr, "    Patterson map or electron difference map, or fits \n");
    fprintf(stderr, "    a model structure to experimental structure factor. \n");  
    fprintf(stderr, "  Usage: %s inputfile\n", name);
    fprintf(stderr, "  The input file should contain the input variables.\n");

    exit(8);
}

/*************************************************************************************
 *      int find_el_nr(char atom_name[20])                                           * 
 *                                                                                   *
 * Returns the element number of the atom with name atom_name                        *
 *************************************************************************************/

int find_el_nr(char atom_name[21]) 
{
#define NR_ELEMENTS 19 		/* element number of highest defined element */
    const char el_name[NR_ELEMENTS][20] =  {"H",   /* List of atomic names sorted by element (atomic) number */
					    "He",  /* Starting by 0 for H */ 
					    "Li", 
					    "Be", 
					    "B",
					    "C",
					    "N",
					    "O",					   
				/*	    "F",
					    "Ne",
					    "Na",
					    "Mg",
					    "Al",	*/
					    "Si",
					    "Si_val",
					    "Si4+",
					    "O-",
					    "P", 
					    "In",
					    "In3+",
					    "Au",
					    "Mn",
					    "Fe",
					    "Pb"
    };

    
    int counter;
    int el_nr = -1;

    for (counter = 0; counter <NR_ELEMENTS; counter++) {
	if (strcmp(atom_name, &(el_name[counter][0])) == 0) {
	    el_nr = counter + 1;
	};
    };

    if (el_nr == -1) {
	fprintf (stderr, "Error: Atomic scattering factor data not available for element %s!\n", atom_name);
    };

    return el_nr; 
}


/*************************************************************************************
 *  void find_el_name](int el_nr, char el_name[20])                                  * 
 *                                                                                   *
 * Puts the element name of the atom with element number el_nr in el_name            *
 *************************************************************************************/

void find_el_name(int el_nr, char *atom_name)
{
#define NR_ELEMENTS 19 		/* element number of highest defined element */
    const char el_name[NR_ELEMENTS][19] =  {"H",   /* List of atomic names sorted by element (atomic) number */
					    "He", /* Starting by 0 for H */ 
    					    "Li", 
					    "Be", 
					    "B",
					    "C",
					    "N",
					    "O",					   
				/*	    "F",
					    "Ne",
					    "Na",
					    "Mg",
					    "Al",	*/
					    "Si",	
					    "Si_val",
					    "Si4+",
					    "O-",
				            "P", 
					    "In",
					    "In3+",
					    "Au",
					    "Mn",
					    "Fe",
					    "Pb"
};

    strcpy(atom_name, &(el_name[el_nr - 1][0]));

    return;
}

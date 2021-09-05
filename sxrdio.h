/****************************************************************
 * Contains I/O functions used by sxrdcalc.c			*
 ****************************************************************/

/************************************************************************************************
 * unsigned int read_struc(struct atom_structure *read_atom, struct latt_par_structure *lattice_par, char filename[])							       	*
 * 												*
 * Reads surface structure file into lattice_par and atoms					*
 * Returns number of atoms read.							      	*
 *											       	*
 * filename: Name of structure file							       	*
 * 											       	*
 * Format of structure file:								       	*
 * 1st line: Any comment. is ignored							       	*
 * 2nd line: lattice paramters and angles: a b c alpha beta gamma			       	*
 * other lines: Atomic positions: Atom name x y z 					       	*
 * Atom name is two characters abbrevation, x,y,z in fractional coordinates.		       	* 
 ************************************************************************************************/

extern unsigned int read_struc(struct atom_structure *read_atom, struct latt_par_structure *lattice_par, char filename[]);



/************************************************************************************************
 * unsigned int read_exp_struc_f(char filename[])							*
 * 												*
 * Reads file with experimental structure factors into struc_f_exp				*
 * Returns number of structure factors read							*
 *												*
 * filename: Name of file with experimental structure factors 					*
 * 												*
 * Format of file:									 	*
 * 1st line: Any comment. is ignored								*
 * other lines: h   k   l   structure factor   sigma						*
 *												*
 * sigma is the error of the measured  structure factors					*
 ************************************************************************************************/

extern unsigned int read_exp_struc_f(char filename[]);



/************************************************************************************************
 * int write_structure_factor(char filename[], size_t nr_f, char write_comment[])						       	*
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

extern void write_structure_factor(char filename[], size_t nr_f, char write_comment[]);


/*************************************************************************************
 * void usage(char *program_name)                                                                  * 
 *                                                                                   *
 * Print information how to use the program and exit                                 *
 *************************************************************************************/

extern void usage(char name[]);


/********************************************************************************************************
 * int read_fit_file(struct atom_structure *read_atom, struct latt_par_structure *lattice_par, struct displacements *read_displ, char filename[]) *
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
 * Fixed parameters: Parameters can be fixed by the keyword "fix" and the Nr of the parameter.		*
 *													*
 * 													*
 ********************************************************************************************************/
int read_fit_file(struct atom_structure *read_atom, struct latt_par_structure *lattice_par, struct displacements *read_displ, int *fixed_parameters, double *start_par,  char filename[]);


/****************************************************************************************************************
 * write_fitted_coordinates (struct atom_structure positions, size_t nr_atoms, 					*
 * struct latt_par_structure lattice_par, char *filename[])					*
 *	 													*
 * Write the fitted coordinates to file given by filename							*
 *														*
 * positions: fitted positions											*
 * nr_atoms: Number of atoms											*
 * filename: Filenames for different domains									*
 ****************************************************************************************************************/

void write_fitted_coordinates (struct atom_structure *positions, size_t nr_atoms, struct latt_par_structure latt_par, char write_comment[], char *filename);

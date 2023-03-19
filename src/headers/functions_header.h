
#include "../functions/functions.c"

/*	
	Header to include in main file 
	contains function declarations


*/

// Initialize domain length, bandwidth and the smallest dimensions of a structure in the domain
void initialize_parameters(double *z_max, double *f_max, double *d_min);

// Initialize grid parameters and source function
void initialize_grid(double z_max, double f_max, double d_min, int *Z_steps, int *T_steps, double *dt, double *dz);

// Initialize fields to 0
void initialize_fields(double *E, double *H, int n);

// Create source function
void source_func(int n, double g[n], double dt, double f_max);

// Create dielectric
void dielectric(double *eps_r, double *mu_r, int n, double dz);

// Apply dielectric to update factors
void apply_dielectric(double *f_e, double *f_h, double *eps_r, double *mu_r, double dt, double dz, int n);

// Write array to file
void write_to_file(FILE *file_ptr, char *PATH ,int n, double arr[n]);

// Update fields(one time step)
void update_fields(double *E, double *H, int nz, double *f_h, double *f_e, double ei[3], double hi[3]);

// Insert source to Efield
void insert_source(double* E, double g, int  k, bool soft_source);

// Creates interval with n steps and step ds
void interval(int n, double ds, double *s);
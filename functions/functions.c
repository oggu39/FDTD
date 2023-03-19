
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
/*
	defines functions used in the main scripts
*/
#define c0 299792458.0
#define mu0 4*M_PI*pow(10,-7)
#define eps0 1/(c0*c0*mu0)
// Declarations of misc functions:
// 1.) min function 
double min(double a,double b);

// ---------------- Functions ---------------------------


// 1.) Initializing the parameters
void initialize_parameters(double *z_max, double *f_max, double *d_min){
	printf("###################################\n");
	printf(" Finite difference time domain in 1 dimension \n");
	printf("###################################\n");
	printf("\n\n");
	
	
	printf("Enter the size of the domain[m]: ");
	scanf("%lf", z_max);
	printf("\n");
	
	printf("Enter the bandwidth[MHz]: ");
	scanf("%lf", f_max);
	printf("\n");
	
	printf("Enter the size of the smallest structure in the domain[m]: ");
	scanf("%lf", d_min);

	
	*f_max = pow(10,6)*(*f_max); // Convert to Hz
	
	printf("\nSize of domain: %.3lf m\n", *z_max);
	printf("Maximum frequency: %.3lf MHz \n", (*f_max)/pow(10,6));
	printf("Size of smallest structure in the domain: %.3lf m \n\n", *d_min);	
}


// 2.) Initialize grid and calculate dt, dz, Z_steps, T_steps, source function g_e
void initialize_grid(double z_max, double f_max, double d_min, int *Z_steps, int *T_steps, double *dt, double *dz){
	
	double c = 299792458.0;
	double lambda_min = c/f_max;
	
	int N_d = 10; // We want to resolve the smallest size d_min by approx 5 cells
	float N_lambda = 20; // We want to resolve the smallest wavelength by approx 10 cells
	
	double dz_f = min((double)d_min/N_d, (double)lambda_min/N_lambda); // Space-step dz
	int nz = (int)ceil(z_max/dz_f); // Nr of cells in the domain
	*dz = (double)z_max/nz; // Snap the domain into a integer number of cells
	
	
	double tau = 0.5/f_max; // Width of the gaussian pulse[s]
	double n_min = 1.0; // Smallest refractive index;
	double n_max = 2.0; // Largest refractive indez in the domain. Should be entered by user.
	
	*dt = ((*dz)*n_min)/(2*c); // appropriate time step that satisfies CFL automatically
	double t0 = 6*tau; // Time delay untill the peak of the gaussian pulse ~6 taus
	double t_prop = (2*z_max*n_max)/(c); // maximum propagation time in the domain(overestimate);
	double T_final = 6*tau+12*t_prop;
	
	int nt = (int)ceil(T_final/(*dt)); // number of time steps
	
	// Print parameters ----------------
	printf("Number of grid cells: %d \t Number of time steps: %d\n", nz, nt);
	printf("Grid resolution: %lf [m]  Time step size: %lf e-9 [s] \n",(*dz),pow(10,9)*(*dt));
	printf("\n\t-------- Source --------- \n");
	printf("Gaussian: Amplitude = %.2f[V/m]\t offset = %lfe-9 [s]\t width = %lfe-9 [s]\n\n", 1.00, pow(10,9)*t0, pow(10,9)*tau);

		
	
	// Values to return:
	*T_steps = nt;
	*Z_steps = nz;

	
}


// 3.) Create source function(gaussian with amplitude 1)
void source_func(int n, double g[n], double dt, double f_max){
	double t[n];
	double a;
	double tau = 0.5/f_max;
	double t0 = 6*tau;
	
	for(int i = 0; i < n; i++){
		t[i] = i*dt;
		a = (t[i]-t0)/tau;
		g[i] = 1.0*exp(-pow(a,2));
	}
	
}

// 4.) Initialize fields
void initialize_fields(double *E, double *H, int n){
	for(int i = 0; i < n; i++){
		E[i] = 0;
		H[i] = 0;
	}
}

// 5.) Update fields
void update_fields(double *E, double *H, int nz, double *f_h, double *f_e, double ei[3], double hi[3]){
	// f_h update factor for H
	// f_e update factor for E
	// nz number of cells in the z direction
	int k;
	
	// Update H field
	for(k = 0; k < nz-1; k++){
		H[k] = H[k] + f_h[k]*(E[k+1]-E[k]);
	}
	H[nz-1] = H[nz-1] + f_h[nz-1]*(ei[2]-E[nz-1]); // right boundary
	hi[2] = hi[1]; hi[1] = hi[0]; hi[0] = H[0]; // Absorbing boundary in one dimension
	
	
	
	// Update E field
	E[0] = E[0] + f_e[0]*(H[0]-hi[2]); // Left boundary
	for(k = 1; k < nz;  k++){
		E[k] = E[k] + f_e[k]*(H[k]-H[k-1]);
	}
	ei[2] = ei[1]; ei[1] = ei[0]; ei[0] = E[nz-1]; // Absorbing layer 
}


// 6. insert source
void insert_source(double* E, double g, int  k, bool soft_source){
	// Inserts source g on to the E field at location k*dz;
	if(soft_source){
		E[k] = E[k] + g;
	}else{
		E[k] = g; // Hard source 
	}
	
	
}

// 7.) Place simple dielectric in domain
void dielectric(double *eps_r, double *mu_r, int n, double dz){
	double relative_permittivity, relative_permeability, a, b;
	printf("\n------------ Dielectric ------------\n");
	printf("Permittivity: ");
	scanf("%lf", &relative_permittivity);
	printf("Permeability: ");
	scanf("%lf", &relative_permeability);
	printf("Dielectric left bound[m]: ");
	scanf("%lf", &a);
	printf("Dielectric right bound[m]: ");
	scanf("%lf", &b);
	

	
	int k_a, k_b;
	
	if((a > n*dz) || (a < 0) || (b <= a) || (b > n*dz)){
		printf("Error: Dielectric outside the domain\n");
		exit(1);
	}else{
		k_a = (int)floor(a/dz);
		k_b = (int)ceil(b/dz);
		
		for(int k = 0; k < n; k++){
				if((k>=k_a) && (k <= k_b)){
					eps_r[k] = relative_permittivity;
					mu_r[k] = relative_permeability;
				}else{
					eps_r[k] = 1;
					mu_r[k] = 1;
				}
		}
	}
	
	printf("\n------------ Dielectric parameters ------------\n");
	printf("Location: %.2lf [m]\n", a);
	printf("Dielectric width: %.2lf\n", b-a);
	printf("Relative permittivity: %.2lf\n", relative_permittivity);
	printf("Relative permeability: %.2lf\n", relative_permeability);
}

// Apply dielectric
void apply_dielectric(double *f_e, double *f_h, double *eps_r, double *mu_r, double dt, double dz, int n){
	double C = -dt/(sqrt(mu0*eps0)*dz);
	for(int k = 0; k < n; k++){
		f_e[k] = (1/eps_r[k])*C;
		f_h[k] = (1/mu_r[k])*C;
	}
}



// ----------------------- Misc functions ----------------------------------- 


// n.) Write data in array to file in path
void write_to_file(FILE *file_ptr, char *PATH ,int n, double arr[n]){
	if(file_ptr == NULL){
		printf("File NULL");
		exit(1);
	} else{
		for(int i = 0; i < n; i++){	
			fprintf(file_ptr, "%.15lf \t", arr[i]);
		}
		fprintf(file_ptr, "\n");
		
	}
}

// minimum of two values
double min(double a,double b){
	if(a < b){
		return a;
	} else{
		return b;
	}
}

// Creates an interval with n steps and step ds
void interval(int n, double ds, double *s){
	for(int i = 0; i < n; i++){
		s[i] = i*ds;
	}
}
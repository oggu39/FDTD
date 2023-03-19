#include <stdio.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <stdlib.h>
#include "headers/functions_header.h"
#include <time.h>
#include <unistd.h>


#define c0 299792458.0
#define mu0 4*M_PI*pow(10,-7)
#define eps0 1/(c0*c0*mu0)

int main(){

	double z_max, f_max, d_min; // Problem parameters
	initialize_parameters(&z_max, &f_max, &d_min); // Problem parameters set by user
	
	double dz,dt;
	int Z_steps, T_steps;

	initialize_grid(z_max, f_max, d_min, &Z_steps, &T_steps, &dt, &dz); // Calculate simulation parameters
	
	// Calculate update factors
	double f_e[Z_steps];
	double f_h[Z_steps]; 

	
	// Declare fields, source function and time variable
	double Ex[Z_steps];
	double Hy[Z_steps];
	double g[T_steps]; 
	double t = 0.0;
	double E_reflect[Z_steps];
	double E_transmit[Z_steps];
	double eps_r[Z_steps];
	double mu_r[Z_steps];
	
	// Initialize fields to 0 and create gaussian pulse source
	initialize_fields(Ex,Hy, Z_steps);
	source_func(T_steps, g, dt, f_max); // Creates gaussian source

	// Apply dielectric
	dielectric(eps_r, mu_r, Z_steps, dz);
	apply_dielectric(f_e, f_h, eps_r, mu_r, dt, dz, Z_steps);
	
	
	
	
	// Path and pointers to files for exporting the simulation
	//char *PATH_Efield = "G:\\c\\FDTD\\exported_data\\E_data.txt";
	//char *PATH_Hfield = "G:\\c\\FDTD\\exported_data\\H_data.txt";
	char *PATH_time = "G:\\c\\FDTD\\exported_data\\t_data.txt";
	//char *PATH_z = "G:\\c\\FDTD\\exported_data\\z_data.txt";	
	char *PATH_reflection = "G:\\c\\FDTD\\exported_data\\reflected_field.txt";
	char *PATH_transmission = "G:\\c\\FDTD\\exported_data\\transmitted_field.txt";
	char *PATH_source = "G:\\c\\FDTD\\exported_data\\source_field.txt";
	
	//FILE *E_data = fopen(PATH_Efield, "w");
	FILE *t_data = fopen(PATH_time, "w");
	//FILE *H_data = fopen(PATH_Hfield, "w");
	//FILE *z_data = fopen(PATH_z, "w");
	FILE *reflection_data = fopen(PATH_reflection, "w");
	FILE *transmission_data = fopen(PATH_transmission, "w");
	FILE *source_data = fopen(PATH_source, "w");
	
	for(int i = 0; i < T_steps; i++){
		write_to_file(source_data, PATH_source, 1, &g[i]);
	}
	
	//double z[Z_steps];
	/*
	interval(Z_steps, dz, z);
	for(int j = 0; j < Z_steps; j++){
		write_to_file(z_data, PATH_z, 1, &z[j]);
	}
	fclose(z_data);
	*/
	
	
	int src_location = (int)floor(Z_steps/4);               
	double ei[3] = {0,0,0};
	double hi[3] = {0,0,0};
	// ------------------- main time loop ----------------------------------
	for(int n = 0; n < T_steps; n++){
		
		update_fields(Ex, Hy, Z_steps, f_h, f_e, ei, hi);
		insert_source(Ex, g[n], src_location, true);
		t = n*dt;
		E_reflect[n] = Ex[10]; // Record The field value 10 cells from the left edge
		E_transmit[n] = Ex[Z_steps-10]; // Record field value 10 cell from the right edge
		
		write_to_file(reflection_data, PATH_reflection, 1, &E_reflect[n]);
		write_to_file(transmission_data, PATH_transmission, 1, &E_transmit[n]);
		write_to_file(t_data, PATH_time, 1, &t);
	/* Saving various data
		write_to_file(E_data, PATH_Efield ,Z_steps, Ex);
		write_to_file(H_data, PATH_Hfield, Z_steps, Hy);
		
	*/
	
	}
	
	//fclose(E_data);
	//fclose(H_data);
	fclose(t_data);
	fclose(reflection_data);
	fclose(transmission_data);
	fclose(source_data);
	//fclose(z_data);
	
	/*
	double g_e[T_steps];
	source_func(T_steps, g_e, dt, f_max); // Source function g(t);
	*/

	printf("\t . \n");
	printf("\t . \n");
	printf("\t . \n");
	int nr_seconds = 1;
	usleep(1000000*nr_seconds);
	printf("----------------------------------------------\n");
	printf("--------------- Simulation done --------------\n");
	printf("----------------------------------------------\n");
	system("python G:\\C\\FDTD\\exported_data\\analyse.py");
	return 0;	
}
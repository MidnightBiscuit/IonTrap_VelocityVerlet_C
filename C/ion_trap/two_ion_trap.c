// gcc two_ion_trap.c -o two_ion_trap -lm

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

int main(){
	// maths and phys
	const float pi    = 3.14159265359;
	const float ke    = 8.9875517873681764e9;
	
	// atoms
	const int N_ions  = 2;
	const float m_Ca[2]  = {40*1.66053886e-27,42*1.66053886e-27};
	const float Q_Ca  = 1.60217646e-19;
	
	// trap
	const float U_dc  = 1000;
	const float z0_2  = 0.02*0.02; // m
	const float w_z_2[2] = {2*Q_Ca*U_dc/m_Ca[0]/z0_2, 2*Q_Ca*U_dc/m_Ca[1]/z0_2};
	
	// integration
	const float dt = 5e-9;
	const float dt1 = 0.5*5e-9;
	const float dt2 = 0.5*5e-9 * 0.5*5e-9;
	//~ const int i_free_fly = round(65/sqrtf(w_z_2[0])/dt);
	const int i_free_fly = 4000;
	
	float rjk;
	float r2inv;
	const float softening = 1e-20;
	printf("constants done\n");
	
	float* t_act = (float*)malloc(i_free_fly*sizeof(float));
	float* r_z   = (float*)malloc(N_ions*i_free_fly*sizeof(float));
	float* v_z   = (float*)malloc(N_ions*i_free_fly*sizeof(float));
	float* a_z   = (float*)malloc(N_ions*i_free_fly*sizeof(float));
	printf("%d\n",N_ions*i_free_fly);
	// initialization
	t_act[0] = 0;
	r_z[0]   = 1e-3;
	r_z[0+i_free_fly]   = -1e-3;
	for (int i=0;i<N_ions*i_free_fly+1;i=i+i_free_fly){
	v_z[i]   = 0;
	a_z[i]   = 0;
	}
	printf("initialization done\n");
	// Velocity-Verlet algorithm
	for (int i=1;i<i_free_fly;++i){ // i the time step
		t_act[i] = dt*i;
		for (int j=i;j<N_ions*i_free_fly+1;j=j+i_free_fly){ // j the ion under calculation
		// update position
		r_z[j] = r_z[j-1] + v_z[j-1]*dt + 0.5*a_z[j-1]*dt*dt;
		// update acceleration
		for (int k=i;k<N_ions*i_free_fly;k=k+i_free_fly){ // k the neighbourhood ion
			printf('%k\n',k)
			rjk = r_z[j]-r_z[k];
			r2inv = 1/sqrt(rjk*rjk + softening);
			r2inv = ke*Q_Ca*Q_Ca/m_Ca[0] * r2inv*r2inv*r2inv;
			a_z[j] = -w_z_2[0]*r_z[j] + rjk*r2inv;
		}
		// update velocity
		v_z[j] = v_z[j-1] + dt1 * (a_z[j-1] + a_z[j]);
		}
	//~ printf("%i\n",i);
	}
	
	FILE *fp;

    fp = fopen("two_ions_xva.bin", "w+");
    fwrite(t_act,sizeof(float),i_free_fly,fp);
    fwrite(r_z,sizeof(float),N_ions*i_free_fly,fp);
    fwrite(v_z,sizeof(float),N_ions*i_free_fly,fp);
    fwrite(a_z,sizeof(float),N_ions*i_free_fly,fp);
    fclose(fp);
    
	//~ for( int i = 0; i < 50; i++ ) {
			//~ printf("t : %e r_z: %e a_z: %e v_z: %e\n", t_act[i], r_z[i], a_z[i], v_z[i]);
		//~ }

	
	//~ printf("w_z = %f\n",sqrt(w_z_2));
	//~ printf("%f\n",U_dc);
	//~ printf("sqrtf(w_z_2) = %f\n",1/sqrtf(w_z_2));
	//~ printf("i_free_fly = %i\n",i_free_fly);
	//~ printf("%li %i\n",sizeof(float),i_free_fly);
	//~ printf("w_z = %f\n",sqrt(w_z_2[0])/2/pi);
	//~ printf("%f\n",57590*dt);
	//~ printf("%li\n",sizeof(r_z)/sizeof(*r_z));
	sleep(3);
	return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

int main(){
	const float pi    = 3.14159265359;
	const float m_Ca  = 40*1.66053886e-27;
	const float Q_Ca  = 1.60217646e-19;
	
	const float U_dc  = 1000;
	const float z0_2  = 0.02*0.02; // m
	const float w_z_2 = 2*Q_Ca*U_dc/m_Ca/z0_2;
	
	const float dt = 5e-9;
	const float dt1 = 0.5*5e-9;
	const float dt2 = 0.5*5e-9 * 0.5*5e-9;
	const int i_free_fly = round(65/sqrtf(w_z_2)/dt);
	
	float* t_act = (float*)malloc(i_free_fly*sizeof(float));
	float* r_z   = (float*)malloc(i_free_fly*sizeof(float));
	float* v_z   = (float*)malloc(i_free_fly*sizeof(float));
	float* a_z   = (float*)malloc(i_free_fly*sizeof(float));
	t_act[0] = 0;
	r_z[0]   = 1e-3;
	v_z[0]   = 0;
	a_z[0]   = 0;
	
	// Velocity-Verlet algorithm
	for (int i=1;i<i_free_fly;++i){
		t_act[i] = dt*i;
		// update position
		r_z[i] = r_z[i-1] + v_z[i-1]*dt + 0.5*a_z[i-1]*dt*dt;
		// update acceleration
		a_z[i] = -w_z_2*r_z[i];
		// update velocity
		v_z[i] = v_z[i-1] + dt1 * (a_z[i-1] + a_z[i]);
	}
	
	FILE *fp;

    fp = fopen("xva.bin", "w+");
    fwrite(t_act,sizeof(float),i_free_fly,fp);
    fwrite(r_z,sizeof(float),i_free_fly,fp);
    fwrite(v_z,sizeof(float),i_free_fly,fp);
    fwrite(a_z,sizeof(float),i_free_fly,fp);
    fclose(fp);
    
	for( int i = 0; i < 50; i++ ) {
			printf("t : %e r_z: %e a_z: %e v_z: %e\n", t_act[i], r_z[i], a_z[i], v_z[i]);
		}

	
	//~ printf("w_z = %f\n",sqrt(w_z_2));
	//~ printf("%f\n",U_dc);
	//~ printf("sqrtf(w_z_2) = %f\n",1/sqrtf(w_z_2));
	printf("i_free_fly = %i\n",i_free_fly);
	printf("%li %i\n",sizeof(float),i_free_fly);
	printf("w_z = %f\n",sqrt(w_z_2)/2/pi);
	//~ printf("%f\n",57590*dt);
	//~ printf("%li\n",sizeof(r_z)/sizeof(*r_z));
	sleep(3);
	return 0;
}

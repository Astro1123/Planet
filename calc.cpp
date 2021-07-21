#include "const.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int pn = 0;
int mag = 100;


double radius(double x,double y,double z) {
	return sqrt(x*x+y*y+z*z);
}

double COG(double *x) {
	double p1=0.0;
	double p2=0.0;
	int i;
	for(i=0;i<pn;i++) {
		p1+=data[i].m;
		p2+=data[i].m*x[i];
	}
	return p2/p1;
}

void center(int type) {
	int i;
	double x[P_MAX],r[6];
	if(type >= 0 && type < pn) {
		for(i=0;i<3;i++) {
			r[  i] = data[type].x[i];
			r[i+3] = data[type].v[i];
		}
	} else if(type == -2) {
		return;
	} else {
		for(i=0;i<pn;i++) {
			x[i] = data[i].x[0];
		}
		r[0] = COG(x);
		for(i=0;i<pn;i++) {
			x[i] = data[i].x[1];
		}
		r[1] = COG(x);
		for(i=0;i<pn;i++) {
			x[i] = data[i].x[2];
		}
		r[2] = COG(x);
		for(i=0;i<pn;i++) {
			x[i] = data[i].v[0];
		}
		r[3] = COG(x);
		for(i=0;i<pn;i++) {
			x[i] = data[i].v[1];
		}
		r[4] = COG(x);
		for(i=0;i<pn;i++) {
			x[i] = data[i].v[2];
		}
		r[5] = COG(x);
	}
	for(i=0;i<pn;i++) {
		data[i].x[0] -= r[0]; 
		data[i].x[1] -= r[1];
		data[i].x[2] -= r[2];
		data[i].v[0] -= r[3];
		data[i].v[1] -= r[4];
		data[i].v[2] -= r[5];
	}
	return;
}

int func(int n,double *f) {
	double r;
	int i;
	f[0] = 0.0;
	f[1] = 0.0;
	f[2] = 0.0;
	if(pn <= 1)
		return -1;
	if(n >= pn)
		return -1;
	if(n < 0)
		return -1;
	
	for(i=0;i<pn;i++) {
		if( i == n )
			continue;
		r     = sqrt( pow(data[i].x[0] - data[n].x[0],2) + pow(data[i].x[1] - data[n].x[1],2) + pow(data[i].x[2] - data[n].x[2],2) );
		f[0] += CONST_G * ( data[i].m * data[n].m / pow(r,3) ) * (data[i].x[0] - data[n].x[0]);
		f[1] += CONST_G * ( data[i].m * data[n].m / pow(r,3) ) * (data[i].x[1] - data[n].x[1]);
		f[2] += CONST_G * ( data[i].m * data[n].m / pow(r,3) ) * (data[i].x[2] - data[n].x[2]);
	}
	
	return 0;
}


int VelocityVerlet(int n,double dt,double *r) {
	double f[3],f2[3];
	int i;
	
	func(n,f);	
	for(i=0;i<3;i++) {
		r[i] = 0.5 * dt*dt * f[i] / data[n].m + dt*r[3+i] + data[n].x[i];
	}
	
	func(n,f2);
	for(i=0;i<3;i++) {
		r[3+i] = 0.5 * dt * (f[i] + f2[i]) / data[n].m + data[n].v[i];
	}
	
	return 0;
}

int Leapfrog(int n,double dt,double *r) {
	double f[3],vb[3];
	int i;
	
	func(n,f);
	if(data[n].vhf==-1){	
		for(i=0;i<3;i++) {
			data[n].vh[i] = data[n].v[i] - 0.5 * f[i] / data[n].m * dt;
		}
		data[n].vhf=0;
	}
	
	for(i=0;i<3;i++) {
		vb[i]=data[n].vh[i];
		data[n].vh[i]+=f[i] / data[n].m * dt;
		r[i] = dt * data[n].vh[i] + data[n].x[i];
		r[3+i] = 0.5 * (data[n].vh[i] + vb[i]);
	}
	
	return 0;
}

int EulerRichardson(int n,double dt,double *r) {
	double f[3],vm[3],xm[3];
	int i;
	
	func(n,f);	
	for(i=0;i<3;i++) {
		vm[i] = 0.5 * dt * f[i] / data[n].m + data[n].v[i];
		xm[i] = r[i];
		r[i] = 0.5 * dt * vm[i] + data[n].x[i];
	}
	
	func(n,f);
	for(i=0;i<3;i++) {
		r[3+i] = dt * f[i] / data[n].m + data[n].v[i];
		r[i] = xm[i] + vm[i] * dt;
	}
	
	return 0;
}

int RungeKutta(int n,double dt,double *r) {
	double f[3],k1[3],k2[3],k3[3],k4[3],kx1[3],kx2[3],kx3[3],kx4[3];
	int i;
	
	func(n,f);
	for(i=0;i<3;i++) {
		kx1[i] = dt * data[n].v[i];
		k1[i]  = dt * f[i] / data[n].m;
		data[n].x[i] += 0.5 * kx1[i];
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx2[i] = dt * (data[n].v[i] + k1[i]/2);
		k2[i] = dt * f[i] / data[n].m;
		data[n].x[i] -= 0.5 * kx1[i];
		data[n].x[i] += 0.5 * kx2[i];
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx3[i] = dt * (data[n].v[i] + k2[i]/2);
		k3[i] = dt * f[i] / data[n].m;
		data[n].x[i] -= 0.5 * kx2[i];
		data[n].x[i] += kx3[i];
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx4[i] = dt * (data[n].v[i] + k3[i]);
		k4[i] = dt * f[i] / data[n].m;
		data[n].x[i] -= kx3[i];
		r[3+i] = ( k1[i] + k2[i] * 2.0 + k3[i] * 2.0 + k4[i] ) / 6.0 + data[n].v[i];
		r[i] = ( kx1[i] + kx2[i] * 2.0 + kx3[i] * 2.0 + kx4[i] ) / 6.0 + data[n].x[i];
	}
	
	return 0;
}

int RungeKutta4_38(int n,double dt,double *r) {
	double f[3],k1[3],k2[3],k3[3],k4[3],kx1[3],kx2[3],kx3[3],kx4[3];
	int i;
	
	func(n,f);
	for(i=0;i<3;i++) {
		kx1[i] = dt * data[n].v[i];
		k1[i]  = dt * f[i] / data[n].m;
		data[n].x[i] += kx1[i] / 3;
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx2[i] = dt * (data[n].v[i] + k1[i]/3);
		k2[i] = dt * f[i] / data[n].m;
		data[n].x[i] -= kx1[i] / 3;
		data[n].x[i] += -kx1[i] *2/3 + kx2[i];
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx3[i] = dt * (data[n].v[i] + k2[i] -k1[i]*2/3);
		k3[i] = dt * f[i] / data[n].m;
		data[n].x[i] -= -kx1[i] *2/3 + kx2[i];
		data[n].x[i] += kx1[i]-kx2[i]+kx3[i];
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx4[i] = dt * (data[n].v[i] + k1[i]-k2[i]+k3[i]);
		k4[i] = dt * f[i] / data[n].m;
		data[n].x[i] -= kx1[i]-kx2[i]+kx3[i];
		r[3+i] = ( k1[i] + k2[i] * 3.0 + k3[i] * 3.0 + k4[i] ) / 8.0 + data[n].v[i];
		r[i] = ( kx1[i] + kx2[i] * 3.0 + kx3[i] * 3.0 + kx4[i] ) / 8.0 + data[n].x[i];
	}
	
	return 0;
}

int Euler(int n,double dt,double *r) {
	double f[3];
	int i;
	
	func(n,f);
	
	for(i=0;i<3;i++) {
		r[3+i] = data[n].v[i] + f[i] * dt / data[i].m;
		r[i] = dt*r[3+i] + data[n].x[i];
	}
	
	return 0;
}

int Ralston(int n,double dt,double *r) {
	double f[3];
	int i,k1[3],k2[3],kx1[3],kx2[3];
	
	func(n,f);
	for(i=0;i<3;i++) {
		kx1[i] = dt * data[n].v[i];
		k1[i]  = dt * f[i] / data[n].m;
		data[n].x[i] += kx1[i] *2/3;
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx2[i] = dt * (data[n].v[i] + k1[i] *2/3);
		k2[i]  = dt * f[i] / data[n].m;
		data[n].x[i] -= kx1[i] *2/3;
		r[3+i] = data[n].v[i] + (k1[i]+k2[i]*3)/4;
		r[i] = data[n].x[i] + (kx1[i]+kx2[i]*3)/4;
	}
	
	return 0;
}

int impEuler(int n,double dt,double *r) {
	double f[3];
	int i,k1[3],k2[3],kx1[3],kx2[3];
	
	func(n,f);
	for(i=0;i<3;i++) {
		kx1[i] = dt * data[n].v[i];
		k1[i]  = dt * f[i] / data[n].m;
		data[n].x[i] += kx1[i];
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx2[i] = dt * (data[n].v[i] + k1[i]);
		k2[i]  = dt * f[i] / data[n].m;
		data[n].x[i] -= kx1[i];
		r[3+i] = data[n].v[i] + (k1[i]+k2[i])/2;
		r[i] = data[n].x[i] + (kx1[i]+kx2[i])/2;
	}
	
	return 0;
}

int MidPoint(int n,double dt,double *r) {
	double f[3];
	int i,k1[3],k2[3],kx1[3],kx2[3];
	
	func(n,f);
	for(i=0;i<3;i++) {
		kx1[i] = dt * data[n].v[i];
		k1[i]  = dt * f[i] / data[n].m;
		data[n].x[i] += kx1[i]/2;
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx2[i] = dt * (data[n].v[i] + k1[i]);
		k2[i]  = dt * f[i] / data[n].m;
		data[n].x[i] -= kx1[i]/2;
		r[3+i] = data[n].v[i] + (k1[i]+k2[i])/2;
		r[i] = data[n].x[i] + (kx1[i]+kx2[i])/2;
	}
	
	return 0;
}

int Kutta3o(int n,double dt,double *r) {
	double f[3];
	int i,k1[3],k2[3],k3[3],kx1[3],kx2[3],kx3[3];
	
	func(n,f);
	for(i=0;i<3;i++) {
		kx1[i] = dt * data[n].v[i];
		k1[i]  = dt * f[i] / data[n].m;
		data[n].x[i] += kx1[i] *0.5;
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx2[i] = dt * (data[n].v[i] + k1[i] *0.5);
		k2[i]  = dt * f[i] / data[n].m;
		data[n].x[i] -= kx1[i] *0.5;
		data[n].x[i] += kx2[i]*2-kx1[i];
	}
	for(i=0;i<3;i++) {
		kx3[i] = dt * (data[n].v[i] + k2[i]*2 - k1[i]);
		k3[i]  = dt * f[i] / data[n].m;
		data[n].x[i] -= kx2[i]*2-kx1[i];
		r[3+i] = data[n].v[i] + (k1[i]+k2[i]*4+k3[i])/6;
		r[i] = data[n].x[i] + (kx1[i]+kx2[i]*4+kx3[i])/6;
	}
	
	return 0;
}

int _RungeKuttaFehlberg(int n,double dt,double eps,double *r,double *ar) {
	double f[3],k1[3],k2[3],k3[3],k4[3],k5[3],k6[3],kx1[3],kx2[3],kx3[3],kx4[3],kx5[3],kx6[3];
	double s[6],e[6],em=0.0,ar2[6];
	int i;
	
	func(n,f);
	for(i=0;i<3;i++) {
		kx1[i] = dt * data[n].v[i];
		k1[i]  = dt * f[i] / data[n].m;
		data[n].x[i] += kx1[i] / 4;
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx2[i] = dt * (data[n].v[i] + k1[i]/4);
		k2[i] = dt * f[i] / data[n].m;
		data[n].x[i] -= kx1[i] / 4;
		data[n].x[i] += kx1[i] / 32 * 3 + kx2[i] / 32 * 9;
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx3[i] = dt * (data[n].v[i] + k2[i]/32*9 + k1[i]/32*3);
		k3[i] = dt * f[i] / data[n].m;
		data[n].x[i] -= kx1[i] / 32 * 3 + kx2[i] / 32 * 9;
		data[n].x[i] += kx1[i] / 2197 * 1932 - kx2[i] / 2197 * 7200 + kx3[i] / 2197 * 7296;
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx4[i] = dt * (data[n].v[i] + k1[i] / 2197 * 1932 - k2[i] / 2197 * 7200 + k3[i] / 2197 * 7296);
		k4[i] = dt * f[i] / data[n].m;
		data[n].x[i] -= kx1[i] / 2197 * 1932 - kx2[i] / 2197 * 7200 + kx3[i] / 2197 * 7296;
		data[n].x[i] += kx1[i] / 216 * 439 - kx2[i] * 8 + kx3[i] / 513 * 3680 - kx4[i] / 4104 * 845;
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx5[i] = dt * (data[n].v[i] + k1[i] / 216 * 439 - k2[i] * 8 + k3[i] / 513 * 3680 - k4[i] / 4104 * 845);
		k5[i] = dt * f[i] / data[n].m;
		data[n].x[i] -= kx1[i] / 216 * 439 - kx2[i] * 8 + kx3[i] / 513 * 3680 - kx4[i] / 4104 * 845;
		data[n].x[i] += -kx1[i] / 27 * 8 + kx2[i] * 2 - kx3[i] / 2565 * 3544 + kx4[i] / 4104 * 1859 - kx5[i] / 40 * 11;
	}
	func(n,f);
	for(i=0;i<3;i++) {
		kx6[i] = dt * (data[n].v[i] + (-k1[i]) / 27 * 8 + k2[i] * 2 - k3[i] / 2565 * 3544 + k4[i] / 4104 * 1859 - k5[i] / 40 * 11);
		k6[i] = dt * f[i] / data[n].m;
		data[n].x[i] -= -kx1[i] / 27 * 8 + kx2[i] * 2 - kx3[i] / 2565 * 3544 + kx4[i] / 4104 * 1859 - kx5[i] / 40 * 11;
		ar[i+3]=( k1[i] / 216 * 25 + k3[i] * 1408 / 2565 + k4[i] * 2197 / 4104 - k5[i] / 5 );
		ar[i]  =( kx1[i] / 216 * 25 + kx3[i] * 1408 / 2565 + kx4[i] * 2197 / 4104 - kx5[i] / 5 );
		r[3+i] = ar[3+i] + data[n].v[i];
		r[i] = ar[i] + data[n].x[i];
		s[3+i] = ( k1[i] / 135 * 16 + k3[i] * 6656 / 12825 + k4[i] * 28561 / 56430 - k5[i] / 50 * 9 + k6[i] / 55 * 2 ) + data[n].v[i];
		s[i] = ( kx1[i] / 135 * 16 + kx3[i] * 6656 / 12825 + kx4[i] * 28561 / 56430 - kx5[i] / 50 * 9 + kx6[i] / 55 * 2 ) + data[n].x[i];
		e[3+i] = s[3+i]-r[3+i];
		e[i] = s[i]-r[i];
		if(em<e[3+i])
			em=e[3+i];
		if(em<e[i])
			em=e[i];
	}
	if (abs(em)/dt>eps) {
		_RungeKuttaFehlberg(n,dt/2,eps,r,ar);
		for(i=0;i<3;i++){
			data[n].v[i]+=ar[i+3];
			data[n].x[i]+=ar[i];
		}
		_RungeKuttaFehlberg(n,dt/2,eps,r,ar2);
		for(i=0;i<3;i++){
			data[n].v[i]-=ar[i+3];
			data[n].x[i]-=ar[i];
			ar[i]+=ar2[i];
			ar[i+3]+=ar2[i+3];
		}
	}
	return 0;
}

int RungeKuttaFehlberg(int n,double dt,double eps,double *r) {
	double rc[6],ar[6];
	int i;
	_RungeKuttaFehlberg(n,dt,eps,r,ar);
	for(i=0;i<3;i++){
		r[i+3]=data[n].v[i]+ar[i+3];
		r[i]=data[n].x[i]+ar[i];
	}
	return 0;
}

int PredictorCorrector(int n,double dt,double *r) {
	static int count[P_MAX];
	double f[3],a[P_MAX][4][3],v[P_MAX][4][3],yv,yx;
	int i;
	if(count[n]<4) {
		//RungeKutta(n,dt,*r);
		func(n,f);
		for(i=0;i<3;i++) {
			a[n][count[n]][i] = f[i] / data[n].m;
			v[n][count[n]][i] = r[3+i] - data[n].v[i];
		}
		count[n]++;
	} else {
		for(i=0;i<3;i++) {
			yv = data[n].v[i] + (dt/24)*(55*a[n][3][i]-59*a[n][2][i]+37*a[n][1][i]-9*a[n][0][i]);
			yx = data[n].v[i] + (dt/24)*(55*v[n][3][i]-59*v[n][2][i]+37*v[n][1][i]-9*v[n][0][i]);
		}
	}
	return 0;
}
int modMidPoint(int n,double dt,double *r) {
	return 0;
}

double length(double x) {
	return x / UNIT_AU * mag;
}

double time(double t) {
	return t / UNIT_MIN / UNIT_HOUR / UNIT_DAY;
}

int input(void) {
	char s1[CHAR_MAX]={0};
	char s2[] = ",";
	char *tok;
	int i;
        FILE *fp ;
	
	if((fp=fopen("data.txt","r"))!=NULL){
		while( fgets(s1,CHAR_MAX,fp) != NULL ){
			tok = strtok( s1, s2 );
			data[pn].x[0]=atof(tok);
			tok = strtok( NULL, s2 );
			data[pn].x[1]=atof(tok);
			tok = strtok( NULL, s2 );
			data[pn].x[2]=atof(tok);
			tok = strtok( NULL, s2 );
			data[pn].m=atof(tok);
			tok = strtok( NULL, s2 );
			data[pn].v[0]=atof(tok);
			tok = strtok( NULL, s2 );
			data[pn].v[1]=atof(tok);
			tok = strtok( NULL, s2 );
			data[pn].v[2]=atof(tok);
			tok = strtok( NULL, s2 );
			if(tok == NULL) {
				data[pn].c[0] = 255.0;
				data[pn].c[1] = 255.0;
				data[pn].c[2] = 255.0;
			} else {
				data[pn].c[0]=atof(tok);
				tok = strtok( NULL, s2 );
				data[pn].c[1]=atof(tok);
				tok = strtok( NULL, s2 );
				data[pn].c[2]=atof(tok);
			}
			for(i=0;i<3;i++){
				data[pn].vh[i]=0.0;
			}
			data[pn].vhf=-1;
			if(pn+1>=P_MAX)
				break;
			pn++;
		}
		fclose(fp);
	} else {
		fclose(fp);
		return -1;
	}
	return 0;
}




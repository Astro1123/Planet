#ifndef _CALC_H_
#define _CALC_H_

void center(int);
double radius(double,double,double);
double length(double);
double time(double);
double COG(double);
int input(void);
int func(int,double*);
int Euler(int,double,double*);
int MidPoint(int,double,double*);
int impEuler(int,double,double*);
int VelocityVerlet(int,double,double*);			// 速度ベルレ法（Velocity Verlet method）
int Leapfrog(int,double,double*);			// 蛙飛び法（Leap-frog scheme）
int RungeKutta(int,double,double*);
int RungeKuttaFehlberg(int,double,double,double*);	// Runge-Kutta-Fehlberg method
int PredictorCorrector(int,double,double*);		// predictor-corrector method (Adams-Bashforth-Moulton)
int modMidPoint(int,double,double*);			// 修正中点法（modified midpoint method）
int Ralston(int,double,double*);
int Kutta3o(int,double,double*);
int EulerRichardson(int,double,double*);
int RungeKutta4_38(int,double,double*);

extern int pn;
extern int mag;

#endif

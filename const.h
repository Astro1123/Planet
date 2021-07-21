#ifndef _CONST_H_
#define _CONST_H_

#define P_MAX		1000
#define CHAR_MAX	2048
#define CONST_G		0.0000000000667408	// [ m^3 kg^(-1) s^(-2) ]
#define UNIT_AU		149597870700		// [ m ]
#define UNIT_KM		1000			// [ m ]
#define UNIT_MIN	60			// [ s ]
#define UNIT_HOUR	60			// [ min ]
#define UNIT_DAY	24			// [ h ]



struct planet {
	double x[3];
	double m;
	double v[3];
	double vh[3];
	int vhf;
	int c[3];
	double r;
	char name[64]; 
};

extern struct planet data[];
#endif

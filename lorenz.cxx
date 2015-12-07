#include <iostream>
#include <cmath>

using namespace std;

void calculation(const int a, const int b, const double c, const double x, const double y, const double z, double k[3]);


int main() {
	const int a = 10;
	const int b = 28;
	const double c = 8.0/3.0;

	double r[3];
	r[0] = 1.0;			// x-coordinate
	r[1] = 1.0;			// y-coordinate
	r[2] = 2.0;			// z-coordinate

	const double dt = 0.001;
	const int t_end = 100;
	const int N = t_end/dt + 1;

	double k1[3];
	double k2[3];
	double k3[3];
	double k4[3];


	for(int i = 0; i < N; i++) {
		calculation(a, b, c, r[0], r[1], r[2], k1);
		calculation(a, b, c, r[0] + 0.5*dt*k1[0], r[1] + 0.5*dt*k1[1], r[2] + 0.5*dt*k1[2], k2);
		calculation(a, b, c, r[0] + 0.5*dt*k2[0], r[1] + 0.5*dt*k2[1], r[2] + 0.5*dt*k2[2], k3);
		calculation(a, b, c, r[0] + dt*k3[0], r[1] + dt*k3[1], r[2] + dt*k3[2], k4);

		for(int j = 0; j < 3; j++) {
			r[j] += dt/6.0 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
		}

		cout << i*dt << "\t" << r[0] << "\t" << r[1] << "\t" << r[2] << endl;
	}
	
	
	return 0;
}


void calculation(const int a, const int b, const double c, const double x, const double y, const double z, double k[3]) {
	k[0] = a*(y - x);
	k[1] = x*(b - z) - y;
	k[2] = x*y - c*z;
}

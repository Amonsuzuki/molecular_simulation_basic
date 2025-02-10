#include <bits/stdc++.h>
#include <matplot/matplot.h>
using namespace std;
using namespace matplot;

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

int GROUP = 4;
int X = 4;
int Y = 4;
int Z = 4;
int N = GROUP * X * Y * Z;
int STEP = 2000;
double RHO = 1.05;
double TEMP = 1;
double A = pow((4.0 / RHO), 1.0 / 3.0);
double LX = A * X;
double LY = A * Y;
double LZ = A * Z;
double Dr = 0.15;

double initial[4][3] = {{0, 0, 0}, {A/2, A/2, 0}, {0, A/2, A/2}, {A/2, 0, A/2}};

struct molecule{
	int property;
	double x;
	double y;
	double z;
	double newx;
	double newy;
	double newz;
	double mas = 1;
};

double calculate_force_and_potential(vector<molecule>& fcc){
	double eps = 1.0;
	double sigma = 1.0;
	double ce12 = 4 * eps * pow(sigma, 12);
	double ce6 = 4 * eps * pow(sigma, 6);
	double cf12 = ce12 * 12.0;
	double cf6 = ce6 * 6.0;
	double engp = 0;
	for (int i=0; i<N-1; i++) {
		for (int j=i+1; j<N; j++) {
			double x = fcc[j].x - fcc[i].x;
			double y = fcc[j].y - fcc[i].y;
			double z = fcc[j].z - fcc[i].z;
			if (x > LX / 2.0)
				x -= LX;
			if (x < -(LX / 2.0))
				x += LX;
			if (y > LY / 2.0)
				y -= LY;
			if (y < -(LY/ 2.0))
				y += LY;
			if (z > LZ / 2.0)
				z -= LZ;
			if (z < -(LZ / 2.0))
				z += LZ;
			double r2 = pow(x, 2) + pow(y, 2) + pow(z, 2);
			double r2i = 1 / r2;
			double r6i = pow(r2i, 3);
			double r12i = pow(r6i, 2);
			double ep = ce12 * r12i - ce6 *r6i;
			engp += ep;
		}
	}
	return(engp);
}

double calculate_du(vector<molecule>& fcc, int n_move){
	double eps = 1.0;
	double sigma = 1.0;
	double ce12 = 4 * eps * pow(sigma, 12);
	double ce6 = 4 * eps * pow(sigma, 6);
	double cf12 = ce12 * 12.0;
	double cf6 = ce6 * 6.0;
	double eng1 = 0;
	for (int i=0; i<N; i++) {
		if (i != n_move) {
			double x = fcc[n_move].x - fcc[i].x;
			double y = fcc[n_move].y - fcc[i].y;
			double z = fcc[n_move].z - fcc[i].z;
			if (x > LX / 2.0)
				x -= LX;
			if (x < -(LX / 2.0))
				x += LX;
			if (y > LY / 2.0)
				y -= LY;
			if (y < -(LY/ 2.0))
				y += LY;
			if (z > LZ / 2.0)
				z -= LZ;
			if (z < -(LZ / 2.0))
				z += LZ;
			double r2 = pow(x, 2) + pow(y, 2) + pow(z, 2);
			double r2i = 1 / r2;
			double r6i = pow(r2i, 3);
			double r12i = pow(r6i, 2);
			double ep = ce12 * r12i - ce6 *r6i;
			eng1 += ep;
		}
	}
	double eng2 = 0;
	for (int i=0; i<N; i++) {
		if (i != n_move) {
			double x = fcc[n_move].newx - fcc[i].x;
			double y = fcc[n_move].newy - fcc[i].y;
			double z = fcc[n_move].newz - fcc[i].z;
			if (x > LX / 2.0)
				x -= LX;
			if (x < -(LX / 2.0))
				x += LX;
			if (y > LY / 2.0)
				y -= LY;
			if (y < -(LY/ 2.0))
				y += LY;
			if (z > LZ / 2.0)
				z -= LZ;
			if (z < -(LZ / 2.0))
				z += LZ;
			double r2 = pow(x, 2) + pow(y, 2) + pow(z, 2);
			double r2i = 1 / r2;
			double r6i = pow(r2i, 3);
			double r12i = pow(r6i, 2);
			double ep = ce12 * r12i - ce6 *r6i;
			eng2 += ep;
		}
	}
	return(eng2 - eng1);
}

double initialize(vector<molecule>& fcc) {
	for (int i=0; i < N / GROUP; i++) {
		for (int j=0; j < GROUP; j++) {
			fcc[i * GROUP + j].property = j;
			fcc[i * GROUP + j].x = initial[j][0] + (i % X) * A;
			fcc[i * GROUP + j].y = initial[j][1] + (i / X) % Y * A;
			fcc[i * GROUP + j].z = initial[j][2] + (i / (X * Y)) % Z * A;
		}
	}
	return (calculate_force_and_potential(fcc));
}

double update_positions(vector<molecule>& fcc) {
	random_device rd;
	mt19937_64 mt64(rd());
	uniform_real_distribution<double> uni(0, 1);
	uniform_real_distribution<double> ikura(-1, 1);
	for (int i=0; i<N; i++) {
		int n_move = uni(mt64) * N;
		fcc[n_move].newx = fcc[n_move].newx + Dr * ikura(mt64);
		fcc[n_move].newy = fcc[n_move].newy + Dr * ikura(mt64);
		fcc[n_move].newz = fcc[n_move].newz + Dr * ikura(mt64);
		double du = calculate_du(fcc, n_move);
		if (du < 0) {
			fcc[n_move].x = fcc[n_move].newx;
			fcc[n_move].y = fcc[n_move].newy;
			fcc[n_move].z = fcc[n_move].newz;
			if(fcc[n_move].x > LX)
				fcc[n_move].x -= LX;
			else if (fcc[n_move].x < 0)
				fcc[n_move].x += LX;
			if(fcc[n_move].y > LY)
				fcc[n_move].y -= LY;
			else if (fcc[n_move].y < 0)
				fcc[n_move].y += LY;
			if(fcc[n_move].z > LZ)
				fcc[n_move].z -= LZ;
			else if (fcc[n_move].z < 0)
				fcc[n_move].z += LZ;
		}
		else {
			double p = exp(-du / TEMP);
			if (uni(mt64) < p) {
				fcc[n_move].x = fcc[n_move].newx;
				fcc[n_move].y = fcc[n_move].newy;
				fcc[n_move].z = fcc[n_move].newz;
				if(fcc[n_move].x > LX)
					fcc[n_move].x -= LX;
				else if (fcc[n_move].x < 0)
					fcc[n_move].x += LX;
				if(fcc[n_move].y > LY)
					fcc[n_move].y -= LY;
				else if (fcc[n_move].y < 0)
					fcc[n_move].y += LY;
				if(fcc[n_move].z > LZ)
					fcc[n_move].z -= LZ;
				else if (fcc[n_move].z < 0)
					fcc[n_move].z += LZ;
			}
		}
	}
	return (calculate_force_and_potential(fcc));
}

int main() {
	vector<molecule> fcc(N);
	vector<double> monte(STEP, 0), potential(STEP);

	potential[0] = initialize(fcc);

	for (int i=1; i<STEP; i++) {
		monte[i] = i;
		potential[i] = update_positions(fcc);
	}
	for (int i=0; i<STEP; i++) {
		cout << potential[i] << endl;
	}

	plt::named_plot("potential", monte, potential);
	plt::xlabel("t");
	plt::ylabel("Energy");

	plt::legend();
	plt::show();
	plt::pause(-1);

	vector<double> x, y, z;
	for (const auto& p : fcc) {
		x.push_back(p.x);
		y.push_back(p.y);
		z.push_back(p.z);
	}
//snapshots
//time dependency

/*
	scatter3(x, y, z);
	show();
	*/
	return (0);
}

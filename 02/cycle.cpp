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
double Dt = 0.005;
int STEP = 40000;
double RHO = 0.01;
double TEMP = 1;
double A = pow((4.0 / RHO), 1.0 / 3.0);
double LX = A * X;
double LY = A * Y;
double LZ = A * Z;

double initial[4][3] = {{0, 0, 0}, {A/2, A/2, 0}, {0, A/2, A/2}, {A/2, 0, A/2}};

struct molecule{
	int property;
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double mas = 1;
	double fx[2] = {0, 0};
	double fy[2] = {0, 0};
	double fz[2] = {0, 0};
};

double abso(double d) {
	return ((d < 0.0) ? -d : d);
}

double calculate_force_and_potential(vector<molecule>& fcc, int sign){
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
			double fc = (cf12 * r12i - cf6 * r6i) * r2i;
			double fx = fc * x;
			double fy = fc * y;
			double fz = fc * z;
			fcc[i].fx[sign] -= fx;
			fcc[i].fy[sign] -= fy;
			fcc[i].fz[sign] -= fz;
			fcc[j].fx[sign] += fx;
			fcc[j].fy[sign] += fy;
			fcc[j].fz[sign] += fz;
		}
	}
	return(engp);
}

double initialize(vector<molecule>& fcc, vector<double>& kinetic) {
	for (int i=0; i < N / GROUP; i++) {
		for (int j=0; j < GROUP; j++) {
			fcc[i * GROUP + j].property = j;
			fcc[i * GROUP + j].x = initial[j][0] + (i % X) * A;
			fcc[i * GROUP + j].y = initial[j][1] + (i / X) % Y * A;
			fcc[i * GROUP + j].z = initial[j][2] + (i / (X * Y)) % Z * A;
		}
	}
	mt19937_64 mt64(0);
	uniform_real_distribution<double> uni(-1, 1);
	double total_vx = 0;
	double total_vy = 0;
	double total_vz = 0;
	for (int i=0; i<N; i++) {
		fcc[i].vx = uni(mt64);
		fcc[i].vy = uni(mt64);
		fcc[i].vz = uni(mt64);
		total_vx += fcc[i].vx;
		total_vy += fcc[i].vy;
		total_vz += fcc[i].vz;
	}
	for (int i=0; i<N; i++) {
		fcc[i].vx -= total_vx / N;
		fcc[i].vy -= total_vy / N;
		fcc[i].vz -= total_vz / N;
	}
	double ke = 0;
	for (int i=0; i<N; i++)
		ke += 0.5 * (pow(fcc[i].vx, 2) + pow(fcc[i].vy, 2) + pow(fcc[i].vz, 2));
	ke /= N;
	double temp = ke / 1.5;
	double engk = 0;
	for (int i=0; i<N; i++) {
		fcc[i].vx *= sqrt(TEMP/temp);
		fcc[i].vy *= sqrt(TEMP/temp);
		fcc[i].vz *= sqrt(TEMP/temp);
		engk += 0.5 * pow(fcc[i].vx, 2);
		engk += 0.5 * pow(fcc[i].vy, 2);
		engk += 0.5 * pow(fcc[i].vz, 2);
	}
	kinetic[0] = engk;
	return (calculate_force_and_potential(fcc, 0));
}

double update_positions(vector<molecule>& fcc) {
	for (int i=0; i<N; i++) {
		fcc[i].x += Dt * fcc[i].vx + 0.5 * Dt * Dt / fcc[i].mas * fcc[i].fx[0];
		fcc[i].y += Dt * fcc[i].vy + 0.5 * Dt * Dt / fcc[i].mas * fcc[i].fy[0];
		fcc[i].z += Dt * fcc[i].vz + 0.5 * Dt * Dt / fcc[i].mas * fcc[i].fz[0];
		if(fcc[i].x > LX)
			fcc[i].x -= LX;
		else if (fcc[i].x < 0)
			fcc[i].x += LX;
		if(fcc[i].y > LY)
			fcc[i].y -= LY;
		else if (fcc[i].y < 0)
			fcc[i].y += LY;
		if(fcc[i].z > LZ)
			fcc[i].z -= LZ;
		else if (fcc[i].z < 0)
			fcc[i].z += LZ;
	}
	return (calculate_force_and_potential(fcc, 1));
}

double scaling_update_velocity(vector<molecule>& fcc) {
	double ke = 0;
	for (int i=0; i<N; i++)
		ke += 0.5 * (pow(fcc[i].vx, 2) + pow(fcc[i].vy, 2) + pow(fcc[i].vz, 2));
	ke /= N;
	double temp = ke / 1.5;
	double engk = 0;
	for (int i=0; i<N; i++) {
		fcc[i].vx += 0.5 * Dt / fcc[i].mas * (fcc[i].fx[0] + fcc[i].fx[1]);
		fcc[i].vy += 0.5 * Dt / fcc[i].mas * (fcc[i].fy[0] + fcc[i].fy[1]);
		fcc[i].vz += 0.5 * Dt / fcc[i].mas * (fcc[i].fz[0] + fcc[i].fz[1]);
		fcc[i].vx *= sqrt(TEMP/temp);
		fcc[i].vy *= sqrt(TEMP/temp);
		fcc[i].vz *= sqrt(TEMP/temp);
		engk += 0.5 * pow(fcc[i].vx, 2);
		engk += 0.5 * pow(fcc[i].vy, 2);
		engk += 0.5 * pow(fcc[i].vz, 2);
	}
	return (engk);
}

double update_velocity(vector<molecule>& fcc, double& temp_equil) {
	double ke = 0;
	for (int i=0; i<N; i++)
		ke += 0.5 * (pow(fcc[i].vx, 2) + pow(fcc[i].vy, 2) + pow(fcc[i].vz, 2));
	ke /= N;
	double temp = ke / 1.5;
	temp_equil += temp;
	double engk = 0;
	for (int i=0; i<N; i++) {
		fcc[i].vx += 0.5 * Dt / fcc[i].mas * (fcc[i].fx[0] + fcc[i].fx[1]);
		fcc[i].vy += 0.5 * Dt / fcc[i].mas * (fcc[i].fy[0] + fcc[i].fy[1]);
		fcc[i].vz += 0.5 * Dt / fcc[i].mas * (fcc[i].fz[0] + fcc[i].fz[1]);
		engk += 0.5 * pow(fcc[i].vx, 2);
		engk += 0.5 * pow(fcc[i].vy, 2);
		engk += 0.5 * pow(fcc[i].vz, 2);
	}
	return (engk);
}

void reset_forces(vector<molecule>& fcc) {
	for (int i=0; i<N; i++) {
		fcc[i].fx[0] = fcc[i].fx[1];
		fcc[i].fy[0] = fcc[i].fy[1];
		fcc[i].fz[0] = fcc[i].fz[1];
		fcc[i].fx[1] = 0;
		fcc[i].fy[1] = 0;
		fcc[i].fz[1] = 0;
	}
}

int main() {
	vector<molecule> fcc(N);
	vector<double> t(STEP, 0), kinetic(STEP), potential(STEP), total(STEP);
	double temp_equil = 0;
	double potential_equil = 0;

	potential[0] = initialize(fcc, kinetic);
	total[0] = kinetic[0] + potential[0];

	for (int i=1; i<STEP; i++) {
		t[i] = i * Dt;
		potential[i] = update_positions(fcc);
		if (i < STEP / 2)
			kinetic[i] = scaling_update_velocity(fcc);
		else {
			kinetic[i] = update_velocity(fcc, temp_equil);
			potential_equil += potential[i];
		}
		total[i] = kinetic[i] + potential[i];
		reset_forces(fcc);
	}
	for (int i=0; i<STEP; i++) {
		cout << kinetic[i] << ' ' << potential[i] << ' ' << total[i] << endl;
	}
	cout << "temp_equil: " << temp_equil * 2.0 / STEP << endl;
	cout << "potential_equil: " << potential_equil * 2.0 / STEP << endl;

	plt::named_plot("kinetic", t, kinetic);
	plt::named_plot("potential", t, potential);
	plt::named_plot("total", t, total);
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
//potential_equil
//snapshots
//time dependency

/*
	scatter3(x, y, z);
	show();
	*/
	return (0);
}

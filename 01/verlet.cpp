#include <bits/stdc++.h>
#include "matplotlibcpp.h"
using namespace std;
namespace plt = matplotlibcpp;

double M = 1.0;
double K = 1.0;
double N = 70;
double H = 0.1;

struct plot{
	double t;
	double x = 1;
	double v = 0;
	double ans_x = 1;
	double ans_v = 0;
	double kinetic = 1.0 / 2.0 * M * pow(v, 2);
	double urgentia = 1.0 / 2.0 * K * pow(x, 2);
	double total = 1.0 / 2.0 * K * pow(x, 2);
	void update() {
		kinetic = 1.0 / 2.0 * M * pow(v, 2);
		urgentia = 1.0 / 2.0 * K * pow(x, 2);
		total = kinetic + urgentia;
	}
};

int main() {
	vector<plot> verlet(N);
	for (int i=1; i<N; i++) {
		verlet[i].t = i * H;
		verlet[i].x = verlet[i - 1].x + H * verlet[i - 1].v - pow(H, 2) / 2.0 * (K/M) * verlet[i - 1].x;
		verlet[i].v = verlet[i - 1].v - H / 2 * (K/M) * (verlet[i - 1].x + verlet[i].x);
		verlet[i].update();
	}
	for (int i=1; i<N; i++) {
		verlet[i].ans_x = cos(sqrt(K/M) * verlet[i].t);
		verlet[i].ans_v = -sin(sqrt(K/M) * verlet[i].t);
	}
	vector<double> t(N), x(N), v(N), ans_x(N), ans_v(N), kinetic(N), urgentia(N), total(N);
	transform(verlet.begin(), verlet.end(), t.begin(), [](const plot& p) { return p.t; }); 
	transform(verlet.begin(), verlet.end(), x.begin(), [](const plot& p) { return p.x; }); 
	transform(verlet.begin(), verlet.end(), v.begin(), [](const plot& p) { return p.v; }); 
	transform(verlet.begin(), verlet.end(), ans_x.begin(), [](const plot& p) { return p.ans_x; }); 
	transform(verlet.begin(), verlet.end(), ans_v.begin(), [](const plot& p) { return p.ans_v; }); 
	transform(verlet.begin(), verlet.end(), kinetic.begin(), [](const plot& p) { return p.kinetic; }); 
	transform(verlet.begin(), verlet.end(), urgentia.begin(), [](const plot& p) { return p.urgentia; }); 
	transform(verlet.begin(), verlet.end(), total.begin(), [](const plot& p) { return p.total; }); 
/*
	plt::scatter(t, x, 10.0);
	plt::scatter(t, v, 10.0);
	plt::plot(t, ans_x, "C2-");
	plt::plot(t, ans_v, "C3-");
	plt::named_plot("x", t, x);
	plt::named_plot("v", t, v);
	plt::named_plot("x answer", t, ans_x);
	plt::named_plot("v answer", t, ans_v);
	plt::xlabel("t");
	plt::ylabel("Value");
*/
	plt::scatter(t, kinetic, 10.0);
	plt::scatter(t, urgentia, 10.0);
	plt::scatter(t, total, 10.0);
	plt::named_plot("kinetic", t, kinetic);
	plt::named_plot("urgentia", t, urgentia);
	plt::named_plot("total", t, total);
	plt::xlabel("t");
	plt::ylabel("Energy");

	plt::legend();
	plt::show();
	plt::pause(-1);
	return (0);
}

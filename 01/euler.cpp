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
	vector<plot> euler(N);
	for (int i=1; i<N; i++) {
		euler[i].t = i * H;
		euler[i].x = euler[i - 1].x + H * euler[i - 1].v;
		euler[i].v = euler[i - 1].v - H * K / M * euler[i - 1].x;
		euler[i].update();
	}
	for (int i=1; i<N; i++) {
		euler[i].ans_x = cos(sqrt(K/M) * euler[i].t);
		euler[i].ans_v = -sin(sqrt(K/M) * euler[i].t);
	}
	vector<double> t(N), x(N), v(N), ans_x(N), ans_v(N), kinetic(N), urgentia(N), total(N);
	transform(euler.begin(), euler.end(), t.begin(), [](const plot& p) { return p.t; }); 
	transform(euler.begin(), euler.end(), x.begin(), [](const plot& p) { return p.x; }); 
	transform(euler.begin(), euler.end(), v.begin(), [](const plot& p) { return p.v; }); 
	transform(euler.begin(), euler.end(), ans_x.begin(), [](const plot& p) { return p.ans_x; }); 
	transform(euler.begin(), euler.end(), ans_v.begin(), [](const plot& p) { return p.ans_v; }); 
	transform(euler.begin(), euler.end(), kinetic.begin(), [](const plot& p) { return p.kinetic; }); 
	transform(euler.begin(), euler.end(), urgentia.begin(), [](const plot& p) { return p.urgentia; }); 
	transform(euler.begin(), euler.end(), total.begin(), [](const plot& p) { return p.total; }); 
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

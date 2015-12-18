#pragma once
#include <emmintrin.h>

class Particle {
public:
	Particle(double x, double y, double velocity_x, double velocity_y, double mass, double acceleration_x, double acceleration_y) :
		x_(x), y_(y), velocity_x_(velocity_x), velocity_y_(velocity_y), mass_(mass), acceleration_x_(acceleration_x), acceleration_y_(acceleration_y) {};

	void add_acceleration_pairwise(Particle& interacting_particle);
	void add_acceleration(const Particle& interacting_particle);
	void advance(double time_stamp);
	double x_;
	double y_;
	double mass_;
	double velocity_x_, velocity_y_;
	double acceleration_x_, acceleration_y_;


	// TODO: check intrinsics
	/*__m128d position;

	__m128d interacting_particle = _mm_setr_pd(interacting_particle.m128d_f64[0], interacting_particle.y_);
	auto interacting_xy = _mm_setr_pd(1.5, 2.5);
	auto addition3 = _mm_add_sd(xy, interacting_xy);
	auto added_double = addition3.m128d_f64[0];*/
};

#pragma once
class Particle {
public:
	Particle(double x, double y, double velocity_x, double velocity_y, double mass, double force_x, double force_y) :
		x_(x), y_(y), velocity_x_(velocity_x), velocity_y_(velocity_y), mass_(mass), force_x_(force_x), force_y_(force_y) {};
	void add_force(const Particle& input_particle);
	void advance(double time_stamp);
	void reset_forces();

	double x_, y_;
	double velocity_x_, velocity_y_;
	double mass_;
	double force_x_;
	double force_y_;
};
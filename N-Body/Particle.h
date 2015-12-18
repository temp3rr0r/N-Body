#pragma once
class Particle {
public:
	Particle(double x, double y, double velocity_x, double velocity_y, double mass, double acceleration_x, double acceleration_y) :
		x_(x), y_(y), velocity_x_(velocity_x), velocity_y_(velocity_y), mass_(mass), acceleration_x_(acceleration_x), acceleration_y_(acceleration_y) {};
	void add_acceleration(const Particle& interacting_particle);
	void advance(double time_stamp);

	double x_, y_;
	double velocity_x_, velocity_y_;
	double acceleration_x_, acceleration_y_;
	double mass_;
};
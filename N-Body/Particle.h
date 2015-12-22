#pragma once

class Particle {
public:
	double x_;
	double y_;
	double mass_;
	double velocity_x_, velocity_y_;
	double acceleration_x_, acceleration_y_;

	Particle() {
		x_ = 0.0;
		y_ = 0.0;
		velocity_x_ = 0.0;
		velocity_y_ = 0.0;
		mass_ = 0.0;
		acceleration_x_ = 0.0;
		acceleration_y_ = 0.0;
	}

	Particle(double x, double y, double velocity_x, double velocity_y, double mass, double acceleration_x, double acceleration_y) :
		x_(x), y_(y), mass_(mass), velocity_x_(velocity_x), velocity_y_(velocity_y), acceleration_x_(acceleration_x), acceleration_y_(acceleration_y) {};

	void add_acceleration_pairwise(Particle& interacting_particle);
	double get_distance(const Particle& second_particle) const;
	void add_acceleration(const Particle& interacting_particle);
	void advance(double time_stamp);
	Particle operator+(const Particle& r) const;
	Particle operator-(const Particle& r) const;
	Particle operator*(float r) const;
};

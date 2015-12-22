#pragma once

// Class that stores information for every particle
class Particle {
public:
	float x_;
	float y_;
	float mass_;
	float velocity_x_, velocity_y_;
	float acceleration_x_, acceleration_y_;

	// Default constructor
	Particle() {
		x_ = 0.0;
		y_ = 0.0;
		velocity_x_ = 0.0;
		velocity_y_ = 0.0;
		mass_ = 0.0;
		acceleration_x_ = 0.0;
		acceleration_y_ = 0.0;
	}

	// Constructor useful for center of mass particles
	Particle(float x, float y, float mass) :
		x_(x), y_(y), mass_(mass) {
		velocity_x_ = 0.0;
		velocity_y_ = 0.0;
		acceleration_x_ = 0.0;
		acceleration_y_ = 0.0;
	};

	// Full constructor
	Particle(float x, float y, float velocity_x, float velocity_y, float mass,
		float acceleration_x, float acceleration_y) :
		x_(x), y_(y), mass_(mass), velocity_x_(velocity_x), velocity_y_(velocity_y),
		acceleration_x_(acceleration_x), acceleration_y_(acceleration_y) {};

	void add_acceleration_pairwise(Particle& interacting_particle);
	float get_distance(const Particle& second_particle) const;
	void add_acceleration(float total_mass, float center_of_mass_x, float center_of_mass_y);
	void add_acceleration(const Particle& interacting_particle);
	void advance(float time_stamp);
	Particle operator+(const Particle& r) const;
	Particle operator-(const Particle& r) const;
	Particle operator*(float r) const;
};

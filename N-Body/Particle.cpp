
#include "Particle.h"
#include <cmath>
#include "Settings.h"

// Apply acceleration on both particles in one sweep
void Particle::add_acceleration_pairwise(Particle& interacting_particle) {

	// We calculate the acceleration instead of the forces for faster calculation

	// Get distances
	double dx = interacting_particle.x_ - x_;
	double dy = interacting_particle.y_ - y_;

	// Square of distances
	double distance_square = dx * dx + dy * dy;

	// Keep a minimum square of distance
	if (distance_square < MIN_DISTANCE)
		distance_square = MIN_DISTANCE;

	double distance = sqrt(distance_square);

	velocity_x_ = dx / distance;
	velocity_y_ = dy / distance;

	// Use the mass of the other particle
	double acceleration_factor = (GRAVITATIONAL_CONSTANT / distance_square) * interacting_particle.mass_;
	double interacting_acceleration_factor = (GRAVITATIONAL_CONSTANT / distance_square) * mass_;

	// Apply accelerations
	acceleration_x_ -= acceleration_factor * velocity_x_;
	acceleration_y_ -= acceleration_factor * velocity_y_;
	
	interacting_particle.acceleration_x_ += interacting_acceleration_factor * velocity_x_;
	interacting_particle.acceleration_y_ += interacting_acceleration_factor * velocity_y_;
}

// Aquire the distance of the particle from another particle
double Particle::get_distance(const Particle& second_particle) const {
	// Get distances
	double dx = x_ - second_particle.x_;
	double dy = y_ - second_particle.y_;

	// Square of distances
	double distance_square = dx * dx + dy * dy;

	// Keep a minimum square of distance
	if (distance_square < MIN_DISTANCE)
		distance_square = MIN_DISTANCE;

	return sqrt(distance_square);
}

// Apply acceleration on one particle from a center of mass
void Particle::add_acceleration(double total_mass, double center_of_mass_x, double center_of_mass_y) {

	// Get distances
	double dx = center_of_mass_x - x_;
	double dy = center_of_mass_y - y_;

	// Square of distances
	double distance_square = dx * dx + dy * dy;

	// Keep a minimum square of distance
	if (distance_square < MIN_DISTANCE)
		distance_square = MIN_DISTANCE;

	double distance = sqrt(distance_square);

	velocity_x_ = dx / distance;
	velocity_y_ = dy / distance;

	double acceleration_factor = (GRAVITATIONAL_CONSTANT / distance_square) * total_mass;

	// Apply accelerations
	acceleration_x_ -= acceleration_factor * velocity_x_;
	acceleration_y_ -= acceleration_factor * velocity_y_;
}

// Apply acceleration on one particle from forces of another particle
void Particle::add_acceleration(const Particle& interacting_particle) {

	// Get distances
	double dx = interacting_particle.x_ - x_;
	double dy = interacting_particle.y_ - y_;

	// Square of distances
	double distance_square = dx * dx + dy * dy;
	
	// Keep a minimum square of distance
	if (distance_square < MIN_DISTANCE)
		distance_square = MIN_DISTANCE;
	
	double distance = sqrt(distance_square);

	velocity_x_ = dx / distance;
	velocity_y_ = dy / distance;		

	double acceleration_factor = (GRAVITATIONAL_CONSTANT / distance_square) * interacting_particle.mass_;

	// Apply accelerations
	acceleration_x_ -= acceleration_factor * velocity_x_;
	acceleration_y_ -= acceleration_factor * velocity_y_;
}

// Moves the current particle for a specific time
void Particle::advance(double time_step) {

	// Add accelerations on velocities
	velocity_x_ += time_step * acceleration_x_;
	velocity_y_ += time_step * acceleration_y_;

	// Get the new position
	x_ += velocity_x_ * time_step;
	y_ += velocity_y_ * time_step;

	// If out of grid limits, reverse direction and set valid position
	if (x_ < 0) {
		velocity_x_ *= -1;
		x_ *= -1;
	} else if (x_ > UNIVERSE_SIZE_X) {
		velocity_x_ *= -1;
		x_ -= UNIVERSE_SIZE_X;
	}
	
	if (y_ < 0) {
		velocity_y_ *= -1;
		y_ *= -1;
	} else if (y_ > UNIVERSE_SIZE_Y) {
		velocity_y_ *= -1;
		y_ -= UNIVERSE_SIZE_Y;
	}

	// Reset accelerations
	acceleration_x_ = 0.0;
	acceleration_y_ = 0.0;
}

// Combine two particles into one by adding them
Particle Particle::operator+(const Particle& r) const {
	return Particle(x_ + r.x_, y_ + r.y_, velocity_x_ + r.velocity_x_, velocity_y_ + r.velocity_y_, mass_ + r.mass_,
	                acceleration_x_ + r.acceleration_x_, acceleration_y_ + r.acceleration_y_);
}

// Combine two particles into one by substracting them
Particle Particle::operator-(const Particle& r) const {
	return Particle(x_ - r.x_, y_ - r.y_, velocity_x_ - r.velocity_x_, velocity_y_ - r.velocity_y_, mass_ - r.mass_,
	                acceleration_x_ - r.acceleration_x_, acceleration_y_ - r.acceleration_y_);
}

// Multiply a particle by a scalar number
Particle Particle::operator*(double r) const {
	return Particle(x_ * r, y_ * r, velocity_x_ * r, velocity_y_ * r, mass_ * r, acceleration_x_ * r, acceleration_y_ * r);
}

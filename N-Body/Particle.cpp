
#include "Particle.h"
#include <cmath>
#include "Settings.h"

void Particle::add_acceleration(Particle& interacting_particle) {

	// Get distance
	double dx = interacting_particle.x_ - x_;
	double dy = interacting_particle.y_ - y_;

	// Squares of distances
	double distance_square = dx * dx + dy * dy;
	
	// Keep a minimum square of distance
	if (distance_square < MIN_DISTANCE)
		distance_square = MIN_DISTANCE;

	double distance = sqrt(distance_square);

	velocity_x_ = dx / distance;
	velocity_y_ = dy / distance;
		
	double Gdivd = GRAVITATIONAL_CONSTANT / distance_square;

	double acceleration_factor = Gdivd * mass_;
	double interacting_acceleration_factor = Gdivd * interacting_particle.mass_;

	// Apply accelerations

	acceleration_x_ -= acceleration_factor * velocity_x_;
	acceleration_y_ -= acceleration_factor * velocity_y_;

	interacting_particle.acceleration_x_ += interacting_acceleration_factor * velocity_x_;
	interacting_particle.acceleration_y_ += interacting_acceleration_factor * velocity_y_;
}

void Particle::advance(double time_step) {

	// Add accelerations on velocities
	velocity_x_ += time_step * acceleration_x_;
	velocity_y_ += time_step * acceleration_y_;

	// Get the new position
	x_ += velocity_x_ * time_step;
	y_ += velocity_y_ * time_step;

	// If out of points, reverse direction
	if (x_ < 0) {
		velocity_x_ *= -1;
		x_ *= -1;
	} else if (x_ > UNIVERSE_SIZE_X) {
		velocity_x_ *= -1;
		x_ = UNIVERSE_SIZE_X - x_;
	}
	
	if (y_ < 0) {
		velocity_y_ *= -1;
		y_ *= -1;
	} else if (y_ > UNIVERSE_SIZE_Y) {
		velocity_y_ *= -1;
		y_ = UNIVERSE_SIZE_Y - y_;
	}

	// Reset accelerations
	acceleration_x_ = 0.0;
	acceleration_y_ = 0.0;
}

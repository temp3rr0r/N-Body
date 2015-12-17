
#include "Particle.h"
#include <cmath>
#include "Settings.h"

void Particle::add_force(const Particle& interacting_particle) {

	double EPS = 3E4;
	double dx = interacting_particle.x_ - x_;
	double dy = interacting_particle.y_ - y_;
	double distance = sqrt(dx*dx + dy*dy);
	double F = (GRAVITY * mass_ * interacting_particle.mass_) / (distance * distance + EPS*EPS);
	force_x_ += F * dx / distance;
	force_y_ += F * dy / distance;
}

void Particle::advance(double time_step) {
	velocity_x_ += time_step * force_x_ / mass_;
	velocity_y_ += time_step * force_y_ / mass_;

	x_ += time_step * force_x_ / mass_;
	if (x_ < 0 || x_ > UNIVERSE_SIZE_X) {
		velocity_x_ *= -1;
		x_ *= -1;
	} else if (x_ > UNIVERSE_SIZE_X) {
		velocity_x_ *= -1;
		x_ -= UNIVERSE_SIZE_X;
	}
	
	y_ += time_step * force_y_ / mass_;	
	if (y_ < 0 || y_ > UNIVERSE_SIZE_Y) {
		velocity_y_ *= -1;
		y_ *= -1;
	} else if (y_ > UNIVERSE_SIZE_Y) {
		velocity_y_ *= -1;
		y_ -= UNIVERSE_SIZE_Y;
	}
}

void Particle::reset_forces() {
	force_x_ = 0.0;
	force_y_ = 0.0;
}

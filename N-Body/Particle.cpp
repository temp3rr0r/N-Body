
#include "Particle.h"
#include <cmath>
#include "Settings.h"

inline void Particle::add_force(Particle& interacting_particle) {

	double EPS = 3E4;
	double dx = interacting_particle.x_ - x_;
	double dy = interacting_particle.y_ - y_;
	double distance = sqrt(dx*dx + dy*dy);
	double F = (GRAVITY * mass_ * interacting_particle.mass_) / (distance * distance + EPS*EPS);
	force_x_ += F * dx / distance;
	force_y_ += F * dy / distance;
}

inline void Particle::advance(double time_step) {
	velocity_x_ += time_step * force_x_ / mass_;
	velocity_y_ += time_step * force_y_ / mass_;

	x_ += time_step * force_x_ / mass_;
	y_ += time_step * force_y_ / mass_;
}

inline void Particle::reset_forces() {
	force_x_ = 0.0;
	force_y_ = 0.0;
}

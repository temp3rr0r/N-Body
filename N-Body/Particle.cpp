
#include "Particle.h"
#include <cmath>
#include "Settings.h"

void Particle::add_acceleration(Particle& interacting_particle) {

//	double EPS = 3E4;
//	double dx = interacting_particle.x_ - x_;
//	double dy = interacting_particle.y_ - y_;
//	double distance = sqrt(dx*dx + dy*dy);
//	double F = (GRAVITY * mass_ * interacting_particle.mass_) / (distance * distance + EPS*EPS);
//	force_x_ += F * dx / distance;
//	force_y_ += F * dy / distance;

	double dx = interacting_particle.x_ - x_;
	double dy = interacting_particle.y_ - y_;
	double distsq = dx*dx + dy*dy;
	
	if (distsq < MIN_DISTANCE)
		distsq = MIN_DISTANCE;
	double dist = sqrt(distsq);

	velocity_x_ = dx / dist;
	velocity_y_ = dy / dist;

	double GFORCE = 1.0;
	double Gdivd = GFORCE / distsq;

	double ai = Gdivd * mass_;
	double aj = Gdivd * interacting_particle.mass_;

	acceleration_x_ -= ai * velocity_x_;
	acceleration_y_ -= ai * velocity_y_;

	interacting_particle.acceleration_x_ += aj * velocity_x_;
	interacting_particle.acceleration_y_ += aj * velocity_y_;
}

void Particle::advance(double time_step) {

	// Add accelerations on velocities
	velocity_x_ += time_step * acceleration_x_;
	velocity_y_ += time_step * acceleration_y_;

	x_ += velocity_x_ * time_step;
	y_ += velocity_y_ * time_step;

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

	acceleration_x_ = 0.0;
	acceleration_y_ = 0.0;
}

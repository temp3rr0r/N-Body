#pragma once
#include "Particle.h"

// A class wrapper of typical particles for a Quad tree particle collection
class TreeParticle {
	Particle particle_;
public:
	TreeParticle() { }
	TreeParticle(const Particle& input_particle) : particle_(input_particle) { }
	
	const Particle& get_particle() const {
		return particle_;
	}
	
	void set_particle(const Particle& input_particle) {
		particle_ = input_particle;
	}

	float get_mass() const {
		return particle_.mass_;
	}
};
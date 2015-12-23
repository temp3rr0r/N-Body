#pragma once

#include <vector>
#include "TreeParticle.h"
#include "Particle.h"

// A Quad tree that stores collections of particles
class QuadParticleTree {
	static const uint8_t NUM_CHILDREN = 4; // Quad trees have 4 children
	Particle origin;	// The bounding box of this node
	Particle halfDimension; // Helps in spliting to children quadrants
	QuadParticleTree *children[NUM_CHILDREN]; // Pointers to the children quadrets
	TreeParticle *data; // Temporary placeholder

	float center_of_mass_x_ = 0.0;
	float center_of_mass_y_ = 0.0;
	float total_mass_ = 0.0;
	uint8_t depth = 0; // The depth of the node compared to the root
public:
	QuadParticleTree(const Particle& origin, const Particle& halfDimension);
	float get_side_size() const;
	float get_total_mass() const;
	~QuadParticleTree();
	int get_quadrant_containing_point(const Particle& point) const; // Find the child node quadrant
	bool isLeafNode() const; // Check if it is a leaf
	void insert(TreeParticle* point); // Insert point in the node
	void apply_acceleration(Particle& input_particle) const;
};
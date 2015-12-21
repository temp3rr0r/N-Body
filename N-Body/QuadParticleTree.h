#pragma once

#include <vector>
#include "TreeParticle.h"
#include "Particle.h"

class QuadParticleTree {
	static const uint8_t NUM_CHILDREN = 4;
	Particle origin;
	Particle halfDimension;
	QuadParticleTree *children[NUM_CHILDREN];
	TreeParticle *data;
	bool hasChildren;

	float center_of_mass_x_ = 0.0;
	float center_of_mass_y_ = 0.0;
	float total_mass_ = 0.0;
	uint8_t depth = 0;
public:
	QuadParticleTree(const Particle& origin, const Particle& halfDimension);
	QuadParticleTree(const QuadParticleTree& copy);
	float get_total_mass() const;
	~QuadParticleTree();
	int getQuadrantContainingPoint(const Particle& point) const;
	bool isLeafNode() const;
	void insert(TreeParticle* point);
	void getPointsInsideBox(const Particle& bmin, const Particle& bmax, std::vector<TreeParticle*>& results) const;
};
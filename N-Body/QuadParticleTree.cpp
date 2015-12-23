#pragma once
#include "Particle.h"
#include "QuadParticleTree.h"
#include "Settings.h"

inline QuadParticleTree::QuadParticleTree(const Particle& origin, const Particle& halfDimension) : origin(origin), halfDimension(halfDimension), data(nullptr) {
	// Ensure that the children are empty on the new node
	for (int i = 0; i < NUM_CHILDREN; ++i)
		children[i] = nullptr;
}

inline QuadParticleTree::QuadParticleTree(const QuadParticleTree& copy) : origin(copy.origin), halfDimension(copy.halfDimension), data(copy.data) {
}

float QuadParticleTree::get_side_size() const {
	return origin.get_distance(halfDimension);
}

inline float QuadParticleTree::get_total_mass() const {
	return total_mass_;
}

// Delete quadrants recursively
QuadParticleTree::~QuadParticleTree() {
	for (int i = 0; i < NUM_CHILDREN; ++i)
		delete children[i];
}

// Find which quandrant contains the point
//x : --++
//y : -+-+
inline int QuadParticleTree::getQuadrantContainingPoint(const Particle& point) const {
	int quadrant = 0;
	if (point.x_ >= origin.x_)
		quadrant |= 2;
	if (point.y_ >= origin.y_)
		quadrant |= 1;
	return quadrant;
}

inline bool QuadParticleTree::isLeafNode() const {
	// If this node is a leaf, then at least the first child will be null
	return children[0] == nullptr;
}

void QuadParticleTree::insert(TreeParticle* point) {

	if (isLeafNode()) {
		// If the node is a leaf, we don't have to "dive" any depper in the tree

		// Apply center of mass calculations
		// TODO: what if there are already children in the leaf?
		this->total_mass_ += point->get_mass();
		this->center_of_mass_x_ += point->get_particle().x_;
		this->center_of_mass_y_ += point->get_particle().y_;

		if (data == nullptr) {
			// This leaf has no data already, just store the data
			data = point;
		} else {
			// We are at a leaf but there are data already stored

			if (this->depth < MAX_TREE_DEPTH) { // TODO: switch to vector, to be able to add extra children on a terminal node
				TreeParticle *oldPoint = data; // Temporary placeholder
				data = nullptr;

				// Split the node to children
				for (uint8_t i = 0; i < NUM_CHILDREN; ++i) {
					// Find the new bounding box for the child
					Particle newOrigin = origin;
					newOrigin.x_ += halfDimension.x_ * (i & 2 ? .5f : -.5f);
					newOrigin.y_ += halfDimension.y_ * (i & 1 ? .5f : -.5f);
					children[i] = new QuadParticleTree(newOrigin, halfDimension *.5f);
					// Increase the node depth
					children[i]->depth = depth + 1;
				}

				// Re-insert the older data of the leaf into the correct child
				children[getQuadrantContainingPoint(oldPoint->get_particle())]->insert(oldPoint);
				// Continue to recursively find were to insert the requested point. We don't have
				// to search from the tree root. We can continue from the current node.
				children[getQuadrantContainingPoint(point->get_particle())]->insert(point);
			}
		}
	} else {
		// We are in an interrior node. We have to go deeper into the appropriate child quadrant

		this->total_mass_ += point->get_mass();

		float center_y = 0.0;
		float center_x = 0.0;
		uint8_t children_count = 0;
		if (children != nullptr) {
			for (uint8_t i = 0; i < NUM_CHILDREN; ++i) {
				if (children[i] != nullptr) {
					if (children[i]->data != nullptr) {
						++children_count;
						center_x += children[i]->data->get_particle().x_ * children[i]->data->get_mass();
						center_y += children[i]->data->get_particle().y_ * children[i]->data->get_mass();
					}
				}
			}
		}
		if (children_count > 0) {
			this->center_of_mass_x_ = center_x / this->total_mass_;
			this->center_of_mass_y_ = center_y / this->total_mass_;
		}

		int quadrant = getQuadrantContainingPoint(point->get_particle());

		children[quadrant]->insert(point);
	}
}

Particle& QuadParticleTree::get_center_of_mass_particle(Particle& input_particle) const {
	// Start from root
	if (isLeafNode()) {
		Particle center_of_mass_particle = Particle(center_of_mass_x_, center_of_mass_y_, total_mass_);
		return center_of_mass_particle;
	} else {
		// Get distances	
		Particle center_of_mass_particle = Particle(center_of_mass_x_, center_of_mass_y_, total_mass_);
		float distance_from_center_of_mass = input_particle.get_distance(center_of_mass_particle);
		float side = get_side_size();

		if (side / distance_from_center_of_mass < THETA) {
			return center_of_mass_particle;
		}
		else {
			// Go deeper in the tree
			int quadtrant = getQuadrantContainingPoint(input_particle);
			children[quadtrant]->get_center_of_mass_particle(input_particle);
		}
	}
}

void QuadParticleTree::apply_acceleration(Particle& input_particle) const {
	// Start from root
	if (isLeafNode()) {
		Particle center_of_mass_particle = Particle(center_of_mass_x_, center_of_mass_y_, total_mass_);
		if (input_particle.x_ != center_of_mass_particle.x_ && input_particle.y_ != center_of_mass_particle.y_ && input_particle.mass_ != total_mass_)
			input_particle.add_acceleration(center_of_mass_particle);	
	} else {
		// Get distances	
		Particle center_of_mass_particle = Particle(center_of_mass_x_, center_of_mass_y_, total_mass_);
		float distance_from_center_of_mass = input_particle.get_distance(center_of_mass_particle);
		float side = get_side_size();

		if (side / distance_from_center_of_mass < THETA) {
			if (input_particle.x_ != center_of_mass_particle.x_ && input_particle.y_ != center_of_mass_particle.y_ && input_particle.mass_ != total_mass_)
				input_particle.add_acceleration(center_of_mass_particle);
		} else {
			// Go deeper in the tree
			int quadtrant = getQuadrantContainingPoint(input_particle);
			children[quadtrant]->apply_acceleration(input_particle);
		}
	}
}

// Recursively gather all the TreeParticles of the node
void QuadParticleTree::get_all_particles(size_t size_x, size_t size_y, std::vector<TreeParticle*>& results) const {

	// Starting and ending bounding box are the Grid size X and size Y
	get_points_inside_box(Particle(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0),
		Particle(static_cast<float>(size_x), static_cast<float>(size_y), 0.0f, 0.0f, 0.0f, 0.0f, 0.0), results);
};


// Recursively gather all the TreeParticles that are within a requested bounding box
inline void QuadParticleTree::get_points_inside_box(const Particle& bmin, const Particle& bmax, std::vector<TreeParticle*>& results) const {

	if (isLeafNode()) {
		if (data != nullptr) {
			const Particle& p = data->get_particle();
			if (p.x_ > bmax.x_ || p.y_ > bmax.y_)
				return;
			if (p.x_ < bmin.x_ || p.y_ < bmin.y_)
				return;
			results.push_back(data);
		}
	}
	else {
		for (uint8_t i = 0; i < NUM_CHILDREN; ++i) {
			Particle cmax = children[i]->origin + children[i]->halfDimension;
			Particle cmin = children[i]->origin - children[i]->halfDimension;

			if (cmax.x_ < bmin.x_ || cmax.y_  <bmin.y_)
				continue;
			if (cmin.x_ > bmax.x_ || cmin.y_ > bmax.y_)
				continue;

			children[i]->get_points_inside_box(bmin, bmax, results);
		}
	}
}

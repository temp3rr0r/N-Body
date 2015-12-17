#include <random>
#include "Settings.h"
#include "ParticleHandler.h"
#include <tbb/concurrent_vector.h>
#include "lodepng.h"

void ParticleHandler::allocate_random_particles(size_t particle_count, std::vector<Particle>& particles, size_t size_x, size_t size_y) {
	if (particle_count > 0) {
		std::random_device random_device_;
		std::mt19937 mt_engine(random_device_());
		std::uniform_real_distribution<> real_value(0, MAX_RANDOM);
		/*std::uniform_real_distribution<> real_position_x(0, size_x);
		std::uniform_real_distribution<> real_position_y(0, size_y);*/
		std::uniform_real_distribution<> real_position_x(0.4 * size_x, 0.6 * size_x);
		std::uniform_real_distribution<> real_position_y(0.4 * size_y, 0.6 * size_y);
		std::uniform_real_distribution<> real_mass(0, MAX_MASS);

		for (size_t i = 0; i < particle_count; ++i)
			particles.push_back(Particle(real_position_x(mt_engine), real_position_y(mt_engine), real_value(mt_engine),
				real_value(mt_engine), real_value(mt_engine), real_value(mt_engine), real_mass(mt_engine)));
	}
}

void ParticleHandler::universe_to_png(const std::vector<Particle>& universe, size_t universe_size_x, size_t universe_size_y, const char* filename) {
	// Generate the image
	uint32_t width = static_cast<uint32_t>(universe_size_y);
	uint32_t height = static_cast<uint32_t>(universe_size_x);
	std::vector<uint8_t> image;
	image.resize(width * height * 4);

	// Paint all black
	for (uint32_t x = 0; x < height; x++) {
		for (uint32_t y = 0; y < width; y++) {
			image[4 * width * x + 4 * y + 0] = 0;
			image[4 * width * x + 4 * y + 1] = 0;
			image[4 * width * x + 4 * y + 2] = 0;
			image[4 * width * x + 4 * y + 3] = 255;
		}
	}

	for (const Particle& current_particle : universe) {
		image[4 * width * static_cast<uint32_t>(current_particle.x_) + 4 * static_cast<uint32_t>(current_particle.y_) + 0] = 255; // red
	}

	lodepng::encode(filename, image, width, height);
}


tbb::concurrent_vector<Particle> ParticleHandler::to_concurrent_vector(const std::vector<Particle>& input_particles) {
	tbb::concurrent_vector<Particle, tbb::cache_aligned_allocator<Particle>> returning_concurrent_vector;

	for (Particle current_particle : input_particles)
		returning_concurrent_vector.push_back(current_particle);
	return returning_concurrent_vector;
}

std::vector<Particle> ParticleHandler::to_vector(const tbb::concurrent_vector<Particle>& input_particles) {
	std::vector<Particle> returning_vector;

	for (Particle current_particle : input_particles)
		returning_vector.push_back(current_particle);

	return returning_vector;
}
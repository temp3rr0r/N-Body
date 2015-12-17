// N-Body.cpp : Defines the entry point for the console application.
//


#include "Settings.h"
#include <ostream>
#include <iostream>
#include <vector>
#include "tbb/tick_count.h"
#include "Particle.h"
#include "tbb/cache_aligned_allocator.h"
#include "tbb/concurrent_vector.h"
#include "ParticleHandler.h"

void simulate_serial(std::vector<Particle>& particles, double total_time_steps, double time_step, size_t particle_count) {

	//for (double iteration = 0; iter)
	//for (int count = 0; count < numberofiterations; ++count) {
	//	for (int i = 0; i < N; ++i) {
	//		particles[i].resetForce();
	//		for (int j = 0; j < N; ++j) {
	//			if (i != j) {
	//				particles[i].addForce(particles[j]);
	//			}
	//		}
	//	}
	//	//loop again to update the time stamp here
	//	for (auto& particle : particles) {
	//		particle.update(TIMESTAMP);
	//	}
	//	//std::copy(std::begin(particles), std::end(particles),
	//	//	std::istream_iterator<Particle>(std::cout));
	//}
}

int main()
{
	int thread_count = DEFAULT_NUMBER_OF_THREADS;
	double gravity = GRAVITY;
	size_t particle_count = DEFAULT_PARTICLE_COUNT;
	double total_time_steps = DEFAULT_TOTAL_TIME_STEPS;
	double time_step = TIME_STEP;
	size_t universe_size_x = UNIVERSE_SIZE_X;
	size_t universe_size_y = UNIVERSE_SIZE_Y;

	particle_count = 500;
	total_time_steps = 50;

	if (total_time_steps > 0.0 & particle_count > 0 && universe_size_x > 0 && universe_size_y > 0) {
		
		// Print calculation info
		std::cout << "= Parallel N-Body simulation serially and with Thread Building Blocks =" << std::endl;
		std::cout << "Number of threads: " << thread_count << std::endl;
		std::cout << "Total time steps: " << total_time_steps << std::endl;
		std::cout << "Time step: " << time_step << std::endl;
		std::cout << "Particle count: " << particle_count << std::endl << std::endl;
		std::cout << "Universe Size: " << universe_size_x << " x " << universe_size_y << std::endl << std::endl;

		// Advance Game Of Life
		tbb::tick_count before, after; // Execution timers

		// Initialize particle container
		std::vector<Particle> particles;

		// Put random live cells
		ParticleHandler::allocate_random_particles(particle_count, particles, universe_size_x, universe_size_y);

		// Show Init vector universe
		if (VERBOSE) {
			std::cout << "Init Universe" << std::endl;
		}

		// Simulate
		for (double current_Time_step = 0.0; current_Time_step < total_time_steps; current_Time_step += time_step) {
			for (size_t i = 0; i < particle_count; ++i) {
				particles[i].reset_forces();
				for (size_t j = 0; j < particle_count; ++j) {
					if (j != i)
						particles[i].add_force(particles[j]);
				}
			}

			//loop again to update the time stamp here
			for (auto& particle : particles) {
				particle.advance(time_step);
			}
		}

		// Copy the particle universes
		std::vector<Particle> particles_serial(particles);
		tbb::concurrent_vector<Particle, tbb::cache_aligned_allocator<Particle>> particles_tbb(ParticleHandler::to_concurrent_vector(particles));

		// Serial execution
		before = tbb::tick_count::now();
		//simulate_serial(particles_serial, universe_size_x, universe_size_y, total_time_steps); // Advance Simulation serially
		after = tbb::tick_count::now();
		std::cout << std::endl << "Serial execution: " << 1000 * (after - before).seconds() << " ms" << std::endl;

		// Thread Building Blocks		
		before = tbb::tick_count::now();
		//simulate_tbb(particles_tbb, universe_size_x, universe_size_y, total_time_steps); // Advance Simulation with TBB
		after = tbb::tick_count::now();
		std::cout << std::endl << "Thread Building Blocks execution: " << 1000 * (after - before).seconds() << " ms" << std::endl;

		// Show universe
		if (VERBOSE) {

			std::cout << "Final Universe Serial" << std::endl;
			//grid_modifier.debug_show_universe(universe_serial, universe_size_x, universe_size_y);

			std::cout << "Final Universe TBB" << std::endl;
			//grid_modifier.debug_show_universe(UniverseModifier::to_vector(universe_tbb), universe_size_x, universe_size_y);
		}

		if (SAVE_PNG) {
			ParticleHandler::universe_to_png(particles, universe_size_x, universe_size_y, "final_serial_universe.png");
			//grid_modifier.universe_to_png(grid_modifier.to_vector(universe_tbb), universe_size_x, universe_size_y, "final_tbb_universe.png");*/
		}

		system("pause");
	}

}


#include "Settings.h"
#include <iostream>
#include <vector>
#include <tbb/tick_count.h>
#include "Particle.h"
#include <tbb/cache_aligned_allocator.h>
#include <tbb/concurrent_vector.h>
#include "ParticleHandler.h"
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <cassert>
#include "QuadParticleTree.h"

// Advance the simulation using Thread Bulding Blocks parallelization
void simulate_tbb(tbb::concurrent_vector<Particle>& particles, double total_time_steps, double time_step, size_t particle_count,
	size_t universe_size_x, size_t universe_size_y) {

	// Do Simulate
	int png_step_counter = 0;
	for (double current_time_step = 0.0; current_time_step < total_time_steps; current_time_step += time_step) {

		parallel_for(tbb::blocked_range<size_t>(0, particle_count), // Get the range for this thread
			[&](const tbb::blocked_range<size_t>& r) {
				Particle current_particle = particles[r.begin()]; // Thread local variable
				for (size_t i = r.begin(); i != r.end(); ++i) {
					current_particle = particles[i]; // Store the current particle locally
					for (size_t j = i + 1; j < particle_count; ++j) { // Calculate pairs of accelerations
						current_particle.add_acceleration_pairwise(particles[j]);					
					}
					particles[i] = current_particle; // Store data back into the shared memory
				}
			}
		); // Implicit barrier for all the points of the simulation

		// Now that all the new accelerations were calculated, advance the particles in time
		parallel_for(tbb::blocked_range<size_t>(0, particle_count), // Get range for this thread
			[&](const tbb::blocked_range<size_t>& r) {  
				for (size_t index = r.begin(); index != r.end(); ++index) { // Using index range
					particles[index].advance(time_step);
				}
			}
		);

		++png_step_counter;
		if (SAVE_INTERMEDIATE_PNG_STEPS && SAVE_PNG && png_step_counter >= SAVE_PNG_EVERY) { // Save the intermediate step as png
			png_step_counter = 0;

			std::string file_name = "universe_tbb_timestep_" + std::to_string(current_time_step) + ".png";

			ParticleHandler::universe_to_png(ParticleHandler::to_vector(particles), universe_size_x, universe_size_y, file_name.c_str());
		}
	}
}

void simulate_serial_barnes_hut(std::vector<Particle>& particles, double total_time_steps, double time_step, size_t particle_count,
 	size_t universe_size_x, size_t universe_size_y) {

	int png_step_counter = 0;
	QuadParticleTree* quad_tree;
	std::vector<TreeParticle*> results_quad_particle_tree_results;

	// Allocate tree
	quad_tree = ParticleHandler::to_quad_tree(particles, universe_size_x, universe_size_y);

	//quad_particle_tree->get_points_inside_box(qmin_particle, quad_tree->origin, results_quad_particle_tree);
	for (double current_time_step = 0.0; current_time_step < total_time_steps; current_time_step += time_step) {

		std::vector<TreeParticle*> results_quad_particle_tree;
		// Get all resulting particles
		quad_tree->get_all_particles(universe_size_x, universe_size_y, results_quad_particle_tree);

		for (TreeParticle* current_particle : results_quad_particle_tree) {
			quad_tree->apply_acceleration(current_particle->get_particle());
		}

		// Apply acceleration for all the N particles by traversing the quad tree
//		for (size_t i = 0; i < results_quad_particle_tree.size(); i++) {
//			Particle current_particle = results_quad_particle_tree[i]->get_particle();
//			Particle center_of_mass_particle = quad_tree->get_center_of_mass_particle(current_particle);
//			current_particle.add_acceleration(center_of_mass_particle);
//			results_quad_particle_tree[i]->set_particle(current_particle);
//			results_quad_particle_tree_instantiated.push_back(current_particle);
//		}		
		
		//std::vector<TreeParticle*> new_results_quad_particle_tree;
		//quad_tree->get_all_particles(universe_size_x, universe_size_y, new_results_quad_particle_tree);

		// Advance particles
		for (TreeParticle* current_tree_particle : results_quad_particle_tree)
			current_tree_particle->advance(time_step);


		if (current_time_step + time_step >= total_time_steps)
			particles = ParticleHandler::to_vector(results_quad_particle_tree, universe_size_x, universe_size_y);

		// TODO: delete tree
		delete quad_tree;
		// Re-Allocate tree with the resulting points
		quad_tree = ParticleHandler::to_quad_tree(results_quad_particle_tree, universe_size_x, universe_size_y);

		++png_step_counter;
		if (SAVE_INTERMEDIATE_PNG_STEPS && SAVE_PNG && png_step_counter >= SAVE_PNG_EVERY) {

			png_step_counter = 0;
			std::string file_name = "universe_serial_barnes_hut_timestep_" + std::to_string(current_time_step) + ".png";

			// TODO: fix the ability to print universe
			ParticleHandler::universe_to_png(particles, universe_size_x, universe_size_y, file_name.c_str());
		}
	}
}

// Advance the simulation using serial execution
void simulate_serial(std::vector<Particle>& particles, double total_time_steps, double time_step, size_t particle_count,
	size_t universe_size_x, size_t universe_size_y) {

	// Do simulate
	int png_step_counter = 0;
	for (double current_time_step = 0.0; current_time_step < total_time_steps; current_time_step += time_step) {

		// Calculate all the applied forces as acceleration on every particle
		for (size_t i = 0; i < particle_count; ++i) {
			for (size_t j = 0; j < particle_count; ++j) {
				if (j != i)
					particles[i].add_acceleration(particles[j]); // Gather and apply force for every point combination
			}
		}

		for (Particle& current_particle : particles)
			current_particle.advance(time_step); // Advance the particle posiitions in time

		++png_step_counter;
		if (SAVE_INTERMEDIATE_PNG_STEPS && SAVE_PNG && png_step_counter >= SAVE_PNG_EVERY) { // Save the intermediate step as png
		
			png_step_counter = 0;
			std::string file_name = "universe_serial_timestep_" + std::to_string(current_time_step) + ".png";

			ParticleHandler::universe_to_png(particles, universe_size_x, universe_size_y, file_name.c_str());
		}
	}
}

// Application entry point
int main()
{
	// Get the default simulation values
	int thread_count = DEFAULT_NUMBER_OF_THREADS;
	double gravity = GRAVITATIONAL_CONSTANT;
	size_t particle_count = DEFAULT_PARTICLE_COUNT;
	double total_time_steps = DEFAULT_TOTAL_TIME_STEPS;
	double time_step = TIME_STEP;
	size_t universe_size_x = UNIVERSE_SIZE_X;
	size_t universe_size_y = UNIVERSE_SIZE_Y;

	// TODO: User input data
	particle_count = 40;
	total_time_steps = 10;
	universe_size_x = 800;
	universe_size_y = 800;
	thread_count = 4;

	tbb::task_scheduler_init init(thread_count); // Set the number of threads on the TBB scheduler
	
	if (total_time_steps > 0.0 && particle_count > 0 && universe_size_x > 0 && universe_size_y > 0) {
		
		// Print calculation info
		std::cout << "= Parallel N-Body simulation serially and with Thread Building Blocks =" << std::endl;
		std::cout << "Number of threads: " << thread_count << std::endl;
		std::cout << "Total time steps: " << total_time_steps << std::endl;
		std::cout << "Time step: " << time_step << std::endl;
		std::cout << "Particle count: " << particle_count << std::endl << std::endl;
		std::cout << "Universe Size: " << universe_size_x << " x " << universe_size_y << std::endl << std::endl;

		tbb::tick_count before, after; // Execution timers

		// Initialize particle container
		std::vector<Particle> particles;

		// Put random particles
		ParticleHandler::allocate_random_particles(particle_count, particles, universe_size_x, universe_size_y);

		// TODO: Show Init vector universe
		if (VERBOSE) {
			std::cout << "Init Universe" << std::endl;
		}

		// Simulate

		// Copy the particle universes into the serial and parallel execution containers
		std::vector<Particle> particles_serial(particles);
		tbb::concurrent_vector<Particle, tbb::cache_aligned_allocator<Particle>> particles_tbb(ParticleHandler::to_concurrent_vector(particles));
		std::vector<Particle> particles_serial_barnes_hut(particles);

		// Benchmark the Serial execution
		std::cout << std::endl << "Serial execution... ";
		before = tbb::tick_count::now();
		simulate_serial(particles_serial, total_time_steps, time_step, particle_count, universe_size_x, universe_size_y); // Advance Simulation serially
		after = tbb::tick_count::now();
		std::cout << 1000 * (after - before).seconds() << " ms" << std::endl;

		// Serial execution
		std::cout << std::endl << "Serial execution (Barnes-Hut)... ";
		before = tbb::tick_count::now();
		simulate_serial_barnes_hut(particles_serial_barnes_hut, total_time_steps, time_step, particle_count, universe_size_x, universe_size_y); // Advance Simulation serially
		after = tbb::tick_count::now();
		std::cout << 1000 * (after - before).seconds() << " ms" << std::endl;
		
		// Benchmark the Thread Building Blocks execution
		std::cout << std::endl << "Thread Building Blocks execution... ";
		before = tbb::tick_count::now();
		simulate_tbb(particles_tbb, total_time_steps, time_step, particle_count, universe_size_x, universe_size_y); // Advance Simulation with TBB
		after = tbb::tick_count::now();
		std::cout << 1000 * (after - before).seconds() << " ms" << std::endl;

		// Assert the equality and validity of the results
		assert(ParticleHandler::are_equal(particles_serial, ParticleHandler::to_vector(particles_tbb)) == true); // compare serial with parallel
		assert(ParticleHandler::are_equal(particles, particles_serial) == false); // compare serial with init
		assert(ParticleHandler::are_equal(particles, ParticleHandler::to_vector(particles_tbb)) == false); // compare parallel with init

		// TODO: Show universe in the console
		if (VERBOSE) {

			std::cout << "Final Universe Serial" << std::endl;
			//grid_modifier.debug_show_universe(universe_serial, universe_size_x, universe_size_y);

			std::cout << "Final Universe TBB" << std::endl;
			//grid_modifier.debug_show_universe(UniverseModifier::to_vector(universe_tbb), universe_size_x, universe_size_y);
		}

		if (SAVE_PNG) { // Save final universes to png
			ParticleHandler::universe_to_png(particles, universe_size_x, universe_size_y, "init_universe.png");
			ParticleHandler::universe_to_png(particles_serial, universe_size_x, universe_size_y, "final_serial_universe.png");
			ParticleHandler::universe_to_png(particles_serial_barnes_hut, universe_size_x, universe_size_y, "final_serial_universe_barnes_hut.png");
			ParticleHandler::universe_to_png(ParticleHandler::to_vector(particles_tbb), universe_size_x, universe_size_y, "final_tbb_universe.png");
		}

		system("pause");
	}

}


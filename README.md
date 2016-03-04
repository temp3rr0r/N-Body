# Parallel N-Body Simulator#

Simulations of dynamical systems of particles are often used in physics to predict behavior of planets, stars, galaxies, gas particles... The interaction between particles is described by physically sound equations and ”integrated” in time to predict the outcome of the simulation. Write a parallel version of the N-body algorithm in C++ using Intel’s Threading Building Blocks library for the parallel code. A straightforward approach is to use the ”naive” approach described above by using the three nested loops (over j, i and t, from inner to outer).

### Documentation
Documentation comparing the speed-up between the serial and parallel versions: http://1drv.ms/1TvdwIo

###References
1. *Intel Developer zone - n-bodies: a parallel TBB solution: computing accelerations? or forces?* - https://software.intel.com/enus/blogs/2009/09/22/n-bodies-a-parallel-tbb-solution-computingaccelerations-or-forces
2. *Wikipedia: Fast Inverse Square Root* - https://en.wikipedia.org/wiki/Fastinversesquareroot
3. *Wikipedia: Barnes-Hut simulation* - https://en.wikipedia.org/wiki/BarnesHutsimulation
4. *Tutorial: Scalable Memory Allocator* - https://www.threadingbuildingblocks.org/tutorial-intel-tbb-scalablememory-allocator
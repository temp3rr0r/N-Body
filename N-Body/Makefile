all:
	$(CXX) N-Body.cpp -O3 -o nbody -std=c++11 -ltbb -lpng -lz
benchmark:
	$(CXX) Benchmark_N-Body.cpp -O3 -o benchmarknbody -std=c++11 -ltbb -lpng -lz
run:
	./nbody
clean:
	rm *.png
	rm *.csv
	rm nbody
	rm benchmarknbody

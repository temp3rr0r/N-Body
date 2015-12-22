#pragma once
#include <cstdint>

// Default N-Body simulation settings

static const bool VERBOSE = false;
static const bool SAVE_PNG = true;

static const bool SAVE_INTERMEDIATE_PNG_STEPS = false;
static const int SAVE_PNG_EVERY = 500;

static const int DEFAULT_NUMBER_OF_THREADS = 4;

static const int DEFAULT_PARTICLE_COUNT = 10;
static const float DEFAULT_TOTAL_TIME_STEPS = 10.0f;
static const float TIME_STEP = 0.01f;
static const float MIN_DISTANCE = 10.0f;

static const float GRAVITATIONAL_CONSTANT = 6.673e-11f;

static const float MAX_RANDOM = 1.0f;

static const float MAX_MASS = 1.0f;

static const float THETA = 0.5f;
static const uint8_t MAX_TREE_DEPTH = 100;

static const int UNIVERSE_SIZE_X = 800;
static const int UNIVERSE_SIZE_Y = 800;




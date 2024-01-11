#ifndef CONSTANTS_H
#define CONSTANTS_H

// Constants for grid dimensions
const int NSTEPS = 10000;        // Number of simulation steps
const int ITERATIONS = 0;        // Number of iterations
const double DX = 0.01;          // Grid spacing in the x-direction
const double DY = 0.01;          // Grid spacing in the y-direction

// Constants for the heat equation simulation
const double DIFFUSION_CONSTANT = 0.5;  // Diffusion constant determines the rate of heat transfer in the material

// Constants for image extraction
const int IMAGE_OUTPUT_INTERVAL = 1000;   // Image output interval determines how often the simulation results are saved as images
const int RESTART_OUTPUT_INTERVAL = 2000; // Restart output interval defines the frequency of creating restart checkpoints for the simulation

// MPI tags for communication
const int TAG_WRITE = 22;   // Tag for MPI communication during writing
const int TAG_READ = 44;    // Tag for MPI communication during reading

// File configs
const int FILENAME_LENGTH = 64;                // Maximum length for filename
const char CHECKPOINT[] = "HEAT.dat";          // Default filename

// Default field dimensions
const int DEFAULT_ROWS = 2000;
const int DEFAULT_COLS = 4000;

// Temperature values
const double INSIDE_DISC_TEMPERATURE = 65.0;   // Temperature inside the disc
const double OUTSIDE_DISC_TEMPERATURE = 5.0;   // Temperature outside the disc

// MPI communication tags for halo exchange
const int ROW_TAG_UP = 11;        // MPI tag for sending/receiving rows to/from the upper neighbor
const int ROW_TAG_DOWN = 12;      // MPI tag for sending/receiving rows to/from the lower neighbor
const int COLUMN_TAG_LEFT = 13;   // MPI tag for sending/receiving columns to/from the left neighbor
const int COLUMN_TAG_RIGHT = 14;  // MPI tag for sending/receiving columns to/from the right neighbor

#endif // CONSTANTS_H

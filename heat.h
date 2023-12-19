#ifndef HEAT_H
#define HEAT_H

// Data structure for representing the temperature field
struct Field {
    int nx;            // Local dimensions of the field in the x-direction
    int ny;            // Local dimensions of the field in the y-direction
    int nx_full;       // Global dimensions of the field in the x-direction
    int ny_full;       // Global dimensions of the field in the y-direction
    double dx;         // Grid spacing in the x-direction
    double dy;         // Grid spacing in the y-direction
    double *data;      // Array containing temperature field data
};

// Structure holding information for parallelizing tasks using MPI
struct ParallelData {
    int size;                       // Total number of MPI tasks
    int rank;                       // Rank of the current MPI task
    int nup, ndown, nleft, nright;  // Ranks of neighboring MPI tasks (up, down, left, right)
    MPI_Comm comm;                  // MPI Cartesian communicator
    MPI_Request requests[8];        // Array of requests for non-blocking communication
    MPI_Datatype rowtype;           // MPI Datatype for communication of rows
    MPI_Datatype columntype;        // MPI Datatype for communication of columns
    MPI_Datatype subarraytype;      // MPI Datatype for communication in text I/O
    MPI_Datatype restarttype;       // MPI Datatype for communication in restart I/O
    MPI_Datatype filetype;          // MPI Datatype for file view in restart I/O
};

// Inline function for indexing 2D arrays
static inline int idx(int i, int j, int width)
{
    return i * width + j;
}

// Function prototypes

/**
 * @brief Allocate memory for a two-dimensional array.
 *
 * This function allocates memory for a two-dimensional array with dimensions nx and ny.
 *
 * @param nx Size of the first dimension.
 * @param ny Size of the second dimension.
 * @return Pointer to the allocated memory.
 */
double *allocate_2d_array(int nx, int ny);

/**
 * @brief Set dimensions of the field based on the parallel configuration.
 *
 * This function calculates the local dimensions of the field based on the overall dimensions (nx, ny) and the parallel configuration.
 *
 * @param temperature Pointer to the field structure.
 * @param nx Size of the first dimension of the global field.
 * @param ny Size of the second dimension of the global field.
 * @param parallel Pointer to the parallel data structure.
 */
void init_field_properties(Field *temperature, int nx, int ny, ParallelData *parallel);

/**
 * @brief Initialize parallel data and set up communication.
 *
 * This function initializes parallel data and sets up communication using MPI. It creates
 * a Cartesian communicator, determines local domain sizes, checks if the grid can be evenly
 * divided among processors, and creates data types for halo exchange, subarray needed in
 * text I/O, and restart I/O.
 *
 * @param parallel Pointer to the ParallelData struct to store parallel information.
 * @param nx Total size of the domain in the x-direction.
 * @param ny Total size of the domain in the y-direction.
 */
void init_parallel_data(ParallelData *parallel, int nx, int ny);

/**
 * @brief Initialize the heat equation solver.
 *
 * This function initializes the heat equation solver by processing command line arguments,
 * setting up initial field conditions, and checking for restart or input file options. It
 * also handles domain decomposition using MPI, reads field data from a file, and sets up
 * the initial conditions based on the provided parameters.
 *
 * @param argc Number of command line arguments.
 * @param argv Array of command line argument strings.
 * @param current Pointer to the current field data.
 * @param previous Pointer to the previous field data.
 * @param nsteps Pointer to the number of time steps.
 * @param parallel Pointer to the ParallelData struct for MPI communication.
 * @param iter0 Pointer to the iteration number for restarts.
 */
void init_simulation(int argc, char *argv[], Field *temperature1, Field *temperature2, int *nsteps, ParallelData *parallel, int *iter0);

/**
 * @brief Generate the initial temperature field.
 *
 * The pattern is a disc with a radius of nx_full / 6 in the center of the grid.
 * Boundary conditions are different constant temperatures outside the grid.
 *
 * @param temperature Pointer to the field structure representing the temperature.
 * @param parallel Pointer to the parallel data structure.
 */
void init_heat_field(Field *temperature, ParallelData *parallel);

/**
 * @brief Exchange the boundary values between neighboring processes.
 *
 * This function uses non-blocking communication to exchange boundary values
 * with neighboring processes in a 2D Cartesian communicator.
 *
 * @param temperature Pointer to the field structure representing the temperature.
 * @param parallel Pointer to the parallel data structure.
 */
void start_halo_exchange(Field *temperature, ParallelData *parallel);

/**
 * @brief Complete the non-blocking communication by waiting for completion.
 *
 * This function waits for the completion of all non-blocking communication requests
 * initiated during the exchange of boundary values between neighboring processes.
 *
 * @param parallel Pointer to the parallel data structure.
 */
void complete_halo_exchange(ParallelData *parallel);

/**
 * @brief Update the temperature values using a five-point stencil.
 *
 * This function applies the five-point stencil to evolve the interior temperature
 * field at the next time step based on the previous temperature field.
 *
 * @param curr Pointer to the field structure representing the current temperature.
 * @param prev Pointer to the field structure representing the previous temperature.
 * @param a Coefficient for the heat equation.
 * @param dt Time step size.
 */
void update_interior_temperature(Field *curr, Field *prev, double a, double dt);

/**
 * @brief Update the temperature values using a five-point stencil for border-dependent regions.
 *
 * This function selectively updates the border-dependent regions of the temperature field
 * at the next time step using the five-point stencil, considering fixed boundary conditions.
 *
 * @param curr Pointer to the field structure representing the current temperature.
 * @param prev Pointer to the field structure representing the previous temperature.
 * @param a Coefficient for the heat equation.
 * @param dt Time step size.
 */
void update_boundary_temperature(Field *curr, Field *prev, double a, double dt);

/**
 * @brief Output routine that prints a picture of the temperature distribution.
 *
 * This function outputs the temperature distribution to a PNG file, considering
 * the parallelization and the presence of ghost layers.
 *
 * @param temperature Pointer to the field structure representing the temperature.
 * @param iter Iteration number for filename distinction.
 * @param parallel Pointer to the parallel data structure.
 */
void write_field_to_file(Field *temperature, int iter, ParallelData *parallel);

/**
 * @brief Read the initial temperature distribution from a file and initialize
 * the temperature fields temperature1 and temperature2 to the same initial state.
 *
 * This function reads the temperature distribution from a file, sets up parallelization,
 * allocates arrays (including ghost layers), and distributes the data among processes.
 * It also sets the boundary values and copies the initial state to another field.
 *
 * @param temperature1 Pointer to the first field structure representing the temperature.
 * @param temperature2 Pointer to the second field structure representing the temperature.
 * @param filename Name of the file containing the initial temperature distribution.
 * @param parallel Pointer to the parallel data structure.
 */
void read_field_from_file(Field *temperature1, Field *temperature2, char *filename, ParallelData *parallel);

/**
 * @brief Write a restart checkpoint that contains field dimensions, current
 * iteration number, and temperature field.
 *
 * This function opens the file, writes field dimensions, current iteration number,
 * and the temperature field to a restart checkpoint file using MPI I/O.
 *
 * @param temperature Pointer to the field structure representing the temperature.
 * @param parallel Pointer to the parallel data structure.
 * @param iter Current iteration number.
 */
void write_restart_data(Field *temperature, ParallelData *parallel, int iter);

/**
 * @brief Read a restart checkpoint that contains field dimensions, current
 * iteration number, and temperature field.
 *
 * This function opens the restart checkpoint file, reads field dimensions, 
 * current iteration number, and temperature field using MPI I/O. It then sets 
 * the correct dimensions to MPI metadata, sets local dimensions, and allocates 
 * memory for the temperature field.
 *
 * @param temperature Pointer to the field structure representing the temperature.
 * @param parallel Pointer to the parallel data structure.
 * @param iter Pointer to the variable holding the current iteration number.
 */
void read_restart_data(Field *temperature, ParallelData *parallel, int *iter);

/**
 * @brief Copy data from temperature1 into temperature2.
 *
 * This function copies the temperature field data from one field structure
 * (temperature1) into another (temperature2).
 *
 * @param temperature1 Pointer to the source field structure.
 * @param temperature2 Pointer to the destination field structure.
 */
void copy_field_data(Field *temperature1, Field *temperature2);

/**
 * @brief Swap the data of two field structures.
 *
 * This function swaps the data pointers of two field structures, 
 * allowing for efficient exchange of data without copying.
 *
 * @param temperature1 Pointer to the first field structure.
 * @param temperature2 Pointer to the second field structure.
 */
void swap_field_data(Field *temperature1, Field *temperature2);

/**
 * @brief Allocate memory for a temperature field and initialize it to zero.
 *
 * This function allocates memory for a temperature field, including ghost layers,
 * and initializes the field to zero.
 *
 * @param temperature Pointer to the field structure.
 */
void allocate_field(Field *temperature);

/**
 * @brief Deallocate the 2D arrays and MPI datatypes of temperature fields.
 *
 * This function deallocates the memory for the 2D arrays of temperature fields,
 * as well as the MPI datatypes used for communication.
 *
 * @param temperature1 Pointer to the first field structure.
 * @param temperature2 Pointer to the second field structure.
 * @param parallel Pointer to the parallel data structure.
 */
void cleanup_resources(Field *temperature1, Field *temperature2, ParallelData *parallel);

#endif  // HEAT_H

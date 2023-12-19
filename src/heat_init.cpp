// Standard C++ and MPI headers for parallel programming
#include <iostream>       // Standard C++ I/O
#include <cstdlib>        // Standard C library functions
#include <cstring>        // String manipulation functions
#include <cassert>        // Assertion macros
#include <cstdio>         // C Standard Input and Output Library
#include <cmath>          // Math functions
#include <unistd.h>       // POSIX operating system API

#include <mpi.h>          // MPI for parallel programming

// Custom libraries for PNG writing and custom heat-related functionality
#include "pngsaver.h"    // Custom library for PNG writing
#include "heat.h"         // Custom library for heat-related functionality
#include "constants.h"    // Constants header file

double* allocate_2d_array(int nx, int ny) {
    double* array = new double[nx * ny];
    return array;
}

void init_field_properties(Field* temperature, int nx, int ny, ParallelData* parallel) {
    int nx_local, ny_local;
    int dims[2], coords[2], periods[2];

    // Get Cartesian coordinates and dimensions
    MPI_Cart_get(parallel->comm, 2, dims, periods, coords);
    nx_local = nx / dims[0];
    ny_local = ny / dims[1];

    // Set field properties
    temperature->dx = DX;
    temperature->dy = DY;
    temperature->nx = nx_local;
    temperature->ny = ny_local;
    temperature->nx_full = nx;
    temperature->ny_full = ny;
}

void init_parallel_data(ParallelData* parallel, int nx, int ny) {
    int nx_local;
    int ny_local;
    int world_size;
    int dims[2] = {0, 0};
    int periods[2] = {0, 0};

    // Get the total number of MPI processes
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Create Cartesian grid dimensions
    MPI_Dims_create(world_size, 2, dims);

    // Calculate local domain sizes
    nx_local = nx / dims[0];
    ny_local = ny / dims[1];

    // Check if the grid can be evenly divided among processors
    if (nx_local * dims[0] != nx) {
        std::cerr << "Cannot divide grid evenly among processors in the x-direction. " << nx_local << " x " << dims[0] << " != " << nx << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -2);
    }
    if (ny_local * dims[1] != ny) {
        std::cerr << "Cannot divide grid evenly among processors in the y-direction. " << ny_local << " x " << dims[1] << " != " << ny << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -2);
    }

    // Create Cartesian communicator
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &parallel->comm);
    MPI_Cart_shift(parallel->comm, 0, 1, &parallel->nup, &parallel->ndown);
    MPI_Cart_shift(parallel->comm, 1, 1, &parallel->nleft, &parallel->nright);

    // Get rank and size
    MPI_Comm_size(parallel->comm, &parallel->size);
    MPI_Comm_rank(parallel->comm, &parallel->rank);

    // Print information (only from rank 0)
    if (parallel->rank == 0) {
        std::cout << "Using domain decomposition " << dims[0] << " x " << dims[1] << std::endl;
        std::cout << "Local domain size " << nx_local << " x " << ny_local << std::endl;
    }

    // Create data types for halo exchange
    MPI_Type_vector(nx_local + 2, 1, ny_local + 2, MPI_DOUBLE, &parallel->columntype);
    MPI_Type_contiguous(ny_local + 2, MPI_DOUBLE, &parallel->rowtype);
    MPI_Type_commit(&parallel->columntype);
    MPI_Type_commit(&parallel->rowtype);

    // Define sizes, subsizes, and offsets for subarray type
    int sizes[2] = {nx_local + 2, ny_local + 2};
    int subsizes[2] = {nx_local, ny_local};
    int offsets[2] = {1, 1};

    // Adjust sizes and offsets for rank 0
    if (parallel->rank == 0) {
        sizes[0] = nx;
        sizes[1] = ny;
        offsets[0] = 0;
        offsets[1] = 0;
    }

    // Create data type for subarray needed in text I/O
    MPI_Type_create_subarray(2, sizes, subsizes, offsets, MPI_ORDER_C, MPI_DOUBLE, &parallel->subarraytype);
    MPI_Type_commit(&parallel->subarraytype);

    // Create data type for restart I/O
    int coords[2];
    MPI_Cart_coords(parallel->comm, parallel->rank, 2, coords);
    sizes[0] = nx + 2;
    sizes[1] = ny + 2;
    offsets[0] = 1 + coords[0] * nx_local;
    offsets[1] = 1 + coords[1] * ny_local;

    // Adjust sizes and offsets for boundary ranks
    if (coords[0] == 0) {
        offsets[0] -= 1;
        subsizes[0] += 1;
    }
    if (coords[0] == dims[0] - 1) {
        subsizes[0] += 1;
    }
    if (coords[1] == 0) {
        offsets[1] -= 1;
        subsizes[1] += 1;
    }
    if (coords[1] == dims[1] - 1) {
        subsizes[1] += 1;
    }

    // Create data type for restart I/O
    MPI_Type_create_subarray(2, sizes, subsizes, offsets, MPI_ORDER_C, MPI_DOUBLE, &parallel->filetype);
    MPI_Type_commit(&parallel->filetype);

    // Create data type for restart I/O
    sizes[0] = nx_local + 2;
    sizes[1] = ny_local + 2;
    offsets[0] = 1;
    offsets[1] = 1;
    
    // Adjust offsets for boundary ranks
    if (coords[0] == 0) {
        offsets[0] = 0;
    }
    if (coords[1] == 0) {
        offsets[1] = 0;
    }

    // Create data type for restart I/O
    MPI_Type_create_subarray(2, sizes, subsizes, offsets, MPI_ORDER_C, MPI_DOUBLE, &parallel->restarttype);
    MPI_Type_commit(&parallel->restarttype);
}

void init_simulation(int argc, char *argv[], Field *current, Field *previous, int *nsteps, ParallelData *parallel, int *initialIteration) {
    // Default field dimensions
    int rows = DEFAULT_ROWS;
    int cols = DEFAULT_COLS;

    // Name of the optional input file
    char input_file[64];

    // Flag to indicate whether to read a file
    int read_file = 0;

    // Default values for time steps and iteration
    *nsteps = NSTEPS;
    *initialIteration = ITERATIONS;

    // Check if checkpoint exists
    if (access(CHECKPOINT, F_OK) == 0) {
        // Read from a restart checkpoint
        read_restart_data(current, parallel, initialIteration);

        // Set field dimensions and allocate memory for the previous field
        init_field_properties(previous, current->nx_full, current->ny_full, parallel);
        allocate_field(previous);

        // Print restart information (only from rank 0)
        if (parallel->rank == 0)
            std::cout << "Restarting from an earlier checkpoint saved at iteration " << *initialIteration << "." << std::endl;

        // Copy data from the current to the previous field
        copy_field_data(current, previous);
    } else if (read_file) {
        // Read field data from the specified input file
        read_field_from_file(current, previous, input_file, parallel);
    } else {
        // Set up parallel domain decomposition
        init_parallel_data(parallel, rows, cols);

        // Set field dimensions for both current and previous fields
        init_field_properties(current, rows, cols, parallel);
        init_field_properties(previous, rows, cols, parallel);

        // Generate initial field conditions for the current field
        init_heat_field(current, parallel);

        // Allocate memory for the previous field and copy data from the current field
        allocate_field(previous);
        copy_field_data(current, previous);
    }
}

void init_heat_field(Field *temperature, ParallelData *parallel) {
    int i, j, ind, width;
    double radius;
    int dx, dy;
    int dims[2], coords[2], periods[2];

    // Allocate the temperature array, including the ghost layers
    temperature->data = allocate_2d_array(temperature->nx + 2, temperature->ny + 2);

    MPI_Cart_get(parallel->comm, 2, dims, periods, coords);

    // Radius of the source disc
    radius = temperature->nx_full / 6.0;

    width = temperature->ny + 2;

    // Iterate over the grid points to set the initial temperature
    for (i = 0; i < temperature->nx + 2; i++) {
        for (j = 0; j < temperature->ny + 2; j++) {
            // Distance of point (i, j) from the origin
            dx = i + coords[0] * temperature->nx - temperature->nx_full / 2 + 1;
            dy = j + coords[1] * temperature->ny - temperature->ny_full / 2 + 1;
            ind = idx(i, j, width);

            // Set temperature based on whether the point is inside the disc
            if (dx * dx + dy * dy < radius * radius) {
                temperature->data[ind] = INSIDE_DISC_TEMPERATURE;  // Inside the disc
            } else {
                temperature->data[ind] = OUTSIDE_DISC_TEMPERATURE; // Outside the disc
            }
        }
    }

    // Boundary conditions
    // Left boundary
    if (coords[1] == 0) {
        for (i = 0; i < temperature->nx + 2; i++) {
            ind = idx(i, 0, width);
            temperature->data[ind] = 20.0;
        }
    }
    // Right boundary
    if (coords[1] == dims[1] - 1) {
        for (i = 0; i < temperature->nx + 2; i++) {
            ind = idx(i, temperature->ny + 1, width);
            temperature->data[ind] = 70.0;
        }
    }
    // Upper boundary
    if (coords[0] == 0) {
        for (j = 0; j < temperature->ny + 2; j++) {
            ind = idx(0, j, width);
            temperature->data[ind] = 85.0;
        }
    }
    // Lower boundary
    if (coords[0] == dims[0] - 1) {
        for (j = 0; j < temperature->ny + 2; j++) {
            ind = idx(temperature->nx + 1, j, width);
            temperature->data[ind] = 5.0;
        }
    }
}

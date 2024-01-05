#include <iostream>
#include <vector>
#include <mpi.h>

#include "heat.h"
#include "constants.h"

/**
 * @brief Main function for the heat equation solver.
 *
 * This program simulates the heat equation using an explicit finite-difference method.
 * It initializes the simulation, performs time evolution, outputs results at specified
 * intervals, and writes restart checkpoints for easy restarting. The program uses MPI
 * for parallelization to handle domain decomposition and communication between processes.
 *
 * @param argc Number of command line arguments.
 * @param argv Array of command line argument strings.
 *
 * @return Exit code (0 for success, non-zero for errors).
 */
int main(int argc, char **argv)
{
    int dims[2], coords[2], periods[2];
    // Current and previous temperature fields
    Field currentField, previousField;

    // Time step, number of time steps
    double timeStep;
    int numTimeSteps;

    // Parallelization info
    ParallelData parallelInfo;

    // Iteration counters
    int iteration, initialIteration;

    // Delta x and y squared
    double deltaX2 = DX * DX;
    double deltaY2 = DY * DY;

    // Time stamps
    double startTimeStamp;
    double endTimeStamp;

    MPI_Init(&argc, &argv);

    init_simulation(argc, argv, &currentField, &previousField, &numTimeSteps, &parallelInfo, &initialIteration);

    //Brauchen coords um Boundary richtig zuaktualisieren
    MPI_Cart_get(parallel->comm, 2, dims, periods, coords);
    
    // Output the initial field
    write_field_to_file(&currentField, initialIteration, &parallelInfo);
    initialIteration++;

    // Calculate the largest stable time step
    timeStep = (deltaY2 * deltaX2) / (2 * DIFFUSION_CONSTANT * (deltaX2 + deltaY2));

    // Get the start time stamp
    startTimeStamp = MPI_Wtime();

    // Time evolve
    for (iteration = initialIteration; iteration < initialIteration + numTimeSteps; iteration++) {
        // ToDo
        // Use your implemented functions in a correct order here.

        // Output field at specified intervals
        if (iteration % IMAGE_OUTPUT_INTERVAL == 0) {
            write_field_to_file(&currentField, iteration, &parallelInfo);
        }

        // Write a checkpoint for easy restarting
        if (iteration % RESTART_OUTPUT_INTERVAL == 0) {
            write_restart_data(&currentField, &parallelInfo, iteration);
        }

        // Swap current and previous fields
        swap_field_data(&currentField, &previousField);
    }

    endTimeStamp = MPI_Wtime();

    // Determine the CPU time used for the iteration
    if (parallelInfo.rank == 0) {
        std::cout << "Iteration took " << (endTimeStamp - startTimeStamp) << " seconds." << std::endl;
        std::cout << "Reference value at 7,3: " << previousField.data[idx(7, 3, currentField.ny + 2)] << std::endl;
    }

    // Output the final field
    write_field_to_file(&currentField, iteration, &parallelInfo);

    cleanup_resources(&currentField, &previousField, &parallelInfo);
    MPI_Finalize();

    return 0;
}

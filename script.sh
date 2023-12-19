#!/bin/bash
#SBATCH -J heat_mpi_job            # Job name
#SBATCH -o heat_mpi_output_%j.out  # Output file name
#SBATCH -e heat_mpi_error_%j.err   # Error file name
#SBATCH -p s_hadoop                # Specified partition
#SBATCH --time=00:20:00            # Max wall time

# Purge all currently loaded modules to start with a clean environment
module purge

# Load necessary modules
module load compiler/gcc/11.2.0
module load mpi/openmpi/4.1.2-gcc-10.2.0

# Compile the MPI program
mpic++ -O3 -Wall -c main.cpp pngsaver.cpp heat_init.cpp heat_update.cpp heat_io.cpp

# Link the compiled objects into an executable

# -O3: Optimization level 3
# -Wall: Enable most warning messages
# -o: Output file name
# -lpng: Link against the libpng library
# -lm: Link against the math library
mpic++ -O3 -Wall -o heat_mpi main.o pngsaver.o heat_init.o heat_update.o heat_io.o -lpng -lm

# Run the MPI program using mpirun
# -np 8: Number of MPI tasks
# --oversubscribe: Allow more processes than available cores
mpirun --oversubscribe -n 4 heat_mpi

# Clean up object files after the run
rm -f *.o

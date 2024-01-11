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

mkdir build
# Compile the MPI program
mpic++ -O3 -Wall -c src/main.cpp -o build/main.o
mpic++ -O3 -Wall -c src/pngsaver.cpp -o build/pngsaver.o 
mpic++ -O3 -Wall -c src/heat_init.cpp -o build/heat_init.o
mpic++ -O3 -Wall -c src/heat_update.cpp -o build/heat_update.o
mpic++ -O3 -Wall -c src/heat_io.cpp -o build/heat_io.o

# Link the compiled objects into an executable

# -O3: Optimization level 3
# -Wall: Enable most warning messages
# -o: Output file name
# -lpng: Link against the libpng library
# -lm: Link against the math library
mpic++ -O3 -Wall -o build/heat_mpi build/main.o build/pngsaver.o build/heat_init.o build/heat_update.o build/heat_io.o -lpng -lm

# Run the MPI program using mpirun
# -np 8: Number of MPI tasks
# --oversubscribe: Allow more processes than available cores
mpirun --oversubscribe -n 16 build/heat_mpi

# Clean up object files after the run
rm -f build/*.o

#include <mpi.h>
#include <iostream>
#include <cstring>
#include <cassert>
#include "pngsaver.h"
#include "heat.h"
#include "constants.h"

void write_field_to_file(Field *temperature, int iter, ParallelData *parallel) {
    char filename[FILENAME_LENGTH];

    // The actual write routine takes only the actual data (without ghost layers),
    // so we need an array for that.
    int height, width;
    double *full_data;

    int coords[2];
    int ix, jy;

    int i, p;

    height = temperature->nx_full;
    width = temperature->ny_full;

    if (parallel->rank == 0) {
        // Copy the inner data
        full_data = allocate_2d_array(height, width);
        for (i = 0; i < temperature->nx; i++)
            memcpy(&full_data[idx(i, 0, width)], &temperature->data[idx(i + 1, 1, temperature->ny + 2)], temperature->ny * sizeof(double));

        // Receive data from other ranks
        for (p = 1; p < parallel->size; p++) {
            MPI_Cart_coords(parallel->comm, p, 2, coords);
            ix = coords[0] * temperature->nx;
            jy = coords[1] * temperature->ny;
            MPI_Recv(&full_data[idx(ix, jy, width)], 1, parallel->subarraytype, p, TAG_WRITE, parallel->comm, MPI_STATUS_IGNORE);
        }

        // Store the data to a PNG file
        sprintf(filename, "%s_%04d.png", "heat", iter);
        savePngImage(full_data, height, width, filename, 'c');
        delete[] full_data;
    } else {
        // Send data
        MPI_Ssend(temperature->data, 1, parallel->subarraytype, 0, TAG_WRITE, parallel->comm);
    }
}

void read_field_from_file(Field *temperature1, Field *temperature2, char *filename, ParallelData *parallel) {
    FILE *fp;
    int nx, ny, i, j;
    double *full_data = nullptr; // Initialize full_data to nullptr

    int coords[2];
    int ix, jy, p;

    int count;

    // Open the file for reading
    fp = fopen(filename, "r");
    
    // Read the header
    count = fscanf(fp, "# %d %d \n", &nx, &ny);
    if (count < 2) {
        fprintf(stderr, "Error while reading the input file!\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // Set up parallelization and field dimensions
    init_parallel_data(parallel, nx, ny);
    init_field_properties(temperature1, nx, ny, parallel);
    init_field_properties(temperature2, nx, ny, parallel);

    // Allocate arrays (including ghost layers)
    temperature1->data = allocate_2d_array(temperature1->nx + 2, temperature1->ny + 2);
    temperature2->data = allocate_2d_array(temperature2->nx + 2, temperature2->ny + 2);

    if (parallel->rank == 0) {
        // Full array
        full_data = allocate_2d_array(nx, ny);

        // Read the actual data
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                count = fscanf(fp, "%lf", &full_data[idx(i, j, ny)]);
            }
        }

        // Copy to own local array
        for (i = 0; i < temperature1->nx; i++) {
            memcpy(&temperature1->data[idx(i + 1, 1, temperature1->ny + 2)], &full_data[idx(i, 0, ny)], temperature1->ny * sizeof(double));
        }

        // Send to other processes
        for (p = 1; p < parallel->size; p++) {
            MPI_Cart_coords(parallel->comm, p, 2, coords);
            ix = coords[0] * temperature1->nx;
            jy = coords[1] * temperature1->ny;
            MPI_Send(&full_data[idx(ix, jy, ny)], 1, parallel->subarraytype, p, TAG_READ, parallel->comm);
        }
    } else {
        // Receive data
        MPI_Recv(temperature1->data, 1, parallel->subarraytype, 0, TAG_READ, parallel->comm, MPI_STATUS_IGNORE);
    }

    // Set the boundary values
    for (i = 0; i < temperature1->nx + 1; i++) {
        temperature1->data[idx(i, 0, temperature1->ny + 2)] = temperature1->data[idx(i, 1, temperature1->ny + 2)];
        temperature1->data[idx(i, temperature1->ny + 1, temperature1->ny + 2)] = temperature1->data[idx(i, temperature1->ny, temperature1->ny + 2)];
    }

    for (j = 0; j < temperature1->ny + 2; j++) {
        temperature1->data[idx(0, j, temperature1->ny + 2)] = temperature1->data[idx(1, j, temperature1->ny + 2)];
        temperature1->data[idx(temperature1->nx + 1, j, temperature1->ny + 2)] = temperature1->data[idx(temperature1->nx, j, temperature1->ny + 2)];
    }

    // Copy the initial state to another field
    copy_field_data(temperature1, temperature2);

    if (parallel->rank == 0) {
        delete[] full_data;
    }

    // Close the file
    fclose(fp);
}

void write_restart_data(Field *temperature, ParallelData *parallel, int iter) {
    MPI_File fp;
    MPI_Offset disp;

    // Open the file and write the dimensions
    MPI_File_open(MPI_COMM_WORLD, CHECKPOINT, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);

    if (parallel->rank == 0) {
        MPI_File_write(fp, &temperature->nx_full, 1, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_write(fp, &temperature->ny_full, 1, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_write(fp, &iter, 1, MPI_INT, MPI_STATUS_IGNORE);
    }

    // Set the displacement for temperature data
    disp = 3 * sizeof(int);

    // Set the file view and write temperature data
    MPI_File_set_view(fp, 0, MPI_DOUBLE, parallel->filetype, "native", MPI_INFO_NULL);
    MPI_File_write_at_all(fp, disp, temperature->data, 1, parallel->restarttype, MPI_STATUS_IGNORE);

    // Close the file
    MPI_File_close(&fp);
}

void read_restart_data(Field *temperature, ParallelData *parallel, int *iter) {
    MPI_File fp;
    MPI_Offset disp;

    int nx, ny;

    // Open the file and read grid size and current iteration
    MPI_File_open(MPI_COMM_WORLD, CHECKPOINT, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp);

    MPI_File_read_all(fp, &nx, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_all(fp, &ny, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_all(fp, iter, 1, MPI_INT, MPI_STATUS_IGNORE);

    // Set correct dimensions to MPI metadata
    init_parallel_data(parallel, nx, ny);

    // Set local dimensions and allocate memory for the data
    init_field_properties(temperature, nx, ny, parallel);
    allocate_field(temperature);

    // Set displacement for temperature data
    disp = 3 * sizeof(int);

    // Set file view and read temperature data
    MPI_File_set_view(fp, 0, MPI_DOUBLE, parallel->filetype, "native",
                      MPI_INFO_NULL);
    MPI_File_read_at_all(fp, disp, temperature->data,
                         1, parallel->restarttype, MPI_STATUS_IGNORE);

    // Close the file
    MPI_File_close(&fp);
}

void copy_field_data(Field *temperature1, Field *temperature2) {
    // Ensure that dimensions match
    assert(temperature1->nx == temperature2->nx);
    assert(temperature1->ny == temperature2->ny);

    // Copy data
    memcpy(temperature2->data, temperature1->data, (temperature1->nx + 2) * (temperature1->ny + 2) * sizeof(double));
}

void swap_field_data(Field *temperature1, Field *temperature2) {
    // Swap data pointers
    double *tmp = temperature1->data;
    temperature1->data = temperature2->data;
    temperature2->data = tmp;
}

void allocate_field(Field *temperature) {
    // Allocate memory for the field, including ghost layers
    temperature->data = allocate_2d_array(temperature->nx + 2, temperature->ny + 2);

    // Initialize the field to zero
    memset(temperature->data, 0.0, (temperature->nx + 2) * (temperature->ny + 2) * sizeof(double));
}

void cleanup_resources(Field *temperature1, Field *temperature2, ParallelData *parallel) {
    // Deallocate the 2D arrays
    delete[] temperature1->data;
    delete[] temperature2->data;

    // Free MPI datatypes
    MPI_Type_free(&parallel->rowtype);
    MPI_Type_free(&parallel->columntype);
    MPI_Type_free(&parallel->subarraytype);
    MPI_Type_free(&parallel->restarttype);
    MPI_Type_free(&parallel->filetype);
}
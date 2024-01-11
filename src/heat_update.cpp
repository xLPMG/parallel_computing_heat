#include <mpi.h>
#include "constants.h"
#include "heat.h"
#include "constants.h"
// Include header files if necessary

void start_halo_exchange(Field *temperature, ParallelData *parallel)
{

    // Width for accessing and navigating through the temperature field
    int width = temperature->ny + 2;

    // (up <-> down)
    // Communication 1: Send data to the upper neighbor and receive from the lower neighbor
    MPI_Isend(&temperature->data[idx(1, 1, width)], 1, parallel->rowtype, parallel->nup, ROW_TAG_UP, parallel->comm, &(parallel->requests[0]));
    MPI_Irecv(&temperature->data[idx(temperature->nx + 1, 1, width)], 1, parallel->rowtype, parallel->ndown, ROW_TAG_UP, parallel->comm, &(parallel->requests[1]));
    
    // (down <-> up)
    // Communication 2: Send data to the lower neighbor and receive from the upper neighbor
    MPI_Isend(&temperature->data[idx(temperature->nx, 1, width)], 1, parallel->rowtype, parallel->ndown, ROW_TAG_DOWN, parallel->comm, &(parallel->requests[2]));
    MPI_Irecv(&temperature->data[idx(0, 1, width)], 1, parallel->rowtype, parallel->nup, ROW_TAG_DOWN, parallel->comm, &(parallel->requests[3]));

    // (left <-> right)
    // Communication 3: Send data to the left neighbor and receive from the right neighbor
    MPI_Isend(&temperature->data[idx(0, 1, width)], 1, parallel->columntype, parallel->nleft, COLUMN_TAG_LEFT, parallel->comm, &(parallel->requests[4]));
    MPI_Irecv(&temperature->data[idx(0, temperature->ny + 1, width)], 1, parallel->columntype, parallel->nright, COLUMN_TAG_LEFT, parallel->comm, &(parallel->requests[5]));

    // (right <-> left)
    // Communication 4: Send data to the right neighbor and receive from the left neighbor
    MPI_Isend(&temperature->data[idx(0, temperature->ny, width)], 1, parallel->columntype, parallel->nright, COLUMN_TAG_RIGHT, parallel->comm, &(parallel->requests[6]));
    MPI_Irecv(&temperature->data[idx(0, 0, width)], 1, parallel->columntype, parallel->nleft, COLUMN_TAG_RIGHT, parallel->comm, &(parallel->requests[7]));
}

void complete_halo_exchange(ParallelData *parallel)
{
    MPI_Status recv_status[8];

    // Wait for the completion of non-blocking communication requests related to halo exchange
    MPI_Waitall(8, parallel->requests, recv_status);
}

void update_interior_temperature(Field *curr, Field *prev, double a, double dt)
{
    //Updated alles außer die äußerste Reihe von Werten
    int i, j;
    int ic, iu, id, il, ir; // Indices for center, up, down, left, right
    int width;
    width = curr->ny + 2;
    double dx2, dy2;
    double x_const, y_const;

    //  Determine the temperature field at the next time step.
    //  As fixed boundary conditions are applied, the outermost grid points are not updated.
    dx2 = prev->dx * prev->dx;
    dy2 = prev->dy * prev->dy;
    x_const = (a * dt / (dx2));
    y_const = (a * dt / (dy2));
    // Loop over the interior grid points for the update
    for (i = 2; i < curr->nx; i++)
    {
        for (j = 2; j < curr->ny; j++)
        {
            ic = idx(i, j, width);
            iu = idx(i, j - 1, width);
            id = idx(i, j + 1, width);
            ir = idx(i + 1, j, width);
            il = idx(i - 1, j, width);

            // Update the temperature using the five-point stencil

            curr->data[ic] = prev->data[ic] 
            + (x_const * (prev->data[ir] - (2 * prev->data[ic]) + prev->data[il])) 
            + (y_const * (prev->data[iu] - (2 * prev->data[ic]) + prev->data[id]));
        }
    }
}

void update_boundary_temperature(Field *curr, Field *prev, double a, double dt)
{
    //Updated nur die äußerste Reihe von Werten
    int i, j;
    int ic, iu, id, il, ir; // Indices for center, up, down, left, right
    int width;
    width = curr->ny + 2;
    double dx2, dy2;
    double x_const, y_const;

    //  Determine the temperature field at the next time step.
    //  As fixed boundary conditions are applied, the outermost grid points are not updated.
    dx2 = prev->dx * prev->dx;
    dy2 = prev->dy * prev->dy;
    x_const = (a * dt / (dx2));
    y_const = (a * dt / (dy2));

    /*Unterschied zu update_interior_temperature ist, dass fuer die verschiedenen Grenzen,
    die Konstantenwerte aus init_heat_field in heat_init.cpp genutzt werden*/

    // Update the left border
    i = 1;
    for (j = 1; j < curr->ny + 1; j++)
    {
        ic = idx(i, j, width);
        iu = idx(i, j - 1, width);
        id = idx(i, j + 1, width);
        ir = idx(i + 1, j, width);
        il = idx(i - 1, j, width);

        // Apply the five-point stencil to update the temperature at the left borders.

        curr->data[ic] = prev->data[ic] 
        + (x_const * (prev->data[ir] - (2 * prev->data[ic]) + prev->data[il])) 
        + (y_const * (prev->data[iu] - (2 * prev->data[ic]) + prev->data[id]));
    }

    // Update the right border
    i = curr->nx;
    for (j = 1; j < curr->ny + 1; j++)
    {
        ic = idx(i, j, width);
        iu = idx(i, j - 1, width);
        id = idx(i, j + 1, width);
        ir = idx(i + 1, j, width);
        il = idx(i - 1, j, width);

        // Apply the five-point stencil to update the temperature at the right borders.

        curr->data[ic] = prev->data[ic] 
        + (x_const * (prev->data[ir] - (2 * prev->data[ic]) + prev->data[il])) 
        + (y_const * (prev->data[iu] - (2 * prev->data[ic]) + prev->data[id]));
    }

    // Update the upper border
    j = 1;
    for (i = 1; i < curr->nx + 1; i++)
    {
        ic = idx(i, j, width);
        iu = idx(i, j - 1, width);
        id = idx(i, j + 1, width);
        ir = idx(i + 1, j, width);
        il = idx(i - 1, j, width);

        // Apply the five-point stencil to update the temperature at the upper boundary.

        curr->data[ic] = prev->data[ic] 
        + (x_const * (prev->data[ir] - (2 * prev->data[ic]) + prev->data[il])) 
        + (y_const * (prev->data[iu] - (2 * prev->data[ic]) + prev->data[id]));
    }

    // Update the lower border
    j = curr->ny;
    for (i = 1; i < curr->nx + 1; i++)
    {
        ic = idx(i, j, width);
        iu = idx(i, j - 1, width);
        id = idx(i, j + 1, width);
        ir = idx(i + 1, j, width);
        il = idx(i - 1, j, width);

        // Apply the five-point stencil to update the temperature at the lower borders.

        curr->data[ic] = prev->data[ic] 
        + (x_const * (prev->data[ir] - (2 * prev->data[ic]) + prev->data[il])) 
        + (y_const * (prev->data[iu] - (2 * prev->data[ic]) + prev->data[id]));
    }
}
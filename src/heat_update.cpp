#include <mpi.h>
#include "constants.h"
#include "heat.h"
#include "constants.h"
// Include header files if necessary

/**
 * @brief Starts the halo exchange.
 *
 * This function should initiate the halo exchange to communicate boundary data between neighboring processes.
 *
 * @param temperature Pointer to the temperature field structure.
 * @param parallel Pointer to the parallel data.
 */
void start_halo_exchange(Field *temperature, ParallelData *parallel)
{
    // Buffer Arrays
    double send_buffer_up[temperature->nx], send_buffer_down[temperature->nx], send_buffer_left[temperature->ny], send_buffer_right[temperature->nx];
    double recv_buffer_up[temperature->nx], recv_buffer_down[temperature->nx], recv_buffer_left[temperature->ny], recv_buffer_right[temperature->nx];

    // Zaehlvariablen, um ueber die Daten zu laufen
    int i, j;

    // Width for accessing and navigating through the temperature field
    int width = temperature->ny + 2;

    // (up <-> down)
    j = 1;
    for (i = 1; i <= temperature->nx; i++)
    {
        send_buffer_up[j - 1] = temperature->data[idx(i, j, width)];
    }
    // Communication 1: Send data to the upper neighbor and receive from the lower neighbor
    MPI_Isend(send_buffer_up, temperature->nx, MPI_DOUBLE, parallel->nup, ROW_TAG_UP, parallel->comm, &(parallel->requests[0]));
    MPI_Irecv(recv_buffer_down, temperature->nx, MPI_DOUBLE, parallel->ndown, ROW_TAG_DOWN, parallel->comm, &(parallel->requests[1]));

    // (down <-> up)
    j = temperature->ny;
    for (i = 1; i <= temperature->nx; i++)
    {
        send_buffer_down[j - 1] = temperature->data[idx(i, j, width)];
    }
    // Communication 2: Send data to the lower neighbor and receive from the upper neighbor
    MPI_Isend(send_buffer_down, temperature->nx, MPI_DOUBLE, parallel->ndown, ROW_TAG_DOWN, parallel->comm, &(parallel->requests[2]));
    MPI_Irecv(recv_buffer_up, temperature->nx, MPI_DOUBLE, parallel->nup, ROW_TAG_UP, parallel->comm, &(parallel->requests[3]));

    // (left <-> right)
    i = 1;
    for (j = 1; j <= temperature->nx; j++)
    {
        send_buffer_left[j - 1] = temperature->data[idx(i, j, width)];
    }
    // Communication 3: Send data to the left neighbor and receive from the right neighbor
    MPI_Isend(send_buffer_left, temperature->ny, MPI_DOUBLE, parallel->nleft, COLUMN_TAG_LEFT, parallel->comm, &(parallel->requests[4]));
    MPI_Irecv(recv_buffer_right, temperature->ny, MPI_DOUBLE, parallel->nright, COLUMN_TAG_RIGHT, parallel->comm, &(parallel->requests[5]));

    // (right <-> left)
    i = temperature->ny;
    for (j = 1; j <= temperature->nx; j++)
    {
        send_buffer_right[j - 1] = temperature->data[idx(i, j, width)];
    }
    // Communication 4: Send data to the right neighbor and receive from the left neighbor
    MPI_Isend(send_buffer_right, temperature->ny, MPI_DOUBLE, parallel->nright, COLUMN_TAG_RIGHT, parallel->comm, &(parallel->requests[6]));
    MPI_Irecv(recv_buffer_left, temperature->ny, MPI_DOUBLE, parallel->nleft, COLUMN_TAG_LEFT, parallel->comm, &(parallel->requests[7]));
}

/**
 * @brief Waits for the completion of the halo exchange.
 *
 * This function is used to wait for the completion of all requests of the halo exchange.
 *
 * @param parallel Pointer to the parallel data.
 */
void complete_halo_exchange(ParallelData *parallel)
{
    MPI_Status recv_status[parallel->size - 1];

    // Wait for the completion of non-blocking communication requests related to halo exchange
    MPI_Waitall(8, parallel->requests, recv_status);
}

/**
 * @brief Updates the interior temperature field.
 *
 * This function should update the interior temperature field based on the five-point stencil.
 *
 * @param curr Pointer to the current field structure.
 * @param prev Pointer to the previous field structure.
 * @param a Thermal diffusivity.
 * @param dt Time step size.
 */
void update_interior_temperature(Field *curr, Field *prev, double a, double dt)
{
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
    for (i = 1; i < curr->nx + 1; i++)
    {
        for (j = 1; j < curr->ny + 1; j++)
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

/**
 * @brief Update the borders of the temperature field.
 *
 * This function should update the border temperature field based on the five-point stencil.
 *
 * @param curr Pointer to the current field structure.
 * @param prev Pointer to the previous field structure.
 * @param a Thermal diffusivity.
 * @param dt Time step size.
 */
void update_boundary_temperature(Field *curr, Field *prev, double a, double dt)
{
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
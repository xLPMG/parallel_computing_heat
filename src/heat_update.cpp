#include <mpi.h>
#include "heat.h"

// Include header files if necessary

void start_halo_exchange(Field *temperature, ParallelData *parallel)
{

    // This function should initiate the halo exchange to communicate boundary data between neighboring processes.

    // Width for accessing and navigating through the temperature field
    int width = temperature->ny + 2;

    // (up <-> down)
    // Communication 1: Send data to the upper neighbor and receive from the lower neighbor
    // This exchanges the ghost cells in the top row of the local temperature field

    // (down <-> up)
    // Communication 2: Send data to the lower neighbor and receive from the upper neighbor
    // This exchanges the ghost cells in the bottom row of the local temperature field

    // (left <-> right)
    // Communication 3: Send data to the left neighbor and receive from the right neighbor
    // This exchanges the ghost cells in the leftmost column of the local temperature field

    // (right <-> left)
    // Communication 4: Send data to the right neighbor and receive from the left neighbor
    // This exchanges the ghost cells in the rightmost column of the local temperature field
}

void complete_halo_exchange(ParallelData *parallel)
{
    // Wait for the completion of non-blocking communication requests related to halo exchange
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

    //  Determine the temperature field at the next time step. 
    //  As fixed boundary conditions are applied, the outermost grid points are not updated.
    dx2 = prev->dx * prev->dx; 
    dy2 = prev->dy * prev->dy;

    // Loop over the interior grid points for the update
    for (i = 1; i < curr->nx + 1; i++)
    {
        for (j = 1; j < curr->ny + 1; j++)
        {
            ic = idx(i, j, width);
            iu = idx(i + 1, j, width);
            id = idx(i - 1, j, width);
            ir = idx(i, j + 1, width);
            il = idx(i, j - 1, width);

            // Update the temperature using the five-point stencil

            // TODO
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

    //  Determine the temperature field at the next time step. 
    //  As fixed boundary conditions are applied, the outermost grid points are not updated.
    dx2 = prev->dx * prev->dx;
    dy2 = prev->dy * prev->dy;

    // Update the left border
    i = 1;
    for (j = 1; j < curr->ny + 1; j++)
    {
        ic = idx(i, j, width);
        iu = idx(i + 1, j, width);
        id = idx(i - 1, j, width);
        ir = idx(i, j + 1, width);
        il = idx(i, j - 1, width);

        // Apply the five-point stencil to update the temperature at the left and right borders.

        // TODO
    }

    // Update the right border
    i = curr->nx;
    for (j = 1; j < curr->ny + 1; j++)
    {
        ic = idx(i, j, width);
        iu = idx(i + 1, j, width);
        id = idx(i - 1, j, width);
        ir = idx(i, j + 1, width);
        il = idx(i, j - 1, width);

        // Apply the five-point stencil to update the temperature at the left and right borders.

        // TODO
    }

    // Update the lower border
    j = 1;
    for (i = 1; i < curr->nx + 1; i++)
    {
        ic = idx(i, j, width);
        iu = idx(i + 1, j, width);
        id = idx(i - 1, j, width);
        ir = idx(i, j + 1, width);
        il = idx(i, j - 1, width);

        // Apply the five-point stencil to update the temperature at the upper and lower borders.

        // TODO
    }

    // Update the upper border
    j = curr->ny;
    for (i = 1; i < curr->nx + 1; i++)
    {
        ic = idx(i, j, width);
        iu = idx(i + 1, j, width);
        id = idx(i - 1, j, width);
        ir = idx(i, j + 1, width);
        il = idx(i, j - 1, width);

        // Apply the five-point stencil to update the temperature at the upper and lower borders.

        // TODO
    }
}
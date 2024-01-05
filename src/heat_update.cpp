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
    double x_const, y_const;
            
    //  Determine the temperature field at the next time step. 
    //  As fixed boundary conditions are applied, the outermost grid points are not updated.
    dx2 = prev->dx * prev->dx; 
    dy2 = prev->dy * prev->dy;
    x_const = (a*dt/ (dx2));
    y_const = (a*dt/ (dx2));
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

            curr->data[ic] = prev->data[ic] 
            + (x_const * (prev->data[iu] - (2 * prev->data[ic]) + prev->data[id])) 
            + (y_const * (prev->data[ir] - (2 * prev->data[ic]) + prev->data[il]));
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
void update_boundary_temperature(Field *curr, Field *prev, double a, double dt, int *coords)
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

    /*Unterschied zu update_interior_temperature ist, dass fuer die verschiedenen Grenzen, 
    die Konstantenwerte aus init_heat_field in heat_init.cpp genutzt werden*/
    
    // Update the left border
    if(coords[1] == 0){
        i = 1;
        for (j = 1; j < curr->ny + 1; j++)
        {
            ic = idx(i, j, width);
            iu = idx(i + 1, j, width);
            id = idx(i - 1, j, width);
            ir = idx(i, j + 1, width);
            il = idx(i, j - 1, width);

            // Apply the five-point stencil to update the temperature at the left borders.

            curr->data[ic] = prev->data[ic] 
                + (x_const * (prev->data[iu] - (2 * prev->data[ic]) + prev->data[id])) 
                + (y_const * (prev->data[ir] - (2 * prev->data[ic]) + 20.0));
        }
    }

    // Update the right border
    if (coords[1] == dims[1] - 1) {
        i = curr->nx;
        for (j = 1; j < curr->ny + 1; j++)
        {
            ic = idx(i, j, width);
            iu = idx(i + 1, j, width);
            id = idx(i - 1, j, width);
            ir = idx(i, j + 1, width);
            il = idx(i, j - 1, width);

            // Apply the five-point stencil to update the temperature at the right borders.

            curr->data[ic] = prev->data[ic] 
                + (x_const * (prev->data[iu] - (2 * prev->data[ic]) + prev->data[id])) 
                + (y_const * (70.0 - (2 * prev->data[ic]) + prev->data[il]));
        }
    }

    // Update the lower border
    if (coords[0] == dims[0] - 1) {
        j = curr->ny;
        for (i = 1; i < curr->nx + 1; i++)
        {
            ic = idx(i, j, width);
            iu = idx(i + 1, j, width);
            id = idx(i - 1, j, width);
            ir = idx(i, j + 1, width);
            il = idx(i, j - 1, width);

            // Apply the five-point stencil to update the temperature at the lower borders.

            curr->data[ic] = prev->data[ic] 
                + (x_const * (prev->data[iu] - (2 * prev->data[ic]) + 5.0)) 
                + (y_const * (prev->data[ir] - (2 * prev->data[ic]) + prev->data[il]));
        }
    }

    // Update the upper border
    if (coords[0] == 0) {
        j = 1;
        for (i = 1; i < curr->nx + 1; i++)
        {
            ic = idx(i, j, width);
            iu = idx(i + 1, j, width);
            id = idx(i - 1, j, width);
            ir = idx(i, j + 1, width);
            il = idx(i, j - 1, width);

            // Apply the five-point stencil to update the temperature at the upper boundary.

            curr->data[ic] = prev->data[ic] 
                + (x_const * (85.0 - (2 * prev->data[ic]) + prev->data[id])) 
                + (y_const * (prev->data[ir] - (2 * prev->data[ic]) + prev->data[il]));
        }
    }
}
#include <mpi.h>
#include "heat.h"

// Include header files if necessary

void start_halo_exchange(Field *temperature, ParallelData *parallel) {
    
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

void complete_halo_exchange(ParallelData *parallel) {
    // Wait for the completion of non-blocking communication requests related to halo exchange
}

void update_interior_temperature(Field *curr, Field *prev, double a, double dt) {

    // This function should update the interior temperature field based on the five-point stencil.

    // Indices for center, up, down, left, right
    // These indices are used for accessing neighboring grid points during the update.

    // Width of the grid (number of columns)
    // The width is used to calculate the indices and navigate through the temperature field.

    // Parameters for the update
    // a: Thermal diffusivity
    // dt: Time step size

    // Determine the temperature field at the next time step
    // As fixed boundary conditions are applied, the outermost grid points are not updated.

    // Loop over the interior grid points for the update
    // Update the temperature using the five-point stencil
}

void update_boundary_temperature(Field *curr, Field *prev, double a, double dt) {
  int i, j;
  int ic, iu, id, il, ir; // Indices for center, up, down, left, right
  int width;
  width = curr->ny + 2;
  double dx2, dy2;

  dx2 = prev->dx * prev->dx; // Determine the temperature field at the next time step. As fixed boundary conditions are applied, the outermost grid points are not updated.
  dy2 = prev->dy * prev->dy;

  // Update the left and right borders
  i = 1;
  for (j = 1; j < curr->ny + 1; j++) {
      ic = idx(i, j, width);
      iu = idx(i + 1, j, width);
      id = idx(i - 1, j, width);
      ir = idx(i, j + 1, width);
      il = idx(i, j - 1, width);

      // Apply the five-point stencil to update the temperature at the left and right borders.
  }

  i = curr->nx;
  for (j = 1; j < curr->ny + 1; j++) {
      // Update indicies using idx
      
      // Apply the five-point stencil to update the temperature at the left and right borders.
  }

  // Update the upper and lower borders
  j = 1;
  for (i = 1; i < curr->nx + 1; i++) {
      // Update indicies using idx

      // Apply the five-point stencil to update the temperature at the upper and lower borders.
  }

  // Update the lower and upper borders
  j = curr->ny;
  for (i = 1; i < curr->nx + 1; i++) {
      // Update indicies using idx

      // Apply the five-point stencil to update the temperature at the upper and lower borders.
  }
}
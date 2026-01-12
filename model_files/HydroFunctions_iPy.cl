//
//  HydroFunctions.cl
//  clMudflat2D
//
//  Created by Johan Van de Koppel on 31/08/2017.
//  Copyright 2017 Johan Van de Koppel. All rights reserved.
//

#ifndef HYDROFUNCTIONS_CL
#define HYDROFUNCTIONS_CL

////////////////////////////////////////////////////////////////////////////////
// Apriori prototyping
////////////////////////////////////////////////////////////////////////////////

float d2_dxy2( __global float* ); // A prototype definition for d2_dxy2
float d_dx(__global float* );
float d_dy( __global float* );
float dO_dx(__global float*, __global float*, __global float* );
float dO_dy( __global float*, __global float*, __global float* );
float d_uh_dx(__global float*, __global float* );
float d_vh_dy( __global float*, __global float* );
void PeriodicBoundaries( __global float* );
void NeumanBoundaries( __global float* );
void DirichletBoundaries( __global float*, float );
void ReflectingBoundaries(__global float*,__global float* );
void PersistentFluxBoundaries(__global float* );

typedef unsigned int Thread_Id;
#define ON              1
#define OFF             0

////////////////////////////////////////////////////////////////////////////////
// Laplacation operator definition, to calculate diffusive terms
////////////////////////////////////////////////////////////////////////////////

float d2_dxy2(__global float* z)
{
    const float dx = dX;  // Forcing dX to become a float
    const float dy = dY;  // Forcing dY to become a float
    
    // Getting thread coordinates on the grid
    const size_t current = get_global_id(0);
    const size_t row     = (size_t)current/Grid_Width;
    const size_t column  = current%Grid_Width;
    
    const size_t left = row * Grid_Width + column-1;
    const size_t right = row * Grid_Width + column+1;
    const size_t top = (row-1) * Grid_Width + column;
    const size_t bottom = (row+1) * Grid_Width + column;
    
    return  (z[left] + z[right ] - 2.0*z[current])/dx/dx +
            (z[top ] + z[bottom] - 2.0*z[current])/dy/dy ;
}

////////////////////////////////////////////////////////////////////////////////
// Gradient operator definitions, to calculate advective fluxes
////////////////////////////////////////////////////////////////////////////////

float d_dx(__global float* z)
{
    const float dx = dX;
    
    const size_t current = get_global_id(0);
    const size_t row     = (size_t)current/(size_t)Grid_Width;
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    const size_t left    = row * Grid_Width + column-1;
    const size_t right   = row * Grid_Width + column+1;
    
    return ( ( z[right] - z[left] )/2.0/dx );
    
}

float d_dy(__global float* pop)
{
    const float dy = dY;
    
    const size_t current = get_global_id(0);
    const size_t row     = (size_t)current/(size_t)Grid_Width;
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    const size_t top     = (row-1) * Grid_Width + column;
    const size_t bottom  = (row+1) * Grid_Width + column;
    
    return ( ( pop[bottom] - pop[top] )/2.0/dy );
    
}

////////////////////////////////////////////////////////////////////////////////
// Gradient operator definitions, to calculate the pressure gradient
////////////////////////////////////////////////////////////////////////////////

float dO_dx(__global float* h, __global float* s,__global float* b)
{
    const float dx = dX;
    
    const size_t current = get_global_id(0);
    const size_t row     = (size_t)current/(size_t)Grid_Width;
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    const size_t left    = row * Grid_Width + column-1;
    const size_t right   = row * Grid_Width + column+1;
    
    return  ( (h[right]+s[right]+b[right])
             - (h[left ]+s[left ]+b[left ]) ) /2.0/dx;
}

float dO_dy(__global float* h, __global float* s,__global float* b)
{
    const float dy = dY;
    
    const size_t current = get_global_id(0);
    const size_t row     = (size_t)current/(size_t)Grid_Width;
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    const size_t top     = (row-1) * Grid_Width + column;
    const size_t bottom  = (row+1) * Grid_Width + column;
    
    return ( (h[bottom]+s[bottom]+b[bottom])
            - (h[top   ]+s[top   ]+b[top] ) ) /2.0/dy;
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
float dO_dx2(__global float* h, __global float* s)
{
    const float dx = dX;
    
    const size_t current = get_global_id(0);
    const size_t row     = (size_t)current/(size_t)Grid_Width;
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    const size_t left    = row * Grid_Width + column-1;
    const size_t right   = row * Grid_Width + column+1;
    
    return  ( (h[right]+s[right])
             - (h[left ]+s[left ]) ) /2.0/dx;
}

float dO_dy2(__global float* h, __global float* s)
{
    const float dy = dY;
    
    const size_t current = get_global_id(0);
    const size_t row     = (size_t)current/(size_t)Grid_Width;
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    const size_t top     = (row-1) * Grid_Width + column;
    const size_t bottom  = (row+1) * Grid_Width + column;
    
    return ( (h[bottom]+s[bottom])
            - (h[top   ]+s[top   ]) ) /2.0/dy;
}
////////////////////////////////////////////////////////////////////////////////
// Definition of the functions that compute the derivatives of uh and vh
////////////////////////////////////////////////////////////////////////////////

float d_uh_dx(__global float* u, __global float* h)
{
    const float dx = dX;
    
    const size_t current = get_global_id(0);
    const size_t row     = (size_t)current/(size_t)Grid_Width;
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    const size_t left    = row * Grid_Width + column-1;
    const size_t right   = row * Grid_Width + column+1;
    
    return ( ( u[right]*h[right] - u[left]*h[left] )/2.0/dx );
}

float d_vh_dy(__global float* v, __global float* h)
{
    const float dy = dY;
    
    const size_t current = get_global_id(0);
    const size_t row     = (size_t)current/(size_t)Grid_Width;
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    const size_t top     = (row-1) * Grid_Width + column;
    const size_t bottom  = (row+1) * Grid_Width + column;
    
    return (( v[bottom]*h[bottom] - v[top]*h[top] )/2.0/dy );
}

////////////////////////////////////////////////////////////////////////////////
// Periodic Boundary conditions function
////////////////////////////////////////////////////////////////////////////////

void PeriodicBoundaries(__global float* z)
{
    const size_t current = get_global_id(0);
    const size_t row     = (size_t)current/(size_t)Grid_Width;
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    if(row==0) // Lower boundary
    {
        z[current] = z[(Grid_Height-2) * Grid_Width+column];
    }
    else if(row==Grid_Height-1) // Upper boundary
    {
        z[current] = z[1*Grid_Width+column];
    }
    else if(column==0) // Left boundary
    {
        z[current] = z[row * Grid_Width + Grid_Width-2];
    }
    else if(column==Grid_Width-1) // Right boundary
    {
        z[current] = z[row * Grid_Width + 1];
    }
}

////////////////////////////////////////////////////////////////////////////////
// Neumann Boundary conditions function, having zero flux on the edge
////////////////////////////////////////////////////////////////////////////////

void NeumanBoundaries(__global float* z)
{
    const size_t current = get_global_id(0);
    const size_t row     = (size_t)current/(size_t)Grid_Width;
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    if(row==0) // Lower boundary
    {
        z[row * Grid_Width + column]=z[1*Grid_Width + column];
    }
    else if(row==Grid_Height-1) // Upper boundary
    {
        z[row * Grid_Width + column]=z[(Grid_Height-2) * Grid_Width + column];
    }
    else if(column==0) // Left boundary
    {
        z[row * Grid_Width + column]=z[row * Grid_Width + 1];
    }
    else if(column==Grid_Width-1) // Right boundary
    {
        z[row * Grid_Width + column]=z[row * Grid_Width + Grid_Width-2];
    }
}

////////////////////////////////////////////////////////////////////////////////
// Reflecting Boundary conditions function, for shallow water equations
////////////////////////////////////////////////////////////////////////////////

void ReflectingBoundaries(__global float* u,__global float* v)
{
    const size_t current = get_global_id(0);
    const size_t row     = (size_t)current/(size_t)Grid_Width;
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    if(row==0) // Lower boundary
    {
        u[row * Grid_Width + column] = u[1*Grid_Width+column];
        v[row * Grid_Width + column] =-v[1*Grid_Width+column]; // [-]
    }
    else if(row==Grid_Height-1) // Upper boundary
    {
        u[row * Grid_Width + column] = u[(Grid_Height-2) * Grid_Width+column];
        v[row * Grid_Width + column] =-v[(Grid_Height-2) * Grid_Width+column]; // [-]
    }
    else if(column==0) // Left boundary
    {
        u[row * Grid_Width + column] =-u[row * Grid_Width + 1];
        v[row * Grid_Width + column] = v[row * Grid_Width + 1];
    }
    else if(column==Grid_Width-1) // Right boundaries
    {
        u[row * Grid_Width + column] =-u[row * Grid_Width + Grid_Width-2];
        v[row * Grid_Width + column] = v[row * Grid_Width + Grid_Width-2];
    }
}

////////////////////////////////////////////////////////////////////////////////
// Persistent Flux Boundary condition function, extrapolating over de boundaries
////////////////////////////////////////////////////////////////////////////////

void PersistentFluxBoundaries(__global float* z)
{
    // Getting thread coordinates on the grid
    const Thread_Id current = get_global_id(0);
    const Thread_Id row     = (Thread_Id)current/(Thread_Id)Grid_Width;
    const Thread_Id column  = (Thread_Id)current%(Thread_Id)Grid_Width;
    
    // current = row * Grid_Width + column;
    
    if(row==0) // Lower boundary
    {
        z[current] = 2*z[(row+1) * Grid_Width + column] - z[ (row+2) * Grid_Width + column];
    }
    else if(row==Grid_Height-1) // Upper boundary
    {
        z[current] = 2*z[(row-1) * Grid_Width + column] - z[(row-2) * Grid_Width + column];
    }
    else if(column==0) // Left boundary
    {
        z[current] = 2*z[row * Grid_Width + column+1] - z[row * Grid_Width + column+2];
    }
    else if(column==Grid_Width-1) // Right boundary
    {
        z[current] = 2*z[row * Grid_Width + column-1] - z[row * Grid_Width + column-2];
    }
}

////////////////////////////////////////////////////////////////////////////////
// Dirichlet Boundary condition function, having fixed values on the edge
////////////////////////////////////////////////////////////////////////////////

void DirichletBoundaries(__global float* z, float Value)
{
    // Getting thread coordinates on the grid
    const size_t current = get_global_id(0);
    const size_t row     = (size_t)current/(size_t)Grid_Width;
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    // current = row * Grid_Width + column;
    
    if(row==0) // Lower boundary
    {
        z[current]=Value;
    }
    else if(row==Grid_Height-1) // Upper boundary
    {
        z[current]=Value;
    }
    else if(column==0) // Left boundary
    {
        z[current]=Value;
    }
    else if(column==Grid_Width-1) // Right boundary
    {
        z[current]=Value;
    }
}

#endif



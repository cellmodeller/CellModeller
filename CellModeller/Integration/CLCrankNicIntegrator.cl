

void userSignalRates(float gridVolume, float area, float volume, const int cellType, float* rates, __global const float* species, __global const float* signals)
{
    %(sigKernel)s
}
void userSpecRates(float gridVolume, float area, float volume, const int cellType, __global float* rates, __global float* species, __global const float* signals)
{
    %(specKernel)s
}

__kernel void gridCells(const float gridOrigx,
                        const float gridOrigy,
                        const float gridOrigz,
                        const float gridSizex,
                        const float gridSizey,
                        const float gridSizez,
                        const int gridDimx,
                        const int gridDimy,
                        const int gridDimz,
                        __global const float4* pos,
                        __global float* weights,
                        __global int* indices)
{
    int id = get_global_id(0);
    int base = id*8;
    int4 idx = 0;
    float4 p = pos[id];
    float ix = floor((p.x - gridOrigx) / gridSizex);
    float iy = floor((p.y - gridOrigy) / gridSizey);
    float iz = floor((p.z - gridOrigz) / gridSizez);
    idx.x = (int)(ix);
    idx.y = (int)(iy);
    idx.z = (int)(iz);



    // Indices of 8 nearest grid nodes
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++) {
            for (int k=0; k<2; k++) {
                int4 offset = {i,j,k,0};
                int4 ind = idx + offset;
                // Is the neighbourhood node in the grid?
                int inGrid = (int)(ind.x>=0 && ind.x<gridDimx && ind.y>=0 && ind.y<gridDimy && ind.z>=0 && ind.z<gridDimz);
                int flatidx = ind.z + ind.y*gridDimz + ind.x*gridDimz*gridDimy;
                // if in grid put index, else zero
                indices[base + k + j*2 + i*4] = inGrid * flatidx;
            }
        }
    }

    // Trilinear weights on 8 corner nodes
    float dx = (p.x - gridOrigx) / gridSizex;
    dx -= floor(dx);
    float dy = (p.y - gridOrigy) / gridSizey;
    dy -= floor(dy);
    float dz = (p.z - gridOrigz) / gridSizez;
    dz -= floor(dz);

    //  Tedious bounds checking on neighbourhood to avoid wrapping and overflow...
    //    Set weight to zero if off grid (index will also be zero)
    
   weights[base+0] /*[0,0,0]*/ = (float)(ix>=0 && ix<gridDimx && iy>=0 && iy<gridDimy && iz>=0 && iz<gridDimz) *  (1.f-dx)*(1.f-dy)*(1.f-dz);
   weights[base+1] /*[0,0,1]*/ = (float)(ix>=0 && ix<gridDimx && iy>=0 && iy<gridDimy && iz>=-1 && iz<gridDimz-1) * (1.f-dx)*(1.f-dy)*dz;
   weights[base+2] /*[0,1,0]*/ = (float)(ix>=0 && ix<gridDimx && iy>=-1 && iy<gridDimy-1 && iz>=0 && iz<gridDimz) * (1.f-dx)*dy*(1.f-dz);
   weights[base+3] /*[0,1,1]*/ = (float)(ix>=0 && ix<gridDimx && iy>=-1 && iy<gridDimy-1 && iz>=-1 && iz<gridDimz-1) * (1.f-dx)*dy*dz;
   weights[base+4] /*[1,0,0]*/ = (float)(ix>=-1 && ix<gridDimx-1 && iy>=0 && iy<gridDimy && iz>=0 && iz<gridDimz) * dx*(1.f-dy)*(1.f-dz);
   weights[base+5] /*[1,0,1]*/ = (float)(ix>=-1 && ix<gridDimx-1 && iy>=0 && iy<gridDimy && iz>=-1 && iz<gridDimz-1) * dx*(1.f-dy)*dz;
   weights[base+6] /*[1,1,0]*/ = (float)(ix>=-1 && ix<gridDimx-1 && iy>=-1 && iy<gridDimy-1 && iz>=0 && iz<gridDimz) * dx*dy*(1.f-dz);
   weights[base+7] /*[1,1,1]*/ = (float)(ix>=-1 && ix<gridDimx-1 && iy>=-1 && iy<gridDimy-1 && iz>=-1 && iz<gridDimz-1) * dx*dy*dz;
    
}

__kernel void setCellSignals(const int numSignals,
                                    const int gridTotalSize,
                                    const int gridDimx,
                                    const int gridDimy,
                                    const int gridDimz,
                                   __global const int* indices,
                                   __global const float* weights,
                                   __global const float* grid,
                                   __global float* levels)
{
    int id = get_global_id(0);
    int idxbase = id*8;
    int lvlbase = id*numSignals;

    for (int s=0; s<numSignals; s++) {
        levels[lvlbase+s] = 0.0;
    }

    // Iterate over 8 nearest grid nodes
    for (int i=0; i<8; i++) {
        int ind = idxbase + i;
        int grid_idx = indices[ind];
        int ix = grid_idx/(gridDimz*gridDimy);
        int iy = (grid_idx - ix*gridDimz*gridDimy)/gridDimz;
        int iz = grid_idx - ix*gridDimz*gridDimy - iy*gridDimz;
        float wt = weights[ind];
        // sum weighted grid values
        for (int s=0; s<numSignals; s++) {
            size_t gsidx = grid_idx + s*gridTotalSize;
            levels[lvlbase+s] += grid[gsidx]*wt;
        }
    }
}


__kernel void speciesRates(const int numSignals,
                           const int numSpecies,
                           const float gridVolume,
                           __global const float* areas,
                           __global const float* volumes,
                           __global const int* celltype,
                           __global const float* cellSpecLevels,
                           __global const float* cellSignalLevels,
                           __global float* specRate)
{
    int id = get_global_id(0);
    int sigbase = id*numSignals;
    int specbase = id*numSpecies;
    int cellType = celltype[id];
    __global const float* species = cellSpecLevels+specbase;
    __global const float* signals = cellSignalLevels+sigbase;
    __global float* rates = specRate+specbase;

    userSpecRates(gridVolume, areas[id], volumes[id], cellType, rates, species, signals);
}



__kernel void signalRates(const int numSignals,
                          const int numSpecies,
                          const float gridVolume,
                          __global const float* areas,
                          __global const float* volumes,
                          __global const int* celltype,
                          __global const float* cellSpecLevels,
                          __global const float* cellSignalLevels,
                          __global const float* weights,
                          __global float* sigRates)
{
    int id = get_global_id(0);
    int base = id*8*numSignals;
    int wbase = id*8;
    int sigbase = id*numSignals;
    int specbase = id*numSpecies;
    int cellType = celltype[id];

    __global const float* species = cellSpecLevels+specbase;
    __global const float* signals = cellSignalLevels+sigbase;

    float cellSigRates[%(nSignals)i];
    userSignalRates(gridVolume, areas[id], volumes[id], cellType, cellSigRates, species, signals);

    // Iterate over 8 nearest grid nodes
    for (int i=0; i<8; i++) {
        float wt = weights[wbase+i];
        int gbase = base + numSignals*i;
        __global float* rates = sigRates+gbase;

        // put cells signal rate
        for (int s=0; s<numSignals; s++) {
            rates[s] = cellSigRates[s]*wt;
        }
    }
}



__kernel void diluteSpecs(const int numSpecies,
                          __global const float* oldVols,
                          __global const float* vols,
                          __global float* specRates)
{
  int id = get_global_id(0);
  float factor = oldVols[id]/vols[id];
  int base = id*numSpecies;
  for (int i = 0; i < numSpecies; i++) {
    specRates[base+i] = specRates[base+i]*factor;
  }
}

__kernel void speciesDT(const int numSpecies,
                        __global float* cellSpecLevels,
                        __global const float* specRate,
                        const float dt)
{
  int id = get_global_id(0); 
  int specbase = id*numSpecies;

  for(int i=0; i < numSpecies; i++)
  { 
    cellSpecLevels[specbase+i] += dt * specRate[specbase+i];
  }
}

__kernel void setCellSigImplicit(__global const int* indices,
                                 __global const float* weights,
                                 __global const float* grid,
                                 __global const float* transport,
                                 __global float* levels,
                                 const int numSignals,
                                 const int gridTotalSize,
                                 const int gridDimx,
                                 const int gridDimy,
                                 const int gridDimz)
{
    int id = get_global_id(0);
    int idxbase = id*8;
    int lvlbase = id*numSignals*2;

    for (int s=0; s<numSignals; s++) {
        levels[lvlbase+s] = 0.0;
        levels[lvlbase+s+numSignals] = 0.0;
    }

    // Iterate over 8 nearest grid nodes
    for (int i=0; i<8; i++) {
        int ind = idxbase + i;
        int grid_idx = indices[ind];
        int ix = grid_idx/(gridDimz*gridDimy);
        int iy = (grid_idx - ix*gridDimz*gridDimy)/gridDimz;
        int iz = grid_idx - ix*gridDimz*gridDimy - iy*gridDimz;
        float wt = weights[ind];
        // sum weighted grid values
        for (int s=0; s<numSignals; s++) {
            size_t gsidx = grid_idx + s*gridTotalSize;
            levels[lvlbase+s] += grid[gsidx]*wt;
            levels[lvlbase+s+numSignals] += transport[gsidx]*wt;
        }
    }
}

__kernel void setCellSignalsImplicit(const int numSignals,
                                     const int gridTotalSize,
                                     const int gridDimx,
                                     const int gridDimy,
                                     const int gridDimz,
                                     __global const int* indices,
                                     __global const float* weights,
                                     __global const float* grid,
                                     __global const float* transport,
                                     __global float* levels)
{
    int id = get_global_id(0);
    int idxbase = id*8;
    int lvlbase = id*numSignals*2;

    for (int s=0; s<numSignals; s++) {
        levels[lvlbase+s] = 0.0;
        levels[lvlbase+s+numSignals] = 0.0;
    }

    // Iterate over 8 nearest grid nodes
    for (int i=0; i<8; i++) {
        int ind = idxbase + i;
        int grid_idx = indices[ind];
        int ix = grid_idx/(gridDimz*gridDimy);
        int iy = (grid_idx - ix*gridDimz*gridDimy)/gridDimz;
        int iz = grid_idx - ix*gridDimz*gridDimy - iy*gridDimz;
        float wt = weights[ind];
        // sum weighted grid values
        for (int s=0; s<numSignals; s++) {
            size_t gsidx = grid_idx + s*gridTotalSize;
            levels[lvlbase+s] += grid[gsidx]*wt;
            levels[lvlbase+s+numSignals] += transport[gsidx]*wt;
        }
    }
}

__kernel void speciesRatesImplicit(const int numSignals,
                                   const int numSpecies,
                                   const float gridVolume,
                                   __global const float* areas,
                                   __global const float* volumes,
                                   __global const int* celltype,
                                   __global float* cellSpecLevels,
                                   __global const float* cellSignalLevels,
                                   __global float* specRate)
{
    int id = get_global_id(0);
    int sigbase = id*numSignals*2;
    int specbase = id*numSpecies;
    int cellType = celltype[id];
    __global float* species = cellSpecLevels+specbase;
    __global const float* signals = cellSignalLevels+sigbase;
    __global float* rates = specRate+specbase;

    userSpecRates(gridVolume, areas[id], volumes[id], cellType, rates, species, signals);
}

__kernel void signalRatesImplicit(const int numSignals,
                          const int numSpecies,
                          const float gridVolume,
                          __global const float* areas,
                          __global const float* volumes,
                          __global const int* celltype,
                          __global const float* cellSpecLevels,
                          __global const float* cellSignalLevels,
                          __global const float* weights,
                          __global float* sigRates)
{
    int id = get_global_id(0);
    int base = id*8*numSignals;
    int wbase = id*8;
    int sigbase = id*numSignals*2;
    int specbase = id*numSpecies;
    int cellType = celltype[id];

    __global const float* species = cellSpecLevels+specbase;
    __global const float* signals = cellSignalLevels+sigbase;

    float cellSigRates[%(nSignals)i];
    userSignalRates(gridVolume, areas[id], volumes[id], cellType, cellSigRates, species, signals);

    // Iterate over 8 nearest grid nodes
    for (int i=0; i<8; i++) {
        float wt = weights[wbase+i];
        int gbase = base + numSignals*i;
        __global float* rates = sigRates+gbase;

        // put cells signal rate
        for (int s=0; s<numSignals; s++) {
            rates[s] = cellSigRates[s]*wt;
        }
    }
}

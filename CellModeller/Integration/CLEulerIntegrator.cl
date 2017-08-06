void userSpecRates(float4 pos, float area, float volume, const int cellType, const float effGrowth, __global float* rates, __global const float* species)
{
    %s
}

__kernel void speciesRates(const int numSpecies,
                           __global const float4* pos,
                           __global const float* areas,
                           __global const float* volumes,
                           __global const int* celltype,
                           __global const float* effGrow,
                           __global const float* cellSpecLevels,
                           __global float* specRate)
{
    int id = get_global_id(0);
    int specbase = id*numSpecies;
    int cellType = celltype[id];
    __global const float* species = cellSpecLevels+specbase;
    __global float* rates = specRate+specbase;

    userSpecRates(pos[id], areas[id], volumes[id], cellType, effGrow[id], rates, species);
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

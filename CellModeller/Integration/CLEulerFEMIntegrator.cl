

void userSignalRates(float gridVolume, float area, float volume, const int cellType, __global float* rates, __global const float* species, __global const float* signals)
{
    %(sigKernel)s
}
void userSpecRates(float gridVolume, float area, float volume, const int cellType, __global float* rates, __global float* species, __global const float* signals)
{
    %(specKernel)s
}

__kernel void speciesRates(const int numSignals,
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
    int sigbase = id*numSignals;
    int specbase = id*numSpecies;
    int cellType = celltype[id];
    __global float* species = cellSpecLevels+specbase;
    __global const float* signals = cellSignalLevels+sigbase;
    __global float* rates = specRate+specbase;

    userSpecRates(gridVolume, areas[id], volumes[id], cellType, rates, species, signals);
    for (int i=0; i<numSpecies; i++)
    {
    	species[i] += rates[i];
    }
}



__kernel void signalRates(const int numSignals,
                          const int numSpecies,
                          const float gridVolume,
                          __global const float* areas,
                          __global const float* volumes,
                          __global const int* celltype,
                          __global const float* cellSpecLevels,
                          __global const float* cellSignalLevels,
			  __global float* cellSignalRates)
{
    int id = get_global_id(0);
    int sigbase = id*numSignals;
    int specbase = id*numSpecies;
    int cellType = celltype[id];

    __global const float* species = cellSpecLevels+specbase;
    __global const float* signals = cellSignalLevels+sigbase;
    __global float* rates = cellSignalRates + sigbase;

    userSignalRates(gridVolume, areas[id], volumes[id], cellType, rates, species, signals);
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

__kernel void combine_grad_components(__global const float* grad_x,
					__global const float* grad_y,
					__global const float* grad_z,
					__global float4* grad)
{
	int i = get_global_id(0);
	grad[i].s0 = grad_x[i];
	grad[i].s1 = grad_y[i];
	grad[i].s2 = grad_z[i];
	grad[i].s3 = 0.f;
}


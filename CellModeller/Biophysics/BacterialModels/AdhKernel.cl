float userAdhLogic(const float adh_str1, const float adh_str2)
{
    %s
}

__kernel void build_Tmatrix(const int max_contacts,
                            const float muA,
                            const float gamma,
                            __global const float4* centers,
                            __global const float4* dirs,
                            __global const float* lens,
                            __global const float* rads,
                            __global const float* adh_str,
                            __global const int* n_cts,
                            __global const int* frs,
                            __global const int* tos,
                            __global const float4* pts,
                            __global const float4* norms,
                            __global float8* ct_tangs_fr,
                            __global float8* ct_tangs_to,
                            __global float* overlap,
                            __global float* ct_adh_str)
{
    int id = get_global_id(0);
    int ct = get_global_id(1);
    
    if (ct >= n_cts[id]) return;
    
    int i = id*max_contacts + ct;
    
    int a = frs[i];
    float4 r_a = pts[i]-centers[a];
    
    float4 Ia[4];
    cyl_inv_inertia_tensor(muA, lens[a]+2.f*rads[a], dirs[a], Ia);
    
    int b = tos[i];
    // plane contacts have no to_ent, and have negative indices - and don't do adhesion
    if (b < 0) {
        return;
    }
    
    float4 r_b = pts[i]-centers[b];
    
    float4 Ib[4];
    cyl_inv_inertia_tensor(muA, lens[b]+2.f*rads[b], dirs[b], Ib);
    
    //tangents to contacts
    
    float8 ct_tang_fr = 0.f;
    float8 ct_tang_to = 0.f;
    
    float4 d1 = 0.f;//norms[i] * dot(norms[i], dirs[a]) - dirs[a]; //along the axis of the fr cell
    float4 z_axis = {0.f, 0.f, 1.f, 0.f};
    //we do need a special case here, in case norm and cell are parallel
    //if(length(d1)<MARGIN){
    d1 = cross(norms[i],z_axis); //this is general for the 2d case, but not otherwise.
    //}
    d1 = normalize(d1);
    
    float4 dxr_a = 0.f;
    float4 dxr_b = 0.f;
    
    dxr_a = cross(d1, r_a);
    dxr_b = cross(d1, r_b);
    
    ct_tang_fr.s012 = d1.s012/(muA*(lens[a]+2.f*rads[a]));
    ct_tang_fr.s3 = -dot(dxr_a, Ia[0]);
    ct_tang_fr.s4 = -dot(dxr_a, Ia[1]);
    ct_tang_fr.s5 = -dot(dxr_a, Ia[2]);
    ct_tang_fr.s6 = (1.f/gamma) * dot(dirs[a], r_a) * dot(dirs[a], d1)/(lens[a]+2.f*rads[a]);
    
    ct_tang_to.s012 = d1.s012/(muA*(lens[b]+2.f*rads[b]));
    ct_tang_to.s3 = -dot(dxr_b, Ib[0]);
    ct_tang_to.s4 = -dot(dxr_b, Ib[1]);
    ct_tang_to.s5 = -dot(dxr_b, Ib[2]);
    ct_tang_to.s6 = (1.f/gamma) * dot(dirs[b], r_b) * dot(dirs[b], d1)/(lens[b]+2.f*rads[b]);
    
    ct_tangs_fr[i] = ct_tang_fr;
    ct_tangs_to[i] = ct_tang_to;
    
    //this sets the user defined adhesion strength of the contact - this should be the place where adhesion should be made proportional to overlap length (if cells are parallel)
    
    ct_adh_str[i] = overlap[i]*userAdhLogic(adh_str[a], adh_str[b]);
}

__kernel void calculate_adhE(const int max_contacts,
                             __global const int* n_cts,
                             __global const int* n_cell_tos,
                             __global const int* cell_tos,
                             __global const float8* ct_tangs_fr,
                             __global const float8* ct_tangs_to,
                             __global const float* dist,
                             __global const float* ct_adh_str,
                             __global float8* adhE)
{
    int i = get_global_id(0);
    int base = i*max_contacts;
    float8 res = 0.f;
    for (int k = base; k < base+n_cts[i]; k++) {
        float8 oldres = res;
        res += ct_tangs_fr[k]*dist[k]*ct_adh_str[k];
    }
    for (int k = base; k < base+n_cell_tos[i]; k++) {
        int n = cell_tos[k];
        if (n < 0) continue;
        res -= ct_tangs_to[n]*dist[n]*ct_adh_str[n];
    }
    adhE[i] = res;
}


#define NFIX 0 

//---
//
__kernel void computeCellAxis(__global const int *cell_walls,
                                __global const uint *wall_nodes, 
                                __global const float4 *nodes,
                                __global const float4 *cell_centre,
                                __global float4 *cell_axis,
                                const int max_cell_walls)
{
    /* Compute principal components of cell vertices 
        - eigenvectors of covariance matrix
        - take vec with smallest eig value as axis (ie. minor axis)
    */

    int c = get_global_id(0);
    int cwid = c*(max_cell_walls+1);
    float4 ctr = cell_centre[c];
    float4 Cmat = 0.f;

    // Covariance matrix
    int nNodes=0;
    for (int i=0, w=cell_walls[cwid]; i<max_cell_walls && w>=0; ++i, w=cell_walls[cwid+i], ++nNodes)
    {
      uint wn1 = wall_nodes[w*2];
      uint wn2 = wall_nodes[w*2+1];
      float4 pos1 = nodes[wn1]-ctr;
      float4 pos2 = nodes[wn2]-ctr;
      Cmat.x += pos1.x*pos1.x + pos2.x*pos2.x;
      Cmat.y += pos1.x*pos1.y + pos2.x*pos2.y;
      Cmat.z = Cmat.y;
      Cmat.w += pos1.y*pos1.y + pos2.y*pos2.y;
    }

    // determinant
    float D = Cmat.x*Cmat.w - Cmat.y*Cmat.z;
    // trace
    float T = Cmat.x + Cmat.w;

    // eigenvalues
    float v1 = T*0.5f + sqrt(T*T/4-D);
    float v2 = T*0.5f - sqrt(T*T/4-D);
    
    // eigenvectors
    float4 e11=0.f;
    e11.x = v1 - Cmat.w;
    e11.y = Cmat.z;
    float4 e12=0.f;
    e12.y=1.f;

    float4 e21=0.f;
    e21.x = v2 - Cmat.w;
    e21.y = Cmat.z;
    float4 e22=0.f;
    e22.x=1.f;

    float4 e1;
    bool b = fabs(Cmat.y)>0.f; //(fabs(Cmat.y*Cmat.y/D)>0.f1);
    e1 = e11*(float)b + e12*(float)(!b);
    float4 e2;
    e2 = e21*(float)b + e22*(float)(!b);
    
    // Choose vector with largest eigenvalue
    b = (fabs(v1)<fabs(v2));
    cell_axis[c] = normalize( e1*(float)b + e2*(float)(!b) ); 
}


//---
//
__kernel void computeCellGeometry(__global const int *cell_walls,
                                __global const uint *wall_nodes, 
                                __global const int *wall_cells, 
                                __global const float4 *nodes,
                                __global float *cell_area,
                                __global float4 *cell_centre,
                                __global float *cell_perim,
                                const int max_cell_walls)
{
    int c = get_global_id(0);
    int cwid = c*(max_cell_walls+1);
    float A = 0.f, L = 0.f;
    float4 ctr = 0.f;

    for (int i=0, w=cell_walls[cwid]; i<max_cell_walls && w>=0; ++i, w=cell_walls[cwid+i])
    {
        uint wn1 = wall_nodes[w*2];
        uint wn2 = wall_nodes[w*2+1];
        float4 p1 = nodes[wn1];
        float4 p2 = nodes[wn2];
        float4 wvec = p2-p1;
        wvec.s23 = 0.f;
        float4 wctr =  0.5f*(p1+p2);
        wctr.s23 = 0.f;

        int wc1 = wall_cells[w*2];
        int wc2 = wall_cells[w*2+1];
        bool inner = (wc1==c);//---cell is inner cell
        // nodes are wrong winding if cell is outside, hence subtract
        A += 0.5f*cross(wctr,wvec).z * (1.f*inner - 1.f*(!inner));
        float dL = length(wvec);
        ctr += wctr*dL; // centre as length weighted average of wall centre pos (~CofM)
        L += dL;
    }
    cell_area[c] = A;
    cell_centre[c] = ctr/L;
    cell_perim[c] = L;
}

//---
//
__kernel void energyDeriv(__global const float4 *nodes, 
                __global const int *node_walls, 
                __global const uint *wall_nodes,
                __global const int *wall_cells, 
                __global const float *cell_area,
                __global const float *cell_area0,
                __global const float *cell_perim,
                __global float4 *deriv,
                const float Karea, const float Kperim,
                const uint max_nbs)
{
    int gid = get_global_id(0);
    float4 x = nodes[gid];
    int nwid = gid*(max_nbs+1); // stride is max_nbs+1 for the -1 at end

    // z unit vector for cross product to get normal to wall
    float4 z = 0.f;
    z.z = -1.f;

    float4 dE = 0.f;
    for (int i=0, w=node_walls[nwid]; i<max_nbs && w>=0; ++i, w=node_walls[nwid+i])
    {
        // Get wall node that is not this node
        // avoid divergence by using bool to test which node to use
        // then blend
        int wn1 = wall_nodes[w*2];
        int wn2 = wall_nodes[w*2+1];
        bool p1 = (wn1!=gid);
        float4 p = nodes[wn1]*(float)p1 + nodes[wn2]*(float)(!p1);

        // wall length and normal
        float L = length(x-p);
        float4 walldir = p-x; //only take position part not angle (w)
        walldir.s23 =0.f;
        walldir = normalize(walldir);
        float4 wvec = nodes[wn2] - nodes[wn1];
        wvec.s23 = 0.f;
        float4 wallnorm = normalize(cross(z,wvec));

        // Cell ids on either side of wall
        int wc1 = wall_cells[w*2];
        int wc2 = wall_cells[w*2+1];

        if (wc1>=0)
        {
            // Approximate area term
            dE += Karea * (cell_area[wc1] - cell_area0[wc1]) * L * wallnorm;
            // Perimeter term
            dE += Kperim * (cell_perim[wc1]) * -walldir; // negative because we want direction toward this node
        }
        if (wc2>=0)
        {
            // Approximate area term
            // wc2 is outside cell, hence -ve
            dE -= Karea * (cell_area[wc2] - cell_area0[wc2]) * L * wallnorm;
            // Perimeter term
            dE += Kperim * (cell_perim[wc2]) * -walldir; // negative because we want direction toward this node
        }
    }
    deriv[gid] = dE;
}

//---
//
__kernel void computeWallLength(__global const uint *wall_nodes, 
                                __global const float4 *nodes,
                                __global float *wall_length)
{
    int w = get_global_id(0);
    uint wn1 = wall_nodes[w*2];
    uint wn2 = wall_nodes[w*2+1];
    float4 d = nodes[wn1] - nodes[wn2];
    d.s23 = 0.f;
    //float4 p1 = nodes[wn2];
    wall_length[w] = length(d); //distance(p0,p1); //sqrt(dot(d,d)); 
}

//---
//
__kernel void Ax(__global const float4 *nodes, __global const float4 *disp,
        __global const int *node_walls, __global const uint *wall_nodes,
        __global const float *wall_stiff, __global const float *wall_length, 
        __global float4 *force, 
        const float S, const float I, const float dt,
        const uint max_nbs, const float alpha)
{
    int gid = get_global_id(0);
    force[gid] = 0.f;
    if (gid<NFIX) return;
    float4 u = disp[gid];
    float4 pos = nodes[gid];
    int nwid = gid*(max_nbs+1); // stride is max_nbs+1 for the -1 at end
    
    // z unit vector for cross product to get normal to wall
    float4 z = 0.f;
    z.z = -1.f;

    float4 f = 0.f;
    for (int i=0, w=node_walls[nwid]; i<max_nbs && w>=0; ++i, w=node_walls[nwid+i])
    {
        // Get wall node that is not this node
        // avoid divergence by using bool to test which node to use
        // then blend
        int wn1 = wall_nodes[w*2];
        int wn2 = wall_nodes[w*2+1];
        bool p1 = (wn1!=gid) && (wn1>=NFIX);
        bool p2 = (wn2!=gid) && (wn2>=NFIX);
        bool pp1 = (wn1!=gid);
        float4 nu = disp[wn1]*(float)p1 + disp[wn2]*(float)(p2);
        float4 otherend = nodes[wn1]*(float)pp1 + nodes[wn2]*(float)(!pp1);
        float4 avg = 0.5f*(otherend + pos);

        /* The divergent version looks like this:
           if (wn1!=gid)
           {
           nu = disp[wn1];
           otherend = nodes[wn1];
           }else
           {
           nu = disp[wn2];
           otherend = nodes[wn2];
           }*/
        float L = wall_length[w];
        float4 walldir = otherend-pos; //only take position part not angle (z)
        walldir.s23 =0.f;
        walldir = normalize(walldir);

        float4 wvec = nodes[wall_nodes[w*2+1]] - nodes[wall_nodes[w*2]];
        wvec.s23 = 0.f;
        float4 wallnorm = normalize(cross(z,walldir));

        // sum wall forces
        float stiff = /*((float)(avg.x>0.f)*0.9 + 0.1f) * */ wall_stiff[w]*S/dt;
        float bend = wall_stiff[w]*I/dt;
        float4 posDispDiff = u-nu;
        posDispDiff.s23=0.f;
        float angSum = nu.z+u.z;
        f += walldir*dot(walldir,posDispDiff)*stiff/L ;
        f += wallnorm * ( dot(wallnorm,posDispDiff)*12.f + 6.f*angSum*L ) * bend/(L*L*L);
        f.z += ( dot(wallnorm,posDispDiff)*6.f + 4.f*u.z*L + 2.f*nu.z*L ) * bend/(L*L);
    }
    force[gid] = f + alpha*u; // (A + \alpha I)x, ie. Ax + Tikhonov regularisation
}

__kernel void rhs(__global const float4 *nodes,
        __global const int *node_walls, __global const int *wall_cells,  
        __global const uint *wall_nodes, __global const float *cell_turgor, 
        __global float4 *force, __global const float4 *disp, const float alpha, const uint max_nbs)
{
    int gid = get_global_id(0);
    force[gid] = 0.f;
    if (gid<NFIX) return;
    float4 pos = nodes[gid];

    // z unit vector for cross product to get normal to wall
    float4 z = 0.f;
    z.z = -1.f;

    // Start index into node_walls array
    int nwid = gid*(max_nbs+1); // stride is max_nbs+1 for the -1 at end


    // Loop over walls
    for (int i=0, w=node_walls[nwid]; i<max_nbs && w>=0; ++i, w=node_walls[nwid+i])
    {
        // Compute length*normal for wall
        float4 wvec = nodes[wall_nodes[w*2+1]] - nodes[wall_nodes[w*2]];
        wvec.s23 = 0.f;
        float4 turg = cross(z,wvec);

        // c1 and c2 are cells on either side of wall
        // and their turgor forces will point in opposite directions
        int c1 = wall_cells[w*2];
        int c2 = wall_cells[w*2+1];

        // Can probably avoid these if statements using bool (see Ax() above)
        // Just need to make sure don't cause illegal access from -ve index
        // (Pad array with a dummy entry at [0], and start from [1]?)
        if (c1>=0)
        {
            force[gid] += turg * cell_turgor[c1];
        }
        if (c2>=0)
        {
            force[gid] -= turg * cell_turgor[c2];
        }
    }
    force[gid] += alpha*disp[gid];
}


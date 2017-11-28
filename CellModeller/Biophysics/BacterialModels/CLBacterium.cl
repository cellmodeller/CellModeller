#define EPSILON 0.1f
#define MARGIN 0.01f
#define RHS_FRAC 1.f 
#define ANG_LIMIT ((float)(5.f*3.14159f/180.f))
#define ISQRT2 ((float)(1.f/sqrt(2.f)))
#define POINT 0.01f

// multiply a matrix and a vector
//  m -- 4x4 matrix
//  v -- 4-vector
float4 matmul(float4 m[], float4 v) {
  return m[0]*v.s0 + m[1]*v.s1 + m[2]*v.s2 + m[3]*v.s3;
}


// transpose a matrix
//  m -- matrix to transpose (vector of columns)
//  t -- result
void transpose(float4 m[], float4 t[]) {
  t[0].s0 = m[0].s0;
  t[0].s1 = m[1].s0;
  t[0].s2 = m[2].s0;
  t[0].s3 = m[3].s0;
  t[1].s0 = m[0].s1;
  t[1].s1 = m[1].s1;
  t[1].s2 = m[2].s1;
  t[1].s3 = m[3].s1;
  t[2].s0 = m[0].s2;
  t[2].s1 = m[1].s2;
  t[2].s2 = m[2].s2;
  t[2].s3 = m[3].s2;
  t[3].s0 = m[0].s3;
  t[3].s1 = m[1].s3;
  t[3].s2 = m[2].s3;
  t[3].s3 = m[3].s3;
}


// multiply two 4x4 matrices
//  a -- left matrix
//  b -- right matrix
//  r -- result
// this is the naive way
void matmulmat(float4 a[], float4 b[], float4 r[]) {
  r[0].s0 = a[0].s0 * b[0].s0 + a[1].s0 * b[0].s1 + a[2].s0 * b[0].s2 + a[3].s0 * b[0].s3;
  r[1].s0 = a[0].s0 * b[1].s0 + a[1].s0 * b[1].s1 + a[2].s0 * b[1].s2 + a[3].s0 * b[1].s3;
  r[2].s0 = a[0].s0 * b[2].s0 + a[1].s0 * b[2].s1 + a[2].s0 * b[2].s2 + a[3].s0 * b[2].s3;
  r[3].s0 = a[0].s0 * b[3].s0 + a[1].s0 * b[3].s1 + a[2].s0 * b[3].s2 + a[3].s0 * b[3].s3;
  r[0].s1 = a[0].s1 * b[0].s0 + a[1].s1 * b[0].s1 + a[2].s1 * b[0].s2 + a[3].s1 * b[0].s3;
  r[1].s1 = a[0].s1 * b[1].s0 + a[1].s1 * b[1].s1 + a[2].s1 * b[1].s2 + a[3].s1 * b[1].s3;
  r[2].s1 = a[0].s1 * b[2].s0 + a[1].s1 * b[2].s1 + a[2].s1 * b[2].s2 + a[3].s1 * b[2].s3;
  r[3].s1 = a[0].s1 * b[3].s0 + a[1].s1 * b[3].s1 + a[2].s1 * b[3].s2 + a[3].s1 * b[3].s3;
  r[0].s2 = a[0].s2 * b[0].s0 + a[1].s2 * b[0].s1 + a[2].s2 * b[0].s2 + a[3].s2 * b[0].s3;
  r[1].s2 = a[0].s2 * b[1].s0 + a[1].s2 * b[1].s1 + a[2].s2 * b[1].s2 + a[3].s2 * b[1].s3;
  r[2].s2 = a[0].s2 * b[2].s0 + a[1].s2 * b[2].s1 + a[2].s2 * b[2].s2 + a[3].s2 * b[2].s3;
  r[3].s2 = a[0].s2 * b[3].s0 + a[1].s2 * b[3].s1 + a[2].s2 * b[3].s2 + a[3].s2 * b[3].s3;
  r[0].s3 = a[0].s3 * b[0].s0 + a[1].s3 * b[0].s1 + a[2].s3 * b[0].s2 + a[3].s3 * b[0].s3;
  r[1].s3 = a[0].s3 * b[1].s0 + a[1].s3 * b[1].s1 + a[2].s3 * b[1].s2 + a[3].s3 * b[1].s3;
  r[2].s3 = a[0].s3 * b[2].s0 + a[1].s3 * b[2].s1 + a[2].s3 * b[2].s2 + a[3].s3 * b[2].s3;
  r[3].s3 = a[0].s3 * b[3].s0 + a[1].s3 * b[3].s1 + a[2].s3 * b[3].s2 + a[3].s3 * b[3].s3;
}


// return the inverse of a quaternion
//  q -- a quaternion
float4 quat_inv(float4 q) {
  float l2 = dot(q, q);
  if (l2==0.f)
    return(q);
  float4 inv = {-q.x/l2, -q.y/l2, -q.z/l2, q.w/l2};
  return inv;
}


// multiply two quaternions
//  a -- a quaternion
//  b -- a quaternion
float4 quat_prod(float4 a, float4 b) {
  float4 res;
  res.x =  a.x*b.w + a.y*b.z - a.z*b.y + a.w*b.x;
  res.y = -a.x*b.z + a.y*b.w + a.z*b.x + a.w*b.y;
  res.z =  a.x*b.y - a.y*b.x + a.z*b.w + a.w*b.z;
  res.w = -a.x*b.x - a.y*b.y - a.z*b.z + a.w*b.w;
  return res;
}


// rotate a vector using a quaternion
//  q -- quaternion
//  v -- vector, w ('real') component is ignored
float4 quat_rot(float4 q, float4 v) {
  v.w = 0.f;
  float4 qi = quat_inv(q);
  float4 v_prime = quat_prod(q, quat_prod(v, qi));
  v_prime.w = 0.f;
  return v_prime;
}


// rotate a vector
//  axis -- axis about which to rotate
//  angle -- amount to rotate (radians)
//  v -- the vector to rotate
float4 rot(float4 axis, float angle, float4 v) {
  float4 q = axis*sin(angle/2.f);
  q.w = cos(angle/2.f);
  return quat_rot(q, v);
}


void print_matrix(float4* m) {
  /*
  printf("\
 % 3.3f  % 3.3f  % 3.3f (% 3.3f)\n\
 % 3.3f  % 3.3f  % 3.3f (% 3.3f)\n\
 % 3.3f  % 3.3f  % 3.3f (% 3.3f)\n\
(% 3.3f)(% 3.3f)(% 3.3f)(% 3.3f)\n\n",
         m[0].x, m[1].x, m[2].x, m[3].x,
         m[0].y, m[1].y, m[2].y, m[3].y,
         m[0].z, m[1].z, m[2].z, m[3].z,
         m[0].w, m[1].w, m[2].w, m[3].w);
  */
}

void cyl_inv_inertia_tensor(float muA, float l, float4 axis, float4 res[]) {
  // first find the inv inertia tensor for a capsule along the x axis
  float diag = 12.f/(muA*l*l*l);

  // now find the matrix that transforms from the x-axis to the axis
  float4 x_axis = {1.f, 0.f, 0.f, 0.f};
  float4 y_axis = {0.f, 1.f, 0.f, 0.f};
  float4 z_axis = {0.f, 0.f, 1.f, 0.f};
  float rot_ang = acos(dot(x_axis, axis));

  float4 y_prime = y_axis;
  float4 z_prime = z_axis;

  if (rot_ang > EPSILON) {
    y_prime = rot(z_prime, rot_ang, y_axis);
    z_prime = normalize(cross(x_axis, axis));
  }

  float4 M[4];
  M[0] = axis;
  M[1] = y_prime;
  M[2] = z_prime;
  M[3] = 0.f;

  float4 MD[4];
  MD[0] = 0.f;
  MD[1] = M[1]*diag;
  MD[2] = M[2]*diag;
  MD[3] = 0.f;

  float4 MDT[4];
  transpose(MD, MDT);

  // I in the world is = M(I_D)M^T
  matmulmat(M, MDT, res);
}



// Set the sq of each cell based on its position.
__kernel void bin_cells(const int grid_x_min,
                        const int grid_x_max,
                        const int grid_y_min,
                        const int grid_y_max,
                        const float grid_spacing,
                        __global const float4* centers,
                        __global int* sqs)
{
  int i = get_global_id(0);
  int x = (int)floor(centers[i].x / grid_spacing) - grid_x_min;
  int y = (int)floor(centers[i].y / grid_spacing) - grid_y_min;
  sqs[i] = y*(grid_x_max-grid_x_min) + x;
}


// find the closest points on two line segments
void closest_points_on_segments(const float4 r_a,  // center of first segment
                                const float4 r_b,  // center of second segment
                                const float4 a,    // direction of first segment (unit)
                                const float4 b,    // direction of second segment (unit)
                                const float len_a, // length of first segment
                                const float len_b, // length of second segment
                                float4* p_a,       // (return) point on first segment
                                float4* p_b,       // (return) point on second segment
                                float4* p_a2,      // (return) 2nd point on first segment
                                float4* p_b2,      // (return) 2nd point on second segment
                                int* two_pts)      // (return) were two points picked?
{
  float hlen_a = len_a / 2.f;
  float hlen_b = len_b / 2.f;
  float4 r = r_b - r_a;
  float a_dot_r = dot(a, r);
  float b_dot_r = dot(b, r);
  float a_dot_b = dot(a, b);
  float denom = 1.f - a_dot_b * a_dot_b;

  float t_a = 0.f;
  float t_b = 0.f;

  *two_pts = 0;
  float t_a2 = 0.f;
  float t_b2 = 0.f;

  if (sqrt(denom) > EPSILON) {
    // non-parallel lines

    // closest points on the same lines if they were infinitely long
    float t_a0 = (a_dot_r - b_dot_r * a_dot_b) / denom;
    float t_b0 = (a_dot_r * a_dot_b - b_dot_r) / denom;

    // there are a few different cases we have to handle...
    bool on_a = fabs(t_a0) < hlen_a;
    bool on_b = fabs(t_b0) < hlen_b;
    if (!on_a && !on_b) {
      // the corner
      float c_a = copysign(hlen_a, t_a0);
      float c_b = copysign(hlen_b, t_b0);

      // signs of partials at the corner
      float dd_dt_a = 2.f*(c_a - a_dot_b*c_b - a_dot_r);
      float dd_dt_b = 2.f*(c_b - a_dot_b*c_a + b_dot_r);

      if (sign(dd_dt_a) == sign(c_a)) {
        // on the other edge
        t_b = c_b;
        t_a = clamp(t_b*a_dot_b + a_dot_r, -hlen_a, hlen_a);
      } else {
        t_a = c_a;
        t_b = clamp(t_a*a_dot_b - b_dot_r, -hlen_b, hlen_b);
      }
    } else if (on_a && !on_b) {
      t_b = copysign(hlen_b, t_b0);  // +/- end of b?
      t_a = clamp(t_b*a_dot_b + a_dot_r, -hlen_a, hlen_a);
    } else if (!on_a && on_b) {
      t_a = copysign(hlen_a, t_a0);  // +/- end of a?
      t_b = clamp(t_a*a_dot_b - b_dot_r, -hlen_b, hlen_b);
    } else {
      t_a = t_a0;
      t_b = t_b0;
    }
  } else {
    // lines are roughly parallel, this case is degenerate
    // start off assuming the lines are in the same direction
    // project a onto b

    // use the same _dot_r for each for consistency
    float x_dot_r = copysign(min(a_dot_r, b_dot_r), a_dot_r);

    // project the ends of a into b coordinates
    float a_l = -x_dot_r - hlen_a;
    float a_r = -x_dot_r + hlen_a;

    // find the intersection of the two on b
    float i_l = max(a_l, -hlen_b);
    float i_r = min(a_r, hlen_b);

    if (i_l > i_r) {
      // they don't intersect
      if (a_l < -hlen_b) {
        t_a = hlen_a;
        t_b = -hlen_b;
      } else {
        t_a = -hlen_a;
        t_b = hlen_b;
      }
    } else {
      // the segments intersect, pick two points
      *two_pts = 1;
      t_b = i_l;
      t_a = t_b + x_dot_r;
      t_b2 = i_r;
      t_a2 = t_b2 + x_dot_r;
    }

    // if they weren't in the same direction, negate
    if (a_dot_b < 0.f) {
      t_b = -t_b;
      t_b2 = -t_b2;
    }
  }
  *p_a = r_a + t_a*a;
  *p_b = r_b + t_b*b;
  if (two_pts) {
    *p_a2 = r_a + t_a2*a;
    *p_b2 = r_b + t_b2*b;
  }
}


// intest distance from a point p to a plane passing through point v
// with unit normal n
float pt_to_plane_dist(float4 v, float4 n, float4 p)
{
  return dot((p-v), n);
}


__kernel void find_plane_contacts(const int max_cells,
                                  const int max_contacts,
                                  const int n_planes,
                                  __global const float4* plane_pts,
                                  __global const float4* plane_norms,
                                  __global const float* plane_coeffs,
                                  __global const float4* centers,
                                  __global const float4* dirs,
                                  __global const float* lens,
                                  __global const float* rads,
                                  __global int* n_cts,
                                  __global int* frs,
                                  __global int* tos,
                                  __global float* dists,
                                  __global float4* pts,
                                  __global float4* norms,
                                  __global float* reldists,
                                  __global float* stiff)
{
  int i = get_global_id(0);

  // collision count
  int k = n_cts[i]; //keep existing contacts

  float4 end1 = centers[i] - 0.5f*lens[i]*dirs[i]; // 'left' end of the cell
  float4 end2 = centers[i] + 0.5f*lens[i]*dirs[i]; // 'right' end of the cell

  for (int n = 0; n < n_planes; n++) { // loop through all planes
    int to1 = -2*n - 1; // 'to' if left end has contact with plane n
    int to2 = to1 - 1;  // 'to' if right end has contact with plane n

    float dist1 = pt_to_plane_dist(plane_pts[n], plane_norms[n], end1)-rads[i];
    float dist2 = pt_to_plane_dist(plane_pts[n], plane_norms[n], end2)-rads[i];

    // check for old contacts with this plane
    int cti1 = -1;
    int cti2 = -1;
    for (int m = i*max_contacts; m < i*max_contacts+n_cts[i]; m++) {
      if (tos[m] == to1) cti1 = m;
      else if (tos[m] == to2) cti2 = m;
    }

    // common to both ends
    float4 norm = -plane_norms[n];

    bool two_pts = ((cti1 >= 0) || (dist1<0.f) ) && ( (cti2 >= 0) || (dist2<0.f) );
    float stiffness = two_pts*ISQRT2 + (!two_pts)*1.0;

    // if we're in contact, or we were in contact, recompute
    if ((cti1 >= 0) || (dist1<0.f) ){
      // need to make a new contact
      if (cti1 < 0) {
        cti1 = i*max_contacts+k;
        k++;
      }

      frs[cti1] = i;
      tos[cti1] = to1;
      dists[cti1] = dist1;
      pts[cti1] = end1; // FIXME: not really the right point
      norms[cti1] = norm;
      reldists[cti1] = stiffness*plane_coeffs[n]*dist1;
      stiff[cti1] = stiffness*plane_coeffs[n];
    }

    if ( (cti2 >= 0) || (dist2<0.f) ){
      if (cti2 < 0) {
        cti2 = i*max_contacts+k;
        k++;
      }

      frs[cti2] = i;
      tos[cti2] = to2;
      dists[cti2] = dist2;
      pts[cti2] = end2;
      norms[cti2] = norm;
      reldists[cti2] = stiffness*plane_coeffs[n]*dist2;
      stiff[cti2] = stiffness*plane_coeffs[n];
    }
  }
  n_cts[i] = k;
}


//
__kernel void find_contacts(const int max_cells,
                            const int n_cells,
                            const int grid_x_min,
                            const int grid_x_max,
                            const int grid_y_min,
                            const int grid_y_max,
                            const int n_sqs,
                            const int max_contacts,
                            __global const float4* centers,
                            __global const float4* dirs,
                            __global const float* lens,
                            __global const float* rads,
                            __global const int* sqs,
                            __global const int* sorted_ids,
                            __global const int* sq_inds,
                            __global int* n_cts,
                            __global int* frs,
                            __global int* tos,
                            __global float* dists,
                            __global float4* pts,
                            __global float4* norms,
                            __global float* reldists,
                            __global float* stiff,
                            __global float* overlap)
{
  // our id
  int i = get_global_id(0);

  // collision count
  int k = n_cts[i]; //keep existing contacts

  // what square are we in?
  int grid_x_range = grid_x_max-grid_x_min;
  int grid_y_range = grid_y_max-grid_y_min;
  int sq_row = sqs[i] / grid_x_range; // square row
  int sq_col = sqs[i] % grid_x_range; // square col

  // loop through our square and the eight squares surrounding it
  // (fewer squares if we're on an edge)
  for (int row = max(0, sq_row-1); row < min((int)(sq_row+2), grid_y_range); row++) {
    for (int col = max(0, sq_col-1); col < min((int)(sq_col+2), grid_x_range); col++) {

      // what square is this?
      int sq = row*grid_x_range + col;

      // loop through all the cell ids in the current square (but
      // don't go past the end of the list)
      for (int n = sq_inds[sq]; n < (sq < n_sqs-1 ? sq_inds[sq+1] : n_cells); n++) {

        int j = sorted_ids[n]; // the neighboring cell

        if (j<=i) continue; // we can't collide with ourself, only find low -> hi contacts

        // Look for any existing contacts, count them, and store idx's
        int n_existing_cts=0;
        int existing_cts_idx[2];
        for (int m=i*max_contacts; m<i*max_contacts+n_cts[i]; m++)
        {
          if (tos[m]==j)
          {
             existing_cts_idx[n_existing_cts++] = m;
          }
        }

        // Are we within possible touching distance, or is there an existing contact to update?
        if (length(centers[i]-centers[j]) > 0.5*(lens[i]+lens[j])+rads[i]+rads[j]+MARGIN && n_existing_cts==0)
        {
          // if not...
          continue;
        }

        // if so, find closest points
        float4 pi, pj; // pi is point on our line seg, pj point on other
        float4 pi2, pj2; // optional second points
        int two_pts = 0; // are there two contacts?
        closest_points_on_segments(centers[i], centers[j],
                                   dirs[i], dirs[j],
                                   lens[i], lens[j],
                                   &pi, &pj, &pi2, &pj2, &two_pts);

        float4 v_ij = pj-pi; // vector between closest points
        float dist = length(v_ij) - (rads[i]+rads[j]);
        float4 norm = normalize(v_ij); // normal on our cell at the contact
        float4 pt = pi + rads[i]*norm; // point on the capsule surface
          
        int ct_i;
        ct_i = i*max_contacts+k; // index of this contact

        float stiffness = two_pts*ISQRT2 + (!two_pts)*1.0;
        
        float overlap_length = two_pts*fabs(dot(pi-pi2,dirs[i]))/(2*POINT) + (!two_pts)*1.0;
          
        if (dist < MARGIN)
        {
            if (n_existing_cts==0)
            {
              // make new contact and compute distance etc.
              k++; // next contact
              frs[ct_i] = i;
              tos[ct_i] = j;
              dists[ct_i] = dist;
              pts[ct_i] = pt;
              norms[ct_i] = norm;
              reldists[ct_i] = stiffness*RHS_FRAC*dist;
              stiff[ct_i] = stiffness;
              overlap[ct_i] = overlap_length;
            }
	}
        if(n_existing_cts>0){
          // recompute dist etc. for existing contact
          int idx = existing_cts_idx[0];
          dists[idx] = dist;
          pts[idx] = pt;
          norms[idx] = norm;
          reldists[idx] = stiffness*RHS_FRAC*dist;
          stiff[idx] = stiffness;
          overlap[idx]=overlap_length;
        }


        if (!two_pts){
	  if(n_existing_cts>1){
	    // Not parallel, but were before - how to deal with this?
	    // Set stiffness and rhs (reldists) to zero so that this row has no effect
	    int idx = existing_cts_idx[1];
	    /*dists[idx] = 0.0;
	    pts[idx] = pt;
	    norms[idx] = norm;*/
	    reldists[idx] = 0.0;
            stiff[idx] = 0.0;
          overlap[idx]=overlap_length;
	  }
	  continue;
	}

        // if we had two contacts, add the second point
        v_ij = pj2-pi2;
        dist = length(v_ij) - (rads[i]+rads[j]);
        norm = normalize(v_ij);
        pt = pi2 + rads[i]*norm;
          
	// Are cells moving together or penetrating?
	if (dist < MARGIN)
	{
            if (n_existing_cts<2)
            {
              ct_i = i*max_contacts+k;
              k++;
              frs[ct_i] = i;
              tos[ct_i] = j;
              dists[ct_i] = dist;
              pts[ct_i] = pt;
              norms[ct_i] = norm;
              reldists[ct_i] = stiffness*RHS_FRAC*dist;
              stiff[ct_i] = stiffness;
              overlap[ct_i]=overlap_length;
            }
	}
        if(n_existing_cts>1){
          // recompute dist etc. for existing contact
          int idx = existing_cts_idx[1];
          dists[idx] = dist;
          pts[idx] = pt;
          norms[idx] = norm;
          reldists[idx] = stiffness*RHS_FRAC*dist;
          stiff[idx] = stiffness;
          overlap[idx]=overlap_length;
        }

      }
    }
  }
  n_cts[i] = k;

  // zero out unused contacts
  // this IS necessary for calculate_Mx to work
  for (int u = k; u < max_contacts; u++) {
    frs[i*max_contacts+u] = 0;
    tos[i*max_contacts+u] = 0;
  }
}



__kernel void collect_tos(const int max_cells,
                          const int n_cells,
                          const int grid_x_min,
                          const int grid_x_max,
                          const int grid_y_min,
                          const int grid_y_max,
                          const int n_sqs,
                          const int max_contacts,
                          __global const int* sqs,
                          __global const int* sorted_ids,
                          __global const int* sq_inds,
                          __global const int* n_cts,
                          __global const int* frs,
                          __global const int* tos,
                          __global int* cell_tos,
                          __global int* n_cell_tos)
{
  // our id
  int i = get_global_id(0);

  // what square are we in?
  int grid_x_range = grid_x_max-grid_x_min;
  int grid_y_range = grid_y_max-grid_y_min;
  int sq_row = sqs[i] / grid_x_range; // square row
  int sq_col = sqs[i] % grid_x_range; // square col

  int k = 0; // how many cts are we the 'to' of?

  // loop through our square and the eight squares surrounding it
  // (fewer squares if we're on an edge)
  for (int row = max(0, sq_row-1); row < min((int)(sq_row+2), grid_y_range); row++) {
    for (int col = max(0, sq_col-1); col < min((int)(sq_col+2), grid_x_range); col++) {

      // what square is this?
      int sq = row*grid_x_range + col;

      // loop through all the cell ids in the current square (but
      // don't go past the end of the list)
      for (int n = sq_inds[sq]; n < (sq < n_sqs-1 ? sq_inds[sq+1] : n_cells); n++) {

        int j = sorted_ids[n]; // the neighboring cell

        if (j>=i || j<0) continue; // only looking at cells with lower ids

        // look through all the other cell's contacts
        for (int m = 0; m < n_cts[j]; m++) {
          int ct_i = (int)(j)*(int)(max_contacts)+m;
          if (tos[ct_i] == i) {
            // if one of them is us, collect it
            cell_tos[i*max_contacts+k] = ct_i;
            k++;
          }
        }
      }
    }
  }
  n_cell_tos[i] = k;

  // zero out unused tos
  // this probably isn't necessary
  for (int u = k; u < max_contacts; u++) {
    cell_tos[i*max_contacts+u] = 0;
  }

}


__kernel void build_matrix(const int max_contacts,
                           const float muA,
                           const float gamma,
                           __global const float4* centers,
                           __global const float4* dirs,
                           __global const float* lens,
                           __global const float* rads,
                           __global const int* n_cts,
                           __global const int* frs,
                           __global const int* tos,
                           __global const float4* pts,
                           __global const float4* norms,
                           __global float8* fr_ents,
                           __global float8* to_ents,
                           __global float* stiff)
{
  int id = get_global_id(0);
  int ct = get_global_id(1);

  if (ct >= n_cts[id]) return;

  int i = id*max_contacts + ct;

  int a = frs[i];
  float4 r_a = pts[i]-centers[a];
  float8 fr_ent = 0.f;

  float4 Ia[4];
  cyl_inv_inertia_tensor(muA, lens[a]+2.f*rads[a], dirs[a], Ia);
  float4 nxr_a = 0.f;
  nxr_a = cross(norms[i], r_a);

  fr_ent.s012 = norms[i].s012/(muA*(lens[a]+2.f*rads[a]));
  fr_ent.s3 = -dot(nxr_a, Ia[0]);
  fr_ent.s4 = -dot(nxr_a, Ia[1]);
  fr_ent.s5 = -dot(nxr_a, Ia[2]);
  fr_ent.s6 = (1.f/gamma) * dot(dirs[a], r_a) * dot(dirs[a], norms[i])/(lens[a]+2.f*rads[a]);
  fr_ents[i] = fr_ent * stiff[i];

  int b = tos[i];

  // plane contacts have no to_ent, and have negative indices
  if (b < 0) {
    to_ents[i] = 0.f;
    return;
  }

  float4 r_b = pts[i]-centers[b];
  float8 to_ent = 0.f;

  float4 Ib[4];
  cyl_inv_inertia_tensor(muA, lens[b]+2.f*rads[b], dirs[b], Ib);
  float4 nxr_b = 0.f;
  nxr_b = cross(norms[i], r_b);

  to_ent.s012 = norms[i].s012/(muA*(lens[b]+2.f*rads[b]));
  to_ent.s3 = -dot(nxr_b, Ib[0]);
  to_ent.s4 = -dot(nxr_b, Ib[1]);
  to_ent.s5 = -dot(nxr_b, Ib[2]);
  to_ent.s6 = (1.f/gamma) * dot(dirs[b], r_b) * dot(dirs[b], norms[i])/(lens[b]+2.f*rads[b]);
  to_ents[i] = to_ent * stiff[i];
}

__kernel void calculate_Mx(const int max_contacts,
                           __global const int* frs,
                           __global const int* tos,
                           __global const float8* fr_ents,
                           __global const float8* to_ents,
                           __global const float8* deltap,
                           __global float* Mx)
{
  int id = get_global_id(0);
  int ct = get_global_id(1);
  int i = id*max_contacts + ct;
  int a = frs[i];
  int b = tos[i];
  if (a == 0 && b == 0) return; // not a contact
  float8 to_ents_i = b < 0 ? 0.f : to_ents[i];
  //my machine can't dot float8s...
  float res0123 = dot(fr_ents[i].s0123, deltap[a].s0123) - dot(to_ents_i.s0123, deltap[b].s0123);
  float res4567 = dot(fr_ents[i].s4567, deltap[a].s4567) - dot(to_ents_i.s4567, deltap[b].s4567);
  Mx[i] = res0123 + res4567;
}


__kernel void calculate_MTMx(const int max_contacts,
                             __global const int* n_cts,
                             __global const int* n_cell_tos,
                             __global const int* cell_tos,
                             __global const float8* fr_ents,
                             __global const float8* to_ents,
                             __global const float* Mx,
                             __global float8* MTMx)
{
  int i = get_global_id(0);
  int base = i*max_contacts;
  float8 res = 0.f;
  for (int k = base; k < base+n_cts[i]; k++) {
    float8 oldres = res;
    res += fr_ents[k]*Mx[k];
  }
  for (int k = base; k < base+n_cell_tos[i]; k++) {
    int n = cell_tos[k];
    if (n < 0) continue;
    res -= to_ents[n]*Mx[n];
  }
  MTMx[i] = res;
}

__kernel void calculate_Minv_x(const float muA,
			       const float gamma,
			       __global const float4* dirs,
			       __global const float* lens,
			       __global const float* rads,
			       __global const float8* x,
			       __global float8* Minvx)
{
  int i = get_global_id(0);

  float8 xi = x[i];
  float8 v = 0.f;
  v.s012 = xi.s012/(muA*(lens[i]+2.f*rads[i]));

  float4 Iinv[4];
  cyl_inv_inertia_tensor(muA, lens[i]+2.f*rads[i], dirs[i], Iinv);
  float4 L = 0.f;
  L.s012 = xi.s345;
  float4 w = matmul(Iinv, L);
  v.s345 = w.s012;

  v.s6 = xi.s6/gamma;

  Minvx[i] = v;
}


__kernel void predict(__global const float4* centers,
                        __global const float4* dirs,
                        __global const float* lens,
                        __global const float4* dcenters,
                        __global const float4* dangs,
                        __global const float* dlens,
                        __global float4* pred_centers,
                        __global float4* pred_dirs,
                        __global float* pred_lens)
{
  int i = get_global_id(0);
  float4 center_i = centers[i];
  float4 dir_i = dirs[i];
  float len_i = lens[i];
  float4 dcenter_i = dcenters[i];
  float4 dang_i = dangs[i];
  float dlen_i = dlens[i];

  pred_centers[i] = center_i + dcenter_i;

  if (length(dang_i)>1e-12)
  {
    float4 rot_axis = normalize(dang_i);
    float rot_angle = length(dang_i);
    pred_dirs[i] = normalize(rot(rot_axis, rot_angle, dir_i));
  } else {
    pred_dirs[i] = dir_i;
  }

  pred_lens[i] = len_i + dlen_i;

}

__kernel void integrate(__global float4* centers,
                        __global float4* dirs,
                        __global float* lens,
                        __global float4* dcenters,
                        __global float4* dangs,
                        __global float* dlens)
{
  int i = get_global_id(0);
  float4 center_i = centers[i];
  float4 dir_i = dirs[i];
  float len_i = lens[i];
  float4 dcenter_i = dcenters[i];
  float4 dang_i = dangs[i];
  float dlen_i = dlens[i];

  centers[i] = center_i + dcenter_i;

  if (length(dang_i)>1e-12)
  {
    float4 rot_axis = normalize(dang_i);
    float rot_angle = length(dang_i);
    dirs[i] = normalize(rot(rot_axis, rot_angle, dir_i));
  }

  lens[i] = len_i + dlen_i;

  // fully damped
  dcenters[i] = 0.f;
  dangs[i] = 0.f;
  dlens[i] = 0.f;
}

__kernel void add_impulse(const float muA,
                          const float gamma,
                          __global const float8* deltap,
                          __global const float4* dirs,
                          __global const float* lens,
                          __global const float* rads,
                          __global float4* dcenters,
                          __global float4* dangs,
                          __global const float* target_dlens,
                          __global float* dlens)
{
  int i = get_global_id(0);

  float4 dir_i = dirs[i];
  float len_i = lens[i];
  float rad_i = rads[i];
  float8 deltap_i = deltap[i];

  float4 dplin = 0.f;
  dplin.s012 = deltap_i.s012;
  dcenters[i] += dplin/(muA*(lens[i]+2.f*rads[i]));

  // FIXME: should probably store these so we don't recompute them
  float4 Iinv[4];
  cyl_inv_inertia_tensor(muA, len_i+2.f*rad_i, dir_i, Iinv);
  float4 dL = 0.f;
  dL.s012 = deltap_i.s345;
  float4 dpang = matmul(Iinv, dL);
  float dpangmag = length(dpang);
  if (dpangmag>ANG_LIMIT)
  {
    dpang *= ANG_LIMIT/dpangmag;
  }
  dangs[i] += dpang;
  
  float dplen = deltap_i.s6;
  //dlens[i] += max(0.f, target_dlens[i] + dplen/gamma);
  dlens[i] = max(0.f, dlens[i]+dplen/gamma);
}



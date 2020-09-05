#define MARGIN 0.0f

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

// intest distance from a point p to a plane passing through point v
// with unit normal n
float pt_to_plane_dist(float4 v, float4 n, float4 p)
{
  return dot((p-v), n);
}

// intest distance from a point p to a sphere at point v
// with radius r
float pt_to_sphere_dist(float4 v, float rad, float4 p)
{
  return length(p - v) - rad;
}

// intest distance from a point p to a sphere at point v
// with radius r
float pt_to_pt_dist(float4 p1, float4 p2)
{
  return length(p1 - p2);
}

__kernel void find_sphere_contacts(const int max_cells,
                                  const int max_contacts,
                                  const int n_spheres,
                                  __global const float4* sphere_pts,
                                  __global const float* sphere_coeffs,
                                  __global const float* sphere_rads,
                                  __global const float* sphere_norms,
                                  __global const float4* centers,
                                  __global const float4* dirs,
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

  for (int n = 0; n < n_spheres; n++) { // loop through all spheres
    int to = -n - 1; // 'to' if left end has contact with sphere n

    float dist = sphere_norms[n] * pt_to_sphere_dist(sphere_pts[n], sphere_rads[n], centers[i])-rads[i];

    // check for old contacts with this sphere
    int cti = -1;
    for (int m = i*max_contacts; m < i*max_contacts+n_cts[i]; m++) {
      if (tos[m] == to) cti = m;
    }

    float stiffness = 1.0;

    // if we're in contact, or we were in contact, recompute
    if ((cti >= 0) || (dist<0.f) ){
      // need to make a new contact
      if (cti < 0) {
        cti = i*max_contacts+k;
        k++;
      }

      frs[cti] = i;
      tos[cti] = to;
      dists[cti] = dist;
      pts[cti] = centers[i]; // FIXME: not really the right point
      norms[cti] = sphere_norms[n] *  normalize(-centers[i] + sphere_pts[n]);
      reldists[cti] = stiffness*sphere_coeffs[n]*dist;
      stiff[cti] = stiffness*sphere_coeffs[n];
    }

  }
  n_cts[i] = k;
}


__kernel void find_plane_contacts(const int max_cells,
                                  const int max_contacts,
                                  const int n_planes,
                                  __global const float4* plane_pts,
                                  __global const float4* plane_norms,
                                  __global const float* plane_coeffs,
                                  __global const float4* centers,
                                  __global const float4* dirs,
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

  for (int n = 0; n < n_planes; n++) { // loop through all planes
    int to = -n - 1; // 'to' if left end has contact with plane n

    float dist = pt_to_plane_dist(plane_pts[n], plane_norms[n], centers[i])-rads[i];

    // check for old contacts with this plane
    int cti = -1;
    for (int m = i*max_contacts; m < i*max_contacts+n_cts[i]; m++) {
      if (tos[m] == to) cti = m;
    }

    // common to both ends
    float4 norm = -plane_norms[n];

    float stiffness = 1.f; 

    // if we're in contact, or we were in contact, recompute
    if ((cti >= 0) || (dist<0.f) ){
      // need to make a new contact
      if (cti < 0) {
        cti = i*max_contacts+k;
        k++;
      }

      frs[cti] = i;
      tos[cti] = to;
      dists[cti] = dist;
      pts[cti] = centers[i]; // FIXME: not really the right point
      norms[cti] = norm;
      reldists[cti] = stiffness*plane_coeffs[n]*dist;
      stiff[cti] = stiffness*plane_coeffs[n];
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
                            const float Ws,
                            const float Wc,
                            __global const float4* centers,
                            __global const float4* dirs,
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
                            __global float* overlap,
			    __global float4* avg_neighbour_dir)
{
  // our id
  int i = get_global_id(0);

  float4 sum_nbr_dir_i = 0.f;
  avg_neighbour_dir[i] = 0.f;

  // collision count
  int k = n_cts[i]; //keep existing contacts
  int num_neighbours = 0; // count number of neighbours (note: not same as contacts since contacts are low-high index

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
        float dist = length(centers[i]-centers[j]) - rads[i]-rads[j];

	int n_existing_cts=0;
	int existing_cts_idx[2];
	if (j>i)
	{
		// Look for any existing contacts, count them, and store idx's
		for (int m=i*max_contacts; m<i*max_contacts+n_cts[i]; m++)
		{
		  if (tos[m]==j)
		  {
		     existing_cts_idx[n_existing_cts++] = m;
		  }
		}
	}

	// Compute sum to find average neighbour direction
	if (dist < MARGIN || n_existing_cts>0)
	{
		num_neighbours++;
		sum_nbr_dir_i = sum_nbr_dir_i + centers[j] - centers[i];
	}

        if (j<=i) continue; // we can't collide with ourself, only find low -> hi contacts


	// This is the model from Smeets et al. (2016)
	float dij = length(centers[i]-centers[j]);
	const float R = 1.f;
	float dist_adh = -(2.f/R) * ( Ws - (Ws + Wc) * (dij - R) / R );
	float4 pt = 0.5f * (centers[i]+centers[j]);
	float4 norm = normalize(centers[j]-centers[i]);

        int ct_i;
        ct_i = i*max_contacts+k; // index of this contact

        float stiffness = 1.f;
        
        if (dist < MARGIN || n_existing_cts>0)
        {
	    if (n_existing_cts==0)
            {
              // make new contact and compute distance etc.
              k++; // next contact
              frs[ct_i] = i;
              tos[ct_i] = j;
              dists[ct_i] = dist_adh;
              pts[ct_i] = pt;
              norms[ct_i] = norm;
              reldists[ct_i] = stiffness*dist_adh;
              stiff[ct_i] = stiffness;
            }
	if(n_existing_cts>0){
	  // recompute dist etc. for existing contact
	  int idx = existing_cts_idx[0];
	  dists[idx] = dist;
	  pts[idx] = pt;
	  norms[idx] = norm;
	  reldists[idx] = stiffness*dist_adh;
	  stiff[idx] = stiffness;
	}

	}
	}
     }
  }
  n_cts[i] = k;
  avg_neighbour_dir[i] = normalize( sum_nbr_dir_i / (float) num_neighbours );

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
                           __global const float4* centers,
                           __global const float4* dirs,
                           __global const float* rads,
                           __global const int* n_cts,
                           __global const int* frs,
                           __global const int* tos,
                           __global const float4* pts,
                           __global const float4* norms,
                           __global float4* fr_ents,
                           __global float4* to_ents,
                           __global float* stiff)
{
  int id = get_global_id(0);
  int ct = get_global_id(1);

  if (ct >= n_cts[id]) return;

  int i = id*max_contacts + ct;

  float4 fr_ent = 0.f;

  fr_ent.s012 = norms[i].s012;
  fr_ents[i] = fr_ent * stiff[i];

  int b = tos[i];

  // plane and sphere contacts have no to_ent, and have negative indices
  if (b < 0) {
    to_ents[i] = 0.f;
    return;
  }

  float4 to_ent = 0.f;
  to_ent.s012 = norms[i].s012;
  to_ents[i] = to_ent * stiff[i];
}

__kernel void calculate_Bx(const int max_contacts,
                           __global const int* frs,
                           __global const int* tos,
                           __global const float4* fr_ents,
                           __global const float4* to_ents,
                           __global const float4* deltap,
			   const float Ws,
			   const float Wc,
                           __global float* Bx)
{
  int id = get_global_id(0);
  int ct = get_global_id(1);
  int i = id*max_contacts + ct;
  int a = frs[i];
  int b = tos[i];
  if (a == 0 && b == 0) return; // not a contact
  float4 to_ents_i = b < 0 ? 0.f : to_ents[i];

  const float R = 1.f;

  float res = dot(fr_ents[i], deltap[a]) - dot(to_ents_i, deltap[b]);
  const float fac = 2.F * (Ws + Wc) / (R*R);
  if (b>=0)
  {
    // Not a plane or sphere contact
    res *= fac;
  }
  Bx[i] = res;
}


__kernel void calculate_BTBx(const int max_contacts,
                             __global const int* n_cts,
                             __global const int* n_cell_tos,
                             __global const int* cell_tos,
                             __global const float4* fr_ents,
                             __global const float4* to_ents,
                             __global const float* Bx,
                             __global float4* BTBx)
{
  int i = get_global_id(0);
  int base = i*max_contacts;
  float4 res = 0.f;
  for (int k = base; k < base+n_cts[i]; k++) {
    float4 oldres = res;
    res += fr_ents[k]*Bx[k];
  }
  for (int k = base; k < base+n_cell_tos[i]; k++) {
    int n = cell_tos[k];
    if (n < 0) continue;
    res -= to_ents[n]*Bx[n];
  }
  BTBx[i] = res;
}

__kernel void calculate_Mx(const float gamma_s,
			       __global const float4* dirs,
			       __global const float* rads,
			       __global const float4* x,
			       __global float4* Mx)
{
  int i = get_global_id(0);

  float4 xi = x[i];
  float4 v = 0.f;
  v.s012 = xi.s012 * gamma_s;

  Mx[i] = v;
}


__kernel void predict(__global const float4* centers,
                        __global const float4* dirs,
                        __global const float4* dcenters,
                        __global float4* pred_centers)
{
  int i = get_global_id(0);
  float4 center_i = centers[i];
  float4 dir_i = dirs[i];
  float4 dcenter_i = dcenters[i];

  pred_centers[i] = center_i + dcenter_i;


}

__kernel void add_impulse(const float gamma_s,
                          __global const float4* deltap,
                          __global const float4* dirs,
                          __global const float* rads,
                          __global float4* dcenters)
{
  int i = get_global_id(0);

  float4 dir_i = dirs[i];
  float rad_i = rads[i];
  float4 deltap_i = deltap[i];

  float4 dplin = 0.f;
  dplin.s012 = deltap_i.s012;
  dcenters[i] += dplin;
}


float angle_between_vectors(float4 vec1, float4 vec2, float4 normal)
{
    // Get the unsigned angle as arccos(dot product)
    float ang = acos( clamp( dot(vec1, vec2), -1.f, 1.f ) );

    // Find the sign of the angle from cross-product with normal
    float4 cross_vec = cross(vec1, vec2);
    float ang_sign = sign(dot(cross_vec, normal));

    return ang * ang_sign;
}

// This kernel integrates the system while at each step rotating the orientation 
// vector onto the tangent to a sphere at (0,0) (if spherical flag is set)
__kernel void integrate(__global float4* centers,
                        __global float4* dirs,
                        __global float4* dcenters,
			__global float4* avg_neighbour_dir,
			__global float4* signal_gradient,
			__global float* noise,
			const float fcil,
			const float ftax,
			const float forg,
			const float D,
			const float dt,
			const int spherical,
			const float sphere_radius,
			const int numSignals)
{
  int i = get_global_id(0);
  float4 center_i = centers[i];
  float4 dir_i = dirs[i];
  float4 dcenter_i = dcenters[i];
  float noise_i = noise[i];
  float4 avg_neighbour_dir_i = avg_neighbour_dir[i];
  float4 signal_gradient_i;
  if (signal_gradient)
  {
	int sigbase = i*numSignals;
  	signal_gradient_i = signal_gradient[sigbase];
  } else
  {
  	signal_gradient_i = 0.f;
  }
  const float4 org_center = {sphere_radius, 0.f, 0.f, 0.f};
  float4 org_center_dir = normalize(org_center - center_i);
  float4 normal = {0.f, 0.f, 1.f, 0.f};
  if (spherical==1) 
  {
	// Find new coordinate system tangent to sphere surface
	normal = normalize(center_i);
	float4 new_y_axis = cross(dir_i, normal);
	
	// Rotate polarity vector onto tangent plane
	dir_i = normalize(cross(normal, new_y_axis));
	
	// Project taxis directions onto tangent plane
	signal_gradient_i = signal_gradient_i - dot(signal_gradient_i, normal)*normal;
	org_center_dir = org_center_dir - dot(org_center_dir, normal)*normal;
	org_center_dir = normalize(org_center_dir);
  } 

  // Rotate by change in angle around z-axis
  float angle = sqrt(2 * D) * noise_i;
  if (length(avg_neighbour_dir_i)>0.f)
  {
	angle += fcil * angle_between_vectors(dir_i, -avg_neighbour_dir_i, normal);
  }
  if (length(signal_gradient_i)>0.f)
  {
	angle += ftax * length(signal_gradient_i) * angle_between_vectors(dir_i, normalize(signal_gradient_i), normal);
  }
  if (length(org_center_dir)>0.f)
  {
	angle += forg * angle_between_vectors(dir_i, org_center_dir, normal);
  }
  dir_i = rot(normal, dt * angle, dir_i);
  dirs[i] = normalize(dir_i);

  centers[i] = center_i + dcenter_i;

  // fully damped
  dcenters[i] = 0.f;
}


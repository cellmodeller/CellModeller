import numpy as np

EPSILON = 0.01

# Multiply a 3x3 matrix and a 3-vector
def matmul(m: np.ndarray, v: np.ndarray) -> np.ndarray:
    # m: shape (3, 3), v: shape (3,)
    return np.dot(m, v)

# Transpose a 3x3 matrix
def transpose(m: np.ndarray) -> np.ndarray:
    # m: shape (3, 3)
    return m.T

# Multiply two 3x3 matrices
def matmulmat(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    # a, b: shape (3, 3)
    return np.dot(a, b)

# Inverse of a unit quaternion (3-vector imaginary part, 1 real part)
def quat_inv(q: np.ndarray) -> np.ndarray:
    # q: shape (4,) -> [x, y, z, w]
    l2 = np.dot(q, q)
    if l2 == 0.0:
        return q
    return np.array([-q[0], -q[1], -q[2], q[3]]) / l2

# Quaternion multiplication
def quat_prod(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    # a, b: shape (4,) -> [x, y, z, w]
    x =  a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0]
    y = -a[0]*b[2] + a[1]*b[3] + a[2]*b[0] + a[3]*b[1]
    z =  a[0]*b[1] - a[1]*b[0] + a[2]*b[3] + a[3]*b[2]
    w = -a[0]*b[0] - a[1]*b[1] - a[2]*b[2] + a[3]*b[3]
    return np.array([x, y, z, w])

# Rotate a 3-vector using a quaternion
def quat_rot(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    # q: shape (4,), v: shape (3,)
    vq = np.array([v[0], v[1], v[2], 0.0])
    qi = quat_inv(q)
    v_prime = quat_prod(quat_prod(q, vq), qi)
    return v_prime[:3]

# Rotate a vector about an axis by an angle (radians)
def rot(axis: np.ndarray, angle: float, v: np.ndarray) -> np.ndarray:
    # axis: shape (3,), v: shape (3,)
    axis = axis / np.linalg.norm(axis)  # ensure unit vector
    s = np.sin(angle / 2.0)
    q = np.append(axis * s, np.cos(angle / 2.0))  # [x, y, z, w]
    return quat_rot(q, v)

def normalize(v: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(v)
    return v / norm if norm > 0 else v

def cross(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.cross(a, b)

def cyl_inertia_tensor(muA: float, l: float, axis: np.ndarray) -> np.ndarray:
    # axis: (3,)
    diag = (muA * l**3) / 12.0

    x_axis = np.array([1.0, 0.0, 0.0])
    y_axis = np.array([0.0, 1.0, 0.0])
    z_axis = np.array([0.0, 0.0, 1.0])

    axis = normalize(axis)
    rot_ang = np.arccos(np.clip(np.dot(x_axis, axis), -1.0, 1.0))

    y_prime = y_axis.copy()
    z_prime = z_axis.copy()

    if rot_ang > EPSILON:
        y_prime = rot(z_prime, rot_ang, y_axis)
        z_prime = normalize(cross(x_axis, axis))

    # M = matrix that aligns x-axis to given axis
    M = np.vstack([axis, y_prime, z_prime])  # (3, 3)
    MT = transpose(M)

    # I: inertia tensor in local frame (x-axis principal)
    I = np.zeros((3, 3))
    I[1, 1] = diag
    I[2, 2] = diag

    # I_world = MT * I * M
    IM = matmulmat(I, M)
    I_world = matmulmat(MT, IM)

    return I_world

def cyl_inv_inertia_tensor(muA: float, l: float, axis: np.ndarray) -> np.ndarray:
    # axis: (3,)
    diag = 12.0 / (muA * l**3)

    x_axis = np.array([1.0, 0.0, 0.0])
    y_axis = np.array([0.0, 1.0, 0.0])
    z_axis = np.array([0.0, 0.0, 1.0])

    axis = normalize(axis)
    rot_ang = np.arccos(np.clip(np.dot(x_axis, axis), -1.0, 1.0))

    y_prime = y_axis.copy()
    z_prime = z_axis.copy()

    if rot_ang > EPSILON:
        y_prime = rot(z_axis, rot_ang, y_axis)
        z_prime = normalize(cross(x_axis, axis))

    M = np.vstack([axis, y_prime, z_prime])  # (3, 3)

    # Diagonal inverse inertia values only on y, z axes
    MD = np.zeros((3, 3))
    MD[1, :] = M[1, :] * diag
    MD[2, :] = M[2, :] * diag

    MDT = transpose(MD)

    I_inv_world = matmulmat(M, MDT)
    return I_inv_world

def closest_points_on_segments(r_a, r_b, a, b, len_a, len_b):
    hlen_a = len_a / 2.0
    hlen_b = len_b / 2.0
    r = r_b - r_a
    a_dot_r = np.dot(a, r)
    b_dot_r = np.dot(b, r)
    a_dot_b = np.dot(a, b)
    denom = 1.0 - a_dot_b * a_dot_b

    t_a = t_b = 0.0
    t_a2 = t_b2 = 0.0
    two_pts = False

    if np.sqrt(np.abs(denom)) > EPSILON:
        # non-parallel lines
        t_a0 = (a_dot_r - b_dot_r * a_dot_b) / denom
        t_b0 = (a_dot_r * a_dot_b - b_dot_r) / denom

        on_a = abs(t_a0) < hlen_a
        on_b = abs(t_b0) < hlen_b

        if not on_a and not on_b:
            c_a = np.copysign(hlen_a, t_a0)
            c_b = np.copysign(hlen_b, t_b0)

            dd_dt_a = 2.0 * (c_a - a_dot_b * c_b - a_dot_r)
            dd_dt_b = 2.0 * (c_b - a_dot_b * c_a + b_dot_r)

            if np.sign(dd_dt_a) == np.sign(c_a):
                t_b = c_b
                t_a = np.clip(t_b * a_dot_b + a_dot_r, -hlen_a, hlen_a)
            else:
                t_a = c_a
                t_b = np.clip(t_a * a_dot_b - b_dot_r, -hlen_b, hlen_b)

        elif on_a and not on_b:
            t_b = np.copysign(hlen_b, t_b0)
            t_a = np.clip(t_b * a_dot_b + a_dot_r, -hlen_a, hlen_a)

        elif not on_a and on_b:
            t_a = np.copysign(hlen_a, t_a0)
            t_b = np.clip(t_a * a_dot_b - b_dot_r, -hlen_b, hlen_b)

        else:
            t_a = t_a0
            t_b = t_b0

    else:
        # lines are roughly parallel
        x_dot_r = np.copysign(min(abs(a_dot_r), abs(b_dot_r)), a_dot_r)

        a_l = -x_dot_r - hlen_a
        a_r = -x_dot_r + hlen_a

        i_l = max(a_l, -hlen_b)
        i_r = min(a_r, hlen_b)

        if i_l > i_r:
            if a_l < -hlen_b:
                t_a = hlen_a
                t_b = -hlen_b
            else:
                t_a = -hlen_a
                t_b = hlen_b
        else:
            # segments intersect: return both ends of intersection
            two_pts = True
            t_b = i_l
            t_a = t_b + x_dot_r
            t_b2 = i_r
            t_a2 = t_b2 + x_dot_r

        if a_dot_b < 0.0:
            t_b = -t_b
            t_b2 = -t_b2

    p_a = r_a + t_a * a
    p_b = r_b + t_b * b
    p_a2 = r_a + t_a2 * a if two_pts else None
    p_b2 = r_b + t_b2 * b if two_pts else None

    return p_a, p_b, p_a2, p_b2, two_pts

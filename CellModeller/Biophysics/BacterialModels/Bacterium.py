from CellModeller.Biophysics.BacterialModels.Geometry import *

import numpy as np
from scipy.spatial import cKDTree

DELTA = 0.01

class Bacterium:
    def __init__(self, sim, gamma, muA, sub_steps=1):
        self.sim = sim
        self.gamma = gamma
        self.muA = muA
        self.sub_steps = sub_steps

    def setRegulator(self, reg):
        self.reg = reg

    def addCell(self, cell, pos=(0,0,0), dir=(1,0,0), rad=0.5, **kwargs):
        cell.pos = np.array(pos, dtype=float)
        cell.dir = np.array(dir, dtype=float)
        cell.radius = rad
        cell.volume = cell.length
        self.compute_ends(cell)

    def set_cells(self):
        pass

    def step(self, dt):
        self.grow_cells(dt)
        self.find_neighbours()
        for _ in range(self.sub_steps):
            self.compute_contacts()
            self.compute_torques()
            self.compute_forces()
            self.compute_compression()
            self.integrate(dt / self.sub_steps)
        return True

    def grow_cells(self, dt):
        cells = self.sim.cellStates
        for cid, cell in cells.items():
            cell.length += cell.growthRate * cell.length * dt

    def find_neighbours(self):
        cells = self.sim.cellStates
        if len(cells) < 2:
            return

        # Create consistent order of cells: list of (cid, cell) pairs
        cid_list = list(cells.keys())
        pos_array = np.array([np.array(cells[cid].pos) for cid in cid_list])

        # Build KD-tree from positions
        tree = cKDTree(pos_array)

        for i, cid in enumerate(cid_list):
            radius = 5 # 4 * cells[cid].radius + 2 * cells[cid].length
            # Find neighbors within given radius
            indices = tree.query_ball_point(pos_array[i], r=radius)
            # Exclude self (i)
            indices = [j for j in indices if j != i]
            # Map indices back to cell IDs
            cells[cid].neighbours = [cid_list[j] for j in indices]

    def compute_contacts(self):
        cells = self.sim.cellStates
        for cid, cell in cells.items():
            cell.contacts = []
            for nbr_cid in cell.neighbours:
                r_a = np.array(cell.pos)
                r_b = np.array(cells[nbr_cid].pos)
                len_a = cell.length
                len_b = cells[nbr_cid].length
                rad_a = cell.radius
                rad_b = cells[nbr_cid].radius
                centre_dist = np.linalg.norm(r_a - r_b)
                if centre_dist > len_a/2 + rad_a + len_b/2 + rad_b:
                    # Cell too far away
                    continue
                a = np.array(cell.dir)
                b = np.array(cells[nbr_cid].dir)
                p_a, p_b, p_a2, p_b2, two_pts = closest_points_on_segments(r_a, r_b, a, b, len_a, len_b)
                dist = np.linalg.norm(p_b - p_a) - cell.radius - cells[nbr_cid].radius
                dist2 = np.linalg.norm(p_b2 - p_a2) - cell.radius - cells[nbr_cid].radius if two_pts else None
                normal = p_b - p_a
                normal = normalize(normal)
                if two_pts:
                    normal2 = p_b2 - p_a2
                    normal2 = normalize(normal2)
                else:
                    normal2 = None
                if dist < DELTA:
                    contact = {
                        "nbr_cid": nbr_cid, 
                        "p_a": p_a, 
                        "p_b": p_b, 
                        "p_a2": p_a2, 
                        "p_b2": p_b2, 
                        "two_pts": two_pts,
                        "dist": dist,
                        "dist2": dist2,
                        "normal": normal,
                        "normal2": normal2
                        }
                    cell.contacts.append(contact)

    def compute_torques(self):
        cells = self.sim.cellStates
        for cid, cell in cells.items():
            total_torque = np.zeros(3)
            for contact in cell.contacts:
                p_a = contact["p_a"]
                p_b = contact["p_b"]
                two_pts = contact["two_pts"]
                p_a2 = contact["p_a2"]
                p_b2 = contact["p_b2"]
                dist = contact["dist"]
                dist2 = contact["dist2"]
                normal = contact["normal"]
                normal2 = contact["normal2"]

                # Vector from contact point on cell A to B (force direction)
                force = self.gamma * normal * dist
                
                # Lever arm from center of cell to point of contact
                r = p_a - np.array(cell.pos)

                # Torque = r × F (cross product of lever arm and force)
                torque = np.cross(r, force)

                total_torque += torque

                # If two contact points exist, compute second torque
                if two_pts:
                    force2 = self.gamma * normal2 * dist2
                    r2 = p_a2 - np.array(cell.pos)
                    torque2 = np.cross(r2, force2)
                    total_torque += torque2

            # Store total torque in the cell (you can add this as an attribute or modify as needed)
            cell.torque = total_torque

    def compute_forces(self):
        cells = self.sim.cellStates
        for cid, cell in cells.items():
            total_force = np.zeros(3)
            for contact in cell.contacts:
                p_a = contact["p_a"]
                p_b = contact["p_b"]
                two_pts = contact["two_pts"]
                p_a2 = contact["p_a2"]
                p_b2 = contact["p_b2"]
                dist = contact["dist"]
                dist2 = contact["dist2"]
                normal = contact["normal"]
                normal2 = contact["normal2"]

                force = self.gamma * dist * normal
                total_force += force

                if two_pts:
                    force2 = self.gamma * dist2 * normal2
                    total_force += force2

            # Store the net force
            cell.force = total_force

    def compute_compression(self):
        cells = self.sim.cellStates
        for cid, cell in cells.items():
            dir_vec = normalize(np.array(cell.dir))
            pos = np.array(cell.pos)
            net_compression = 0.0

            for contact in cell.contacts:
                # Force at first contact point
                f1 = self.gamma * contact["dist"] * contact["normal"]
                r1 = contact["p_a"] - pos
                sign1 = np.sign(np.dot(r1, dir_vec))  # +1 if toward +dir, -1 if toward -dir
                proj1 = np.dot(f1, dir_vec) * sign1
                net_compression += proj1

                # Optional second contact point
                if contact["two_pts"]:
                    f2 = self.gamma * contact["dist2"] * contact["normal2"]
                    r2 = contact["p_a2"] - pos
                    sign2 = np.sign(np.dot(r2, dir_vec))
                    proj2 = np.dot(f2, dir_vec) * sign2
                    net_compression += proj2

            # Store net compressive force along axis
            cell.compression = net_compression

    def integrate(self, dt):
        cells = self.sim.cellStates
        for cid, cell in cells.items():
            # --- Cell compression ---
            cell.length += np.min(cell.compression * dt / self.gamma, 0)

            # --- Linear motion (viscous drag) ---
            force = getattr(cell, 'force', np.zeros(3))
            velocity = force / self.muA / cell.length
            cell.pos = np.array(cell.pos) + dt * velocity

            # --- Rotational motion ---
            torque = getattr(cell, 'torque', np.zeros(3))
            dir_vec = normalize(np.array(cell.dir))
            length = cell.length

            # Inverse inertia tensor in world coordinates
            I_inv = cyl_inv_inertia_tensor(self.muA, length, dir_vec)

            # Angular velocity: ω = I⁻¹ * τ
            ang_vel = matmul(I_inv, torque)
            theta = np.linalg.norm(ang_vel) * dt

            if theta > 1e-8:
                axis = ang_vel / np.linalg.norm(ang_vel)
                new_dir = rot(axis, theta, dir_vec)
                cell.dir = normalize(new_dir)
            self.compute_ends(cell)
            cell.volume = cell.length # * np.pi * cell.radius ** 2

    def compute_ends(self, cell):
        pa = np.array(cell.pos)
        da = np.array(cell.dir)
        cell.ends = (pa - da * cell.length * 0.5, pa + da * cell.length * 0.5)

    def divide_geometry(self, cell, daughter1, daughter2):
        #print(f"Dividing cell {cell.id} into {daughter1.id} and {daughter2.id}, length: {cell.length}")
        parent_pos = np.array(cell.pos)
        parent_len = cell.length
        parent_rad = cell.radius
        parent_dir = np.array(cell.dir)

        daughter_len = parent_len * 0.5 - parent_rad #- 0.025
        daughter_offset = daughter_len * 0.5 + parent_rad
        center_offset = parent_dir * daughter_offset

        cdir1 = parent_dir
        jitter = np.random.uniform(-0.001, 0.001, 2)
        cdir1[0:2] += jitter
        cdir1 /= np.linalg.norm(cdir1)
        cdir2 = parent_dir
        jitter = np.random.uniform(-0.001, 0.001, 2)
        cdir2[0:2] += jitter
        cdir2 /= np.linalg.norm(cdir2)

        daughter1.pos = parent_pos + center_offset
        daughter1.dir = cdir1
        daughter1.length = daughter_len
        daughter1.radius = parent_rad
        daughter1.neighbours = []
        daughter2.pos = parent_pos - center_offset
        daughter2.dir = cdir2
        daughter2.length = daughter_len
        daughter2.radius = parent_rad
        daughter2.neighbours = []

        # Update ends
        self.compute_ends(daughter1)
        self.compute_ends(daughter2)

    def divide(self, pState, d1State, d2State, *args, **kwargs):    
        # Divide the cell in each model
        self.divide_geometry(pState, d1State, d2State)
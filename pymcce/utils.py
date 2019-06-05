"""
Utility module for functions that are frequently used.
"""
import numpy as np
import mdtraj as md


class NeighborSearch(object):
    """
    Class for relatively fast queries of coordinates within a distance
    of specified coordinate.
    """
    def __init__(self, xyz, dist):
        """Initialize a NeighborSearch object by providing an array of
        coordinates and a distance threshold.

        Parameters
        ----------
        xyz : np.ndarray, float, shape=(N, 3)
            A multidmimensional array of three dimensional coordinates
        dist : float
            A distance cutoff to identify points within this distance of the
            query point.
        """
        # create an array of indices around a cubic grid
        self.neighbors = []
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                for k in (-1, 0, 1):
                    self.neighbors.append((i, j, k))
        self.neighbor_array = np.array(self.neighbors, np.int)

        self.min_ = np.min(xyz, axis=0)
        self.cell_size = np.array([dist, dist, dist], np.float)
        cell = np.array((xyz - self.min_) / self.cell_size)  # , dtype=np.int)
        # create a dictionary with keys corresponding to integer representation
        # of transformed XYZ's
        self.cells = {}
        for ix, assignment in enumerate(cell):
            # convert transformed xyz coord into integer index (so coords like
            # 1.1 or 1.9 will go to 1)
            indices = assignment.astype(int)
            # create interger indices
            t = tuple(indices)

            # NOTE: a single index can have multiple coords associated with it
            # if this integer index is already present
            if t in self.cells:
                # obtain its value (which is a list, see below)
                xyz_list, trans_coords, ix_list = self.cells[t]
                # append new xyz to xyz list associated with this entry
                xyz_list.append(xyz[ix])
                # append new transformed xyz to transformed xyz list associated
                # with this entry
                trans_coords.append(assignment)
                # append new array index
                ix_list.append(ix)
            # if this integer index is encountered for the first time
            else:
                # create a dictionary key value pair,
                # key: integer index
                # value: [[list of x,y,z], [list of transformed x,y,z], [list
                # of array indices]]
                self.cells[t] = ([xyz[ix]], [assignment], [ix])

        self.dist_squared = dist * dist

    def query_nbrs_single_point(self, point):
        """
        Given a coordinate point, return all point indexes (0-indexed) that
        are within the threshold distance from it.
        """
        cell0 = np.array((point - self.min_) / self.cell_size, dtype=np.int)
        tuple0 = tuple(cell0)
        near = []
        for index_array in tuple0 + self.neighbor_array:
            t = tuple(index_array)
            if t in self.cells:
                xyz_list, trans_xyz_list, ix_list = self.cells[t]
                for (xyz, ix) in zip(xyz_list, ix_list):
                    diff = xyz - point
                    if np.dot(diff, diff) <= self.dist_squared and float(
                            np.dot(diff, diff)) > 0.0:
                        # near.append(ix)
                        # print ix, np.dot(diff, diff)
                        near.append(ix)
        return near

    def query_point_and_distance(self, point):
        """
        Given a coordinate point, return all point indexes (0-indexed) and
        corresponding distances that are within the threshold distance from it.
        """
        cell0 = np.array((point - self.min_) / self.cell_size, dtype=np.int)
        tuple0 = tuple(cell0)
        near = []
        for index_array in tuple0 + self.neighbor_array:
            t = tuple(index_array)
            if t in self.cells:
                xyz_list, trans_xyz_list, ix_list = self.cells[t]
                for (xyz, ix) in zip(xyz_list, ix_list):
                    diff = xyz - point
                    if np.dot(diff, diff) <= self.dist_squared and float(
                            np.dot(diff, diff)) > 0.0:
                        # near.append(ix)
                        # print ix, np.dot(diff, diff)
                        near.append((ix, np.sqrt(np.dot(diff, diff))))
        return near

    def query_nbrs_multiple_points(self, points):
        """
        Given a coordinate point, return all point indexes (0-indexed) that
        are within the threshold distance from it.
        shape of points has to be (n_lig_atoms, 3)
        """
        near = []
        for point in points:
            cell0 = np.array(
                (point - self.min_) / self.cell_size, dtype=np.int)
            tuple0 = tuple(cell0)

            for index_array in tuple0 + self.neighbor_array:
                t = tuple(index_array)
                if t in self.cells:
                    xyz_list, trans_xyz_list, ix_list = self.cells[t]
                    for (xyz, ix) in zip(xyz_list, ix_list):
                        diff = xyz - point
                        if np.dot(diff, diff) <= self.dist_squared and float(
                                np.dot(diff, diff)) > 0.0:
                            # near.append(ix)
                            # print ix, np.dot(diff, diff)
                            if ix not in near:
                                near.append(ix)
        return near


def write_watpdb_from_coords(filename, coords, full_water_res=False):
    """Summary

    Parameters
    ----------
    traj : TYPE
        Description
    filename : TYPE
        Description
    water_id_list : None, optional
        Description
    wat_coords : None, optional
        Description
    full_water_res : bool, optional
        Description

    Returns
    -------
    TYPE
        Description
    """

    pdb_line_format = "{0:6}{1:>5}  {2:<3}{3:<1}{4:>3} {5:1}{6:>4}{7:1}   {8[0]:>8.3f}{8[1]:>8.3f}{8[2]:>8.3f}{9:>6.2f}{10:>6.2f}{11:>12s}\n"
    ter_line_format = "{0:3}   {1:>5}      {2:>3} {3:1}{4:4} \n"
    pdb_lines = []
    # write form the list of (water, frame) tuples
    # at_index, wat in enumerate(water_id_list):
    at = 0
    res = 0
    wat_i = 0
    with open(filename + ".pdb", 'w') as f:
        f.write("REMARK Initial number of clusters: N/A\n")
        while wat_i < len(coords):
            at_index = at  # % 10000
            res_index = res % 10000
            # wat_coords = md.utils.in_units_of(
            #    coords[wat[0], wat[1], :], "nanometers", "angstroms")
            wat_coords = coords[wat_i]
            # chain_id = possible_chains[chain_id_index]
            chain_id = "A"
            pdb_line = pdb_line_format.format(
                "ATOM", at_index, "O", " ", "WAT", chain_id, res_index, " ", wat_coords, 0.00, 0.00, "O")
            # pdb_lines.append(pdb_line)
            f.write(pdb_line)
            wat_i += 1
            if full_water_res:
                # H1_coords = md.utils.in_units_of(
                #    coords[wat[0], wat[1] + 1, :], "nanometers", "angstroms")
                H1_coords = coords[wat_i]
                pdb_line_H1 = pdb_line_format.format("ATOM", at_index + 1, "H1", " ", "WAT", chain_id, res_index, " ",
                                                     H1_coords, 0.00, 0.00, "H")
                # pdb_lines.append(pdb_line_H1)
                f.write(pdb_line_H1)
                # H2_coords = md.utils.in_units_of(
                #    coords[wat[0], wat[1] + 2, :], "nanometers", "angstroms")
                H2_coords = coords[wat_i + 1]
                pdb_line_H2 = pdb_line_format.format("ATOM", at_index + 2, "H2", " ", "WAT", chain_id, res_index, " ",
                                                     H2_coords, 0.00, 0.00, "H")
                # pdb_lines.append(pdb_line_H2)
                f.write(pdb_line_H2)
                at += 3
                res += 1
                wat_i += 2
            else:
                at += 1
                res += 1
            if res_index == 9999:
                ter_line = ter_line_format.format("TER", at, "WAT", chain_id, res_index)
                at = 1
                # pdb_lines.append(ter_line)


def initialize_grid(center, resolution, dimensions):
    """
    Parameters
    ----------
    center : TYPE
        DESCRIPTION

    """
    # set grid center, res and dimension
    # self.center = np.array(center,dtype=np.float_)
    # self.dims = np.array(dimensions)
    # self.spacing = np.array(resolution,dtype=np.float_)
    print("Creating grid ...")
    center = np.array(center, dtype=np.float_)
    dims = np.array(dimensions, dtype=np.int_)
    spacing = np.array(resolution, dtype=np.float_)
    gridmax = dims * spacing + 1.5
    # set origin
    o = center - (0.5 * dims * spacing)
    origin = np.around(o, decimals=3)
    # set grid size (in terms of total points alog each axis)
    length = np.array(dims / spacing, dtype=np.float_)
    grid_size = np.ceil((length / spacing) + 1.0)
    grid_size = np.cast['uint32'](grid_size)
    # Finally allocate the space for the grid
    grid = np.zeros(dims, dtype=np.int_)

    v_count = 0
    voxel_array = np.zeros((grid.size, 35), dtype="float64")
    # print voxel_quarts_new.shape
    for index, value in np.ndenumerate(grid):
        # point = grid.pointForIndex(index) # get cartesian coords for the
        # grid point
        _index = np.array(index, dtype=np.int32)
        # point = self.spacing * _index + self._origin
        point = _index * spacing + origin + 0.5 * spacing
        voxel_array[v_count, 1] = point[0]
        voxel_array[v_count, 2] = point[1]
        voxel_array[v_count, 3] = point[2]
        voxel_array[v_count, 0] = v_count
        # print voxel_quarts_new[v_count, 0], voxel_quarts_new[v_count, 1],
        # voxel_quarts_new[v_count, 2]
        # create a dictionary key-value pair with voxel index as key and
        # it's coords as
        # voxel_quarts[v_count].append(np.zeros(14, dtype="float64"))
        v_count += 1
    return voxel_array


def get_last_prot_at_index(pdb):
    with open(pdb, "r") as f2:
        l2 = f2.readlines()
        i = 0
        prot_last_at_index = 0
        for l in l2:
            if "EM" in l:
                prot_last_at_index = i
                break
            i += 1
        return prot_last_at_index + 1
        # return len(l2)


def get_ligand_center(ligand):
    lig = md.load_pdb(ligand, no_boxchk=True)
    com = np.zeros((lig.n_frames, 3))
    masses = np.ones(lig.n_atoms)
    masses /= masses.sum()
    com[0, :] = lig.xyz[0, :].astype('float64').T.dot(masses)
    grid_center = com[0, :]

    return grid_center


def assemble_solvated_structure(water_pdb, other_pdb, write_pdb, last_prot_index):
    with open(water_pdb, "r") as f1:
        with open(other_pdb, "r") as f2:
            with open(write_pdb, "w") as f3:
                wat_lines = f1.readlines()
                l2 = f2.readlines()
                prot_lines = l2[:last_prot_index-1]
                mem_lines = l2[last_prot_index-1:]
                #print(prot_lines[-1])
                #print(wat_lines[0])
                #print(mem_lines[0])
                #print(mem_lines[-1])
                all_lines = prot_lines + wat_lines + mem_lines
                for l in all_lines:
                    f3.write(l)

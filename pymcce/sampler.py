"""
Water sampling functions
"""
import numpy as np
import mdtraj as md


def write_watpdb_from_coords_ext(filename, coords, start_index, charge_set):
    """Writes out PDB file from the a coordinate array. Only specific to writing
    out water molecules."
    """

    pdb_line_format = "{0:6}{1:>5}{2:>3}{3:<3}{4:>3} {5:1}{6:>04}{7:1}{8[0]:>8.3f}{8[1]:>8.3f}{8[2]:>8.3f}{9:>8.3f}{10:>12.3f}{11:>16s}\n"
    ter_line_format = "{0:3}   {1:>5}      {2:>3} {3:1}{4:4} \n"
    pdb_lines = []
    at = start_index
    wat_i = 0
    with open(filename + ".pdb", 'w') as f:
        # f.write("REMARK Initial number of clusters: N/A\n")
        while wat_i < len(coords):
            res = coords[wat_i][2]
            at_index = at  # % 10000
            res_index = res % 10000
            wat_coords = coords[wat_i][0]
            conf_id = "_%03d" % coords[wat_i][1]
            chain_id = "W"
            pdb_line = pdb_line_format.format(
                "ATOM", at_index, "O", " ", "HOH", chain_id, res_index, conf_id, wat_coords, 1.600, charge_set[0],
                "01O000H011")
            f.write(pdb_line)
            wat_i += 1
            H1_coords = coords[wat_i][0]
            pdb_line_H1 = pdb_line_format.format("ATOM", at_index + 1, "1H", " ", "HOH", chain_id, res_index, conf_id,
                                                 H1_coords, 1.000, charge_set[1], "01O000H011")
            f.write(pdb_line_H1)
            H2_coords = coords[wat_i + 1][0]
            pdb_line_H2 = pdb_line_format.format("ATOM", at_index + 2, "2H", " ", "HOH", chain_id, res_index, conf_id,
                                                 H2_coords, 1.000, charge_set[2], "01O000H011")
            f.write(pdb_line_H2)
            at += 3
            res += 1
            wat_i += 2
            if res_index == 9999:
                ter_line = ter_line_format.format("TER", at, "HOH", chain_id, res_index)
                at = 1


def generate_uniform_quaternion():
    """Generate a uniform normalized quaternion 4-vector.

    References
    ----------
    [1] K. Shoemake. Uniform random rotations. In D. Kirk, editor,
    Graphics Gems III, pages 124-132. Academic, New York, 1992.
    [2] Described briefly here: http://planning.cs.uiuc.edu/node198.html

    Examples
    --------
    >>> q = MCRotationMove._generate_uniform_quaternion()

    """
    u = np.random.rand(3)
    q = np.array([np.sqrt(1 - u[0]) * np.sin(2 * np.pi * u[1]),
                  np.sqrt(1 - u[0]) * np.cos(2 * np.pi * u[1]),
                  np.sqrt(u[0]) * np.sin(2 * np.pi * u[2]),
                  np.sqrt(u[0]) * np.cos(2 * np.pi * u[2])])
    return q


def rotation_matrix_from_quaternion(q):
    """Compute a 3x3 rotation matrix from a given quaternion (4-vector).

    Parameters
    ----------
    q : 1x4 numpy.ndarray
        Quaterion (need not be normalized, zero norm OK).

    Returns
    -------
    Rq : 3x3 numpy.ndarray
        Orthogonal rotation matrix corresponding to quaternion q.

    Examples
    --------
    >>> q = np.array([0.1, 0.2, 0.3, -0.4])
    >>> Rq = MCRotationMove._rotation_matrix_from_quaternion(q)

    References
    ----------
    [1] http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion

    """

    w, x, y, z = q
    Nq = (q ** 2).sum()  # Squared norm.
    if Nq > 0.0:
        s = 2.0 / Nq
    else:
        s = 0.0

    X = x * s;
    Y = y * s;
    Z = z * s
    wX = w * X;
    wY = w * Y;
    wZ = w * Z
    xX = x * X;
    xY = x * Y;
    xZ = x * Z
    yY = y * Y;
    yZ = y * Z;
    zZ = z * Z

    Rq = np.matrix([[1.0 - (yY + zZ), xY - wZ, xZ + wY],
                    [xY + wZ, 1.0 - (xX + zZ), yZ - wX],
                    [xZ - wY, yZ + wX, 1.0 - (xX + yY)]])

    return Rq


def generate_random_rotation_matrix():
    """Return a random 3x3 rotation matrix.

    Returns
    -------
    Rq : 3x3 numpy.ndarray
        The random rotation matrix.

    """
    q = generate_uniform_quaternion()
    return rotation_matrix_from_quaternion(q)


def rotate_positions(positions, rotation_center_indices=None):
    """Return the positions after applying a random rotation to them.

    Parameters
    ----------
    positions : nx3 numpy.ndarray simtk.unit.Quantity
        The positions to rotate.
    rotation_center_indices : list
        Indices of a subset of atoms to center the rotation around, None by default.

    Returns
    -------
    rotated_positions : nx3 numpy.ndarray simtk.unit.Quantity
        The rotated positions.

    """
    #positions_unit = positions.unit
    x_initial = positions #/ positions_unit

    # Define coordinates for the center of rotation.
    x_rot_centers = x_initial

    # Update coordinates for the center of rotation from the subset defined by indices
    if rotation_center_indices is not None:
        if len(rotation_center_indices) < 1 or len(rotation_center_indices) > x_initial.shape[0]:
            raise IndexError("Length of rotation center indices must be >= % d and <= %d" % (1, x_initial.shape[0]))
        x_rot_centers = positions[rotation_center_indices, :]

    x_initial_mean = x_rot_centers.mean(0)

    # Generate a random rotation matrix.
    rotation_matrix = generate_random_rotation_matrix()

    # Apply rotation.
    x_proposed = (rotation_matrix * np.matrix(x_initial - x_initial_mean).T).T + x_initial_mean
    return x_proposed


def generate_conformers(input_wats, n_conf):

    wat = md.load_pdb(input_wats, no_boxchk=True)
    wat_atom_ids = wat.topology.select("water and name O")
    coords = md.utils.in_units_of(wat.xyz, "nanometers", "angstroms")
    new_coords = []
    res_id = 1
    for wat in wat_atom_ids:
        old_pos = coords[0, range(wat, wat+3), :]
        conf_id = 1
        for i in range(n_conf):
            new_pos = rotate_positions(old_pos, rotation_center_indices=[0])
            O_coords = np.array([new_pos[0, 0], new_pos[0, 1], new_pos[0, 2]])
            new_h1 = np.array([new_pos[1, 0], new_pos[1, 1], new_pos[1, 2]])
            new_h2 = np.array([new_pos[2, 0], new_pos[2, 1], new_pos[2, 2]])
            new_coords.append((O_coords, conf_id, res_id))
            new_coords.append((new_h1, conf_id, res_id))
            new_coords.append((new_h2, conf_id, res_id))
            conf_id += 1
        res_id += 1
    return new_coords


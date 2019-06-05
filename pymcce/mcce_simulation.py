import os, sys
import struct
import numpy as np
import mdtraj as md
from io import StringIO

HEAD3_HEADER = {'FL': 0, 'occ': 1, 'crg': 2, 'Em0': 3, 'pKa0': 4, 'ne': 5,
                'nH': 6, 'vdw0': 7, 'vdw1': 8, 'tors': 9, 'epol': 10,
                'dsolv': 11, 'extra': 12, 'history': 13}


class Simulation(object):
    """A class representing MCCE simulation data."""

    def __init__(self, ms_data_file, head3lst_file, fort38_file):

        self.ms_data_file = ms_data_file
        self.byte_indices = None
        self.total_microstates = 0

        res_list = []



        with open(self.ms_data_file, "rb") as md:
            bytes_n_res = md.read(4)
            n_res = struct.unpack('i', bytes_n_res)[0]
            for i in range(n_res):
                resname = str(md.read(8))
                res_list.append(resname)
            self.n_res = n_res
        self.residue_list = res_list
        # self.residue_hb_matrix = np.zeros((n_res, n_res), dtype=float)

        # with open(ms_gold_file, "r") as ms_gold:
        #    self.n_res = len([res.strip() for res in ms_gold.readlines() if len(res) != 0])

        self.conf_id_name_map = {}
        self.iConf = {}
        self.FL = {}
        self.crg = {}
        self.Em0 = {}
        self.pKa0 = {}
        self.ne = {}
        self.nH = {}
        self.vdw0 = {}
        self.vdw1 = {}
        self.tors = {}
        self.epol = {}
        self.dsolv = {}
        self.extra = {}

        data_holders = [self.iConf, self.FL, self.crg, self.Em0,
                        self.pKa0, self.ne, self.nH, self.vdw0, self.vdw1, self.tors, self.epol,
                        self.dsolv, self.extra]


        with open(head3lst_file, "r") as h3:
            for line in h3.readlines()[1:]:
                data = line.split()
                conf_id, conf_name,  = int(data[0]), data[1]
                self.conf_id_name_map[conf_id] = conf_name
                self.iConf[conf_name] = conf_id
                self.FL[conf_name] = data[2]
                self.crg[conf_name] = float(data[4])
                self.Em0[conf_name] = float(data[5])
                self.pKa0[conf_name] = float(data[6])
                self.ne[conf_name] = float(data[7])
                self.nH[conf_name] = float(data[8])
                self.vdw0[conf_name] = float(data[9])
                self.vdw1[conf_name] = float(data[10])
                self.tors[conf_name] = float(data[11])
                self.epol[conf_name] = float(data[12])
                self.dsolv[conf_name] = float(data[13])
                self.extra[conf_name] = float(data[14])


        # read fort38 data
        conf_occ = {}
        with open(fort38_file, 'r') as f:
            lines = [l.strip().split() for l in f.readlines()]
            header = lines.pop(0)
            for index, l in enumerate(lines):
                conf_id = index + 1
                conf_name = l[0]
                occ = float(l[1])
                if conf_name in self.iConf.keys():
                    #self.iConf[conf_name][1] += occ
                    conf_occ[conf_name] = occ
                else:
                    conf_occ[conf_name] = None

        self.conformer_occupancies = conf_occ

        residue_data = {}

        for k in sorted(self.iConf.keys()):
            conf_name = k
            res_key = conf_name[:3] + conf_name[5:10]
            if res_key not in residue_data.keys():
                residue_data[res_key] = [(k, self.iConf[k])]
            else:
                residue_data[res_key].append((k, self.iConf[k]))
        self.residue_data = residue_data


    def generate_byte_indices(self, sample_frequency, filename):
        """
        Generate byte indices

        Parameters
        ----------
        n_res : TYPE
            Description
        sample_frequency : int, optional
            Description
        filename : None, optional
            Description

        Returns
        -------
        rec_indices : list
            A list of the starting bytes for each record in the microstate data file
        """
        start_byte = 4 + (8 * self.n_res)
        bytes_per_record = (self.n_res * 2) + 20
        file_size = os.path.getsize(filename)
        # n_records = (file_size - start_byte) / bytes_per_record
        rec_indices = list(range(start_byte, file_size, sample_frequency * bytes_per_record))
        return rec_indices

    def parse_trajectory(self, sample_frequency=100):
        """
        Parse ms.dat
        """
        byte_indices = self.generate_byte_indices(sample_frequency, self.ms_data_file)
        total_records = len(byte_indices)
        trajectory = np.zeros([total_records, self.n_res], dtype=int)
        state_counts = np.zeros([total_records], dtype=int)
        energies = np.zeros([total_records], dtype=float)
        # print trajectory
        progress_counter = 0
        # print_progress_bar(progress_counter, self.total_records)
        with open(self.ms_data_file, "rb") as ms:
            for index, record in enumerate(byte_indices):
                ms.seek(record)
                bytes_conf_ids = ms.read(2 * self.n_res)
                bytes_energies_1 = ms.read(8)
                ms.seek(ms.tell() + 8)
                energy = struct.unpack("d", bytes_energies_1)[0]
                bytes_state_count = ms.read(4)
                trajectory[index, :] = np.asarray(struct.unpack(str(self.n_res) + "H", bytes_conf_ids))
                # print(struct.unpack(str(self.n_res) + "H", bytes_conf_ids)[-2:])
                state_count = struct.unpack("i", bytes_state_count)[0]
                self.total_microstates += state_count
                state_counts[index] += state_count
                energies[index] += energy
                progress_counter += 1
                # print_progress_bar(progress_counter, self.total_records)

        self.trajectory = trajectory
        self.state_counts = state_counts
        self.energies = energies

    def get_step2out_data(self, step2out_pdb):

        with open(step2out_pdb, "r") as f:
            atom_lines_step2out = [l.split() for l in f.readlines() if "ATOM" in l]
            conf_names = [l[3] + l[-1][:2] + l[4] for l in atom_lines_step2out]
            atom_names = [l[2] for l in atom_lines_step2out]
            charges = [float(l[-2]) for l in atom_lines_step2out]
            radii = [float(l[-3]) for l in atom_lines_step2out]
        return conf_names, charges, radii, atom_names

    def parse_struct(self, step2out_pdb):
        """
        Obtain structural data for each conformer
        """
        print("Parsing structure ...")
        st = md.load_pdb(step2out_pdb, frame=0, no_boxchk=True)
        conf_names, charges, radii, atom_names = self.get_step2out_data(step2out_pdb)

        conformer_atom_data = {}
        for index, conf in enumerate(conf_names):
            if conf not in conformer_atom_data.keys():
                conformer_atom_data[conf] = [[index, charges[index], radii[index], atom_names[index]]]
            else:
                conformer_atom_data[conf].append([index, charges[index], radii[index], atom_names[index]])

        self.structure = st
        self.conformer_atom_data = conformer_atom_data
        print("Done.")

    def calculate_dipoles(self):

        """

        Returns
        -------

        """
        local_indices = np.array([(a.index, a.residue.atom(0).index) for a in self.structure.top.atoms], dtype='int32')
        local_displacements = md.compute_displacements(self.structure, local_indices, periodic=True)

        molecule_indices = np.array([(a.residue.atom(0).index, 0) for a in self.structure.top.atoms], dtype='int32')
        molecule_displacements = md.compute_displacements(self.structure, molecule_indices, periodic=True)

        xyz = local_displacements + molecule_displacements
        conf_dipoles = {}
        for conf in self.conformer_atom_data:
            conf_at_indices = np.asarray([l[0] for l in self.conformer_atom_data[conf]])
            conf_charges = np.asarray([l[1] for l in self.conformer_atom_data[conf]])
            conf_xyz = xyz[0, conf_at_indices, :]
            moments = conf_xyz.transpose(1, 0).dot(conf_charges)
            moments *= 10.0
            conf_dipoles[conf] = moments

        return conf_dipoles

    def parse_opp_file(self, opp):

        pass
        """
        with open(opp, "r") as f1:
            all_lines = f1.readlines()

            energy_table = [[float(l[20:29]),
                             float(l[29:37]),
                             float(l[37:45]),
                             float(l[45:53])] for l in all_lines]

            return np.asarray(energy_table)
        """
    def parse_energy_data(self, opp_dir):

        energies_dir = opp_dir
        opp_files = [os.path.join(energies_dir, f) for f in os.listdir(energies_dir) if f.endswith(".opp")]
        opp_files = sorted(opp_files)
        n = len(self.iConf.keys())
        energy_matrix_vdw = np.zeros((n, n), dtype=float)
        energy_matrix_elec = np.zeros((n, n), dtype=float)
        conf_id_list = sorted(self.conf_id_name_map.keys())
        conf_names = [self.conf_id_name_map[c] for c in conf_id_list]
        for index, f in enumerate(conf_names):
            if "DM" not in f:
                opp = os.path.join(energies_dir, f + ".opp")
                print(opp)
                a = np.genfromtxt(opp, delimiter=(5, 15, 9, 8, 8, 8))[:, 2:4]
                energy_matrix_elec[index, :] = a[:, 0]
                energy_matrix_vdw[index, :] = a[:, 1]

        return energy_matrix_elec, energy_matrix_vdw
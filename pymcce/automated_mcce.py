#!/usr/local/bin/env python
# -*- coding: utf-8 -*-
"""
Automate MCCE Calculations.
Description
-----------
This module implements a program to set up and run automated mcce calculations.
Notes
-----
This is still in development.

Examples
--------
Coming soon!
"""


import sys
import time
from subprocess import call
import os
import shutil
import re

from collections import OrderedDict
from argparse import ArgumentParser


class MCCEParams(object):
    """
    A class representing an MCCE parameter file

    Parameters
    ----------
    mcce_directory : str
      A string containing full path of the MCCE installation directory.
    calculation_type : str, default=quick
      A string specifying the type of MCCE calculation, if not specified "quick" is used by default.

    Notes
    -----
      This class uses run.prm.quick as a template file to build an MCCE parameter object.
      Should there be an option to initialize with quick or full (?)
    """

    def __init__(self, mcce_directory, calculation_type=None):
        """Initialize parameters for an MCCE calculation.
        """

        self.mcce_directory = mcce_directory
        self.calculation_type = calculation_type
        if calculation_type not in ["quick", "full", "default", None]:
            sys.exit(
            "Unrecognized MCCE calculation type, allowed values are: quick, full or default.")
        self.mcce_params = self.load_params()

    def load_params(self):
        """Loads parameters for MCCE calculation from run.prm file in MCCE installation directory.

        Returns
        -------
        params : dict
          A dictionary of MCCE calculation parameters, each element of the dictionary is a key:
          key = string, containing the parameter name in the prm file, without parentheses
          value = list, first element is the value and second element is the description, both are read from the prm file
        """
        prm_source_file = ""
        if self.calculation_type is None:
            prm_source_file = open(self.mcce_directory + "/" + "run.prm", "r")
        else:
            prm_source_file = open(self.mcce_directory + "/" + "run.prm." + self.calculation_type, "r")
        
        prm_lines = prm_source_file.readlines()
        prm_source_file.close()
        params = OrderedDict()
        for line in prm_lines:
            words = line.split(" ")
            if re.search(r'\((.\w*)\)', words[-1]):
                parameter = words[-1].strip()
                parameter = parameter.strip("()\t\n")
                value = words[0]
                description = line.split(" ")[1:-1]
                if parameter == "MCCE_HOME":
                    value = self.mcce_directory
                if parameter in ["EXTRA", "RENAME_RULES"]:
                    value = self.mcce_directory + "/" + value.split("/")[-1]
                if parameter == "DELPHI_EXE":
                    value = self.mcce_directory + \
                        "/bin/" + value.split("/")[-1]
                if "//" in value:
                    value = value.replace("//", "/")
                params[parameter] = [value, " ".join(description)]
        return params

    def edit_parameters(self, **kwargs):
        """Edits MCCE parameters by updating values in MCCE parameter dictionary keys.

        Parameters
        ----------
        **kwargs : Arbitrary keyword arguments
          Any number of parameters can be specified as MCCE_PARAM="VALUE", where MCCE_PARAM is a valid
          parameter name and "VALUE" is a string with corresponding value.

        Examples
        --------
        >>> prm = MCCEParams("~/mcce/")
        >>> prm.edit_parameters(DO_PREMCCE="t", DO_ROTAMERS="t", DO_ENERGY="t", DO_MONTE="f")
        >>> prm.edit_parameters(DO_PREMCCE="t", DO_ROTAMERS="t")

        Notes
        -----
        During parameter editing, there is no way to check if legal values are passed on. 
        This is probably handled by MCCE initialization step. It is therefore assumed that this function is used responsibly.
        """
        for parameter in kwargs:
            update_value = kwargs[parameter]
            if parameter not in self.mcce_params.keys():
                raise KeyError(
                    "Could not find %s parameter in parameter dictionary, please supply valid parameters" % parameter)
            else:
                self.mcce_params[parameter][0] = update_value

    def write_runprm(self, destination_dir):
        """Writes run.prm file in a sub-directory

        Parameters
        ----------
        destination_dir : string
          String containing path of the directory, where run.prm will be saved.
        """
        runprm_filepath = destination_dir + "run.prm"
        runprm_file = open(runprm_filepath, "w")
        param_line_format = "{0:s} {1:s} ({2:s})\n"
        for parameter in self.mcce_params.keys():
            line_to_write = param_line_format.format(
                self.mcce_params[parameter][0], self.mcce_params[parameter][1], parameter)
            runprm_file.write(line_to_write)
        runprm_file.close()

    def write_submitsh(self, destination_dir, run_name=""):
        submit_text = ["#!/bin/sh\n", "#$ -S /bin/sh\n", "#$ -N mcce\n", "#$ -cwd\n",
                       "#$ -o run.log\n", "#$ -e error.log\n", self.mcce_directory + "/mcce"]
        submit_text[2] = submit_text[2].replace("mcce", "mcce_" + run_name)
        submitsh = open(os.path.join(destination_dir, "submit.sh"), "w")
        for line in submit_text:
            submitsh.write(line)
        submitsh.close()



def automated_run(input_dir, destination_dir, mcce_dir, local=False):
    """Performs an automated mcce calculation on a set of pdb file, located in an input folder.

    Parameters
    ----------
    input_dir : str
      String consisting of a valid path of a directory containing pdb files.
    destination_dir : str
      String consisting of a path representing a target directory where results of each mcce calculation will be stored.
    mcce_dir :  str
      Location of the mcce installation, should contain mcce executable.

    Notes
    -----
    This function will not work if each pdb file that should be processed is present in its own sub-directory. It is suggested that
    all pdb files are collected in the input directory before running this.
    """

    # check if input_dir exists and is not empty
    if not os.path.isdir(input_dir):
        sys.exit("Input directory not found.")
    # check if mcce exists
    input_pdb_files = [pdb_file for pdb_file in os.listdir(
        input_dir) if pdb_file.endswith(".pdb")]
    for pdb_file in input_pdb_files:
        print("Preparing input files for: ", pdb_file[0:-4])
        output_dir = destination_dir + "/mcee_results_" + pdb_file[0:-4]
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_dir = output_dir.replace("//", "/")
        source_pdb_file_path = input_dir + "/" + pdb_file
        target_pdb_file_path = output_dir + "/" + "prot.pdb"
        source_pdb_file_path = source_pdb_file_path.replace("//", "/")
        target_pdb_file_path = target_pdb_file_path.replace("//", "/")
        shutil.copy(source_pdb_file_path, target_pdb_file_path)
        # generate prm file
        os.chdir(output_dir)
        prm = MCCEParams(mcce_dir)
        prm.edit_parameters(DO_PREMCCE="t", DO_ROTAMERS="t",
                            DO_ENERGY="t", DO_MONTE="f")
        prm.write_runprm("")
        prm.write_submitsh("", run_name=pdb_file[0:-4])
        if not local:
            print("This command will be executed for: ", pdb)
            print("qsub submit.sh")
            call("qsub submit.sh", shell=True)
            time.sleep(10)
        else:
            pass


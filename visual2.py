import numpy as np
from math import radians, cos, sin
from ase import Atoms
from ase.io import Trajectory, write
from ase.build import bulk, surface, molecule, add_adsorbate, fcc111
from ase.constraints import ExpCellFilter, StrainFilter, FixAtoms, FixedPlane, FixBondLength
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.visualize import view
from ase.build.rotate import minimize_rotation_and_translation
from ase.io.vasp import read_vasp, write_vasp
import pandas as pd
import ipywidgets as widgets
from IPython.display import display_png, Image as ImageWidget
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
import shutil
import glob
from pathlib import Path
from PIL import Image, ImageDraw
from pymatgen.core import Lattice, Structure, Molecule
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor

import numpy as np
from scipy.spatial.distance import cdist

def visual2(structure, max_size=(100, 100), rotation='15z,-90x', celloff=True):
    def calculate_stretch_factor(cell, angle_x, angle_z):
        """Calculate the stretch factor for visualization based on the cell dimensions and rotation angles."""
        x_length = max(cell[0][0], 1e-5)  # Ensure non-zero length
        y_length = max(cell[1][1], 1e-5)  # Ensure non-zero length
        z_length = max(cell[2][2], 1e-5)  # Ensure non-zero length

        # Calculate effective lengths based on rotation
        cos_angle_z = cos(radians(angle_z))
        sin_angle_z = sin(radians(angle_z))
        cos_angle_x = cos(radians(angle_x))
        sin_angle_x = sin(radians(angle_x))

        effective_x_length = cos_angle_z * x_length + sin_angle_z * z_length
        effective_y_length = cos_angle_x * y_length + sin_angle_x * z_length
        # Calculate stretch factor, ensuring it's greater than zero
        stretch_y = max(np.abs(effective_y_length / effective_x_length), 1e-5)
        return stretch_y

    # Parse rotation angles
    angle_x = float(rotation.split(',')[1].replace('x', ''))
    angle_z = float(rotation.split(',')[0].replace('z', ''))

    # Calculate stretch factor based on the structure's cell and rotation angles
    stretch_y = calculate_stretch_factor(structure.cell, angle_x, angle_z)
    if celloff:
        structure.cell = None

    # Visualize the structure
    renderer = write('./temp.pov', structure, rotation=rotation)
    renderer.render()
    image_path = './temp.png'
    img = Image.open(image_path)

    # Calculate new size with stretch factor
    new_size = (max_size[0], int(max_size[1] * stretch_y))
    img = img.resize(new_size, Image.LANCZOS)
    display(img)

    # Move files to output directory
    destination = './output/'
    files = ['./temp.ini', './temp.pov', './temp.png']

    # Ensure destination directory exists
    os.makedirs(destination, exist_ok=True)

    for file in files:
        if os.path.isfile(os.path.join(destination, os.path.basename(file))):
            os.remove(os.path.join(destination, os.path.basename(file)))

        shutil.move(file, destination)

# Assuming visual is a function you have for visualization
for idx, structure in enumerate(traj):
    visual2(structure, max_size=(200, 200), rotation='0z,-0x',celloff=False)
    visual2(structure, max_size=(200, 200), rotation='0z,-90x',celloff=False)

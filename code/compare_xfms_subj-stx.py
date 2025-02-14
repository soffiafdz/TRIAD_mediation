#!/usr/bin/env python3

import os
from pathlib import Path
import numpy as np
import subprocess
import csv
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed


def invert_stx(xfm_file, out_file):
    """
    Use xfminvert to invert an xfm file.

    Args:
        xfm_file (str): Path to the xfm file.
        out_file (str): Path to the resulting inverted xfm file.

    Returns:
        None

    Raises:
        FileNotFoundError: If either `xfm_file` does not exist.
        subprocess.CalledProcessError: If the xfminvert command fails to run
            successfully.
    """
    # Convert arguments to Path objects
    xfm_file = Path(xfm_file)
    out_file = Path(out_file)

    # Check for file existance
    if not xfm_file.exists():
        raise FileNotFoundError(f"File not found: {xfm_file.name}")


    if not out_file.exists():
        # Create parent directories if unexistent
        if not out_file.parent.exists():
            out_file.parent.mkdir(parents=True, exist_ok=True)

        command = [
            "xfminvert", str(xfm_file), str(out_file)
        ]

        # Execute command
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            raise subprocess.CalledProcessError(
                e.returncode,
                e.cmd,
                output=e.output,
                stderr=e.stderr,
                msg="Error inverting .xfm file",
            )

def concat_stx_stx2(xfm_file1, xfm_file2, out_file):
    """
    Use xfmconcat to concatenate two xfm transformations.

    Args:
        xfm_file1 (str): Path to the first xfm file.
        xfm_file2 (str): Path to the second xfm file.
        out_file (str): Path to the resulting concatenated xfm file.

    Returns:
        None

    Raises:
        FileNotFoundError: If either `xfm_file1` or `xfm_file2` does not exist.
        subprocess.CalledProcessError: If the xfmtool command fails to run
            successfully.
    """
    # Convert arguments to Path objects
    xfm_file1 = Path(xfm_file1)
    xfm_file2 = Path(xfm_file2)
    out_file = Path(out_file)

    # Check for files existance
    if not xfm_file1.exists():
        raise FileNotFoundError(f"File not found: {xfm_file1.name}")

    if not xfm_file2.exists():
        raise FileNotFoundError(f"File not found: {xfm_file2.name}")

    if not out_file.exists():
        # Create parent directories if unexistent
        if not out_file.parent.exists():
            out_file.parent.mkdir(parents=True, exist_ok=True)

        command = [
            "xfmconcat", str(xfm_file1), str(xfm_file2), str(out_file)
        ]

        # Execute command
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            raise subprocess.CalledProcessError(
                e.returncode,
                e.cmd,
                output=e.output,
                stderr=e.stderr,
                msg="Error concatenating .xfm files",
            )

def extract_matrix_from_xfm(xfm_file):
    """
    Extract the 4x4 transformation matrix from the .xfm file.

    Args:
        xfm_file (str): Path to the .xfm file.

    Returns:
        np.ndarray: a 4x4 transformation matrix.
    """
    xfm_file = Path(xfm_file)
    with xfm_file.open('r') as file:
        lines = file.readlines()

    # Find the line containing "Linear_Transform"
    # and extract the following three lines
    for i, line in enumerate(lines):
        if "Linear_Transform" in line:
            # Read the next 3 lines containing the matrix
            matrix_lines = lines[i+1:i+4]
            matrix = []
            for matrix_line in matrix_lines:
                # Convert the values in the line into a list of floats
                matrix.append([
                    float(x.replace(';', ''))
                    for x in matrix_line.split()
                ])
            # Add homogeneous coordinates adjustment
            matrix.append([0, 0, 0, 1])
            # Return the 4x4 matrix as a NumPy array
            return np.array(matrix)


def apply_xfm(cube, xfm):
    """
    Apply the transformation to the cube corners.

    Args:
        cube (np.ndarray): The coordinates of the cube's corners
            (in homogeneous form).
        xfm (np.ndarray): The transformation matrix (4x4).

    Returns:
        np.ndarray: Transformed coordinates after applying the transformation.
    """
    # Apply the transformation
    transformed_cube = np.dot(cube, xfm.T)
    return(transformed_cube)


def compute_displacement(orig_cube, final_cube):
    """
    Compute the Euclidean displacement between the original and transformed
    cube corners.

    Args:
        orig_cube (np.ndarray): Original 3D coordinates of the cube corners.
        final_cube (np.ndarray): Final 3D coordinates after transformations.

    Returns:
        np.ndarray: Array of displacements for each corner.
    """
    # Ensure that the shapes of the matrices match before the comparison
    return np.linalg.norm(orig_cube[:, :3] - final_cube[:, :3], axis=1)


def process_subject_dir(subj_dir, cube_dims=(200, 200, 200)):
    """
    Process a subject directory with
    stx.xfm (native->subj.average) & stx2.xfm (native->subj.avg->stx.template)
    to 1) inverst stx.xfm, 2) concatenate stx.xfm & stx.xfm2, and 3) displace
    an arbitrary cube transformed by the resulting concatenated transformation.

    Args:
        subj_dir (str): The path to the subject directory containing
        the two xfm files.
        cube_dims (tuple): Tuple of length 3 with x, y, z dimensions for the
        cube.

    Returns:
        tuple: A tuple containing the subject's EID and
            an array (8x1) of displacements for each cube's corner.

    Raises:
        FileNotFoundError: If either stx.xfm or stx2.xfm does not exist.
        Exception: If there are issues reading or processing the .xfm files.
    """
    # Cube creation
    # Sanity checks
    if not isinstance(cube_dims, tuple):
        raise TypeError("cube_dims must be a tuple.")

    if len(cube_dims) != 3:
        raise ValueError("cube_dims must have exactly three elements (x,y,z).")

    if not all(isinstance(dim, (int, float)) for dim in cube_dims):
        raise TypeError("All elements in cube_dims must be integer or float.")

    x, y, z = (dim/2 for dim in cube_dims)

    cube_init = np.array([
        [-x, -y, -z, 1],
        [ x, -y, -z, 1],
        [-x,  y, -z, 1],
        [ x,  y, -z, 1],
        [-x, -y,  z, 1],
        [ x, -y,  z, 1],
        [-x,  y,  z, 1],
        [ x,  y,  z, 1]
    ])

    # Extract subject data from directory structure:
    # {Origin}/Subject-code/Visit-date/{files}
    subj_dir = Path(subj_dir)
    subj = subj_dir.name
    visits = [visit for visit in subj_dir.iterdir() if visit.is_dir()]
    for work_path in visits:
        visit = work_path.name

        # Find stx & stx2 files
        stx_file = work_path / f"stx_{subj}_{visit}_t1.xfm"
        stx2_file = work_path / f"stx2_{subj}_{visit}_t1.xfm"

        stx_files = [stx_file, stx2_file]

        for stx_file in [stx_file, stx2_file]:
            if not stx_file.exists():
                raise FileNotFoundError(f"File not found: {stx_file.name}")


        try:
            # Invert stx.xfm
            stx_inv = work_path / f"stx_{subj}_{visit}_t1_inv.xfm"

            if not stx_inv.exists():
                invert_stx(stx_file, stx_inv)

            # Concatenate stx_inv & stx2
            stx_corr = work_path / f"corr_{subj_dir}/{subj}_{visit}_t1.xfm"

            if not stx_corr.exists():
                concat_stx_stx2(stx_inv, stx2_file, stx_corr)

            # Extract the transformation matrix from the corrected .xfm file
            matrix_corr = extract_matrix_from_xfm(stx_corr)

            # Apply transformations
            cube_corr = apply_xfm(cube_init, matrix_corr)

            # Compute the displacement (Euclidean distance) between the corners
            displacement = compute_displacement(cube_init, cube_corr)
            return (subj, visit, displacement)

        except Exception as e:
            print(f"Error processing {subj}, {visit}: {e}")
            return None


# Define the source directory
code_dir = Path(__file__).resolve().parent
data_dir = code_dir.parent / 'data' / 'adni'

# ADNI1,2&Go & ADNI3
adni_dirs = [directory for directory in data_dir.iterdir()
             if "adni" in directory.name]

to_run_dirs = []
for adni_dir in adni_dirs:
    # Ensure that the directory with the xfms exists
    if not adni_dir.exists():
        raise FileNotFoundError(f"ADNI directory {adni_dir.name} not found.")

    # There's Tesla directories in ADNI1,2&GO

    if "1_2_go" in adni_dir.name:
        # Add the two Tesla subdirectories separately
        tesla_dirs = [subdir
                      for subdir in adni_dir.iterdir()
                      if subdir.is_dir() and "t" in subdir.name.lower()]
        to_run_dirs.extend(tesla_dirs)
    else:
        to_run_dirs.append(adni_dir)

# Iterate through the directories
results = {}
for work_dir in to_run_dirs:
    if "adni" in work_dir.name.lower():
        tesla = None
        adni_prot = work_dir.name.upper()
        current_list = results[adni_prot] = []
    else:
        tesla = work_dir.name.upper()
        adni_prot = work_dir.parent.name.upper()
        current_list = results[f"{adni_prot}—{tesla}"] = []

    # Get all subject directories
    subj_dirs = [entry.path
                 for entry in os.scandir(work_dir)
                 if entry.is_dir()]

    # Initialize tqdm progress bar with the total number of subject directories
    with tqdm(total=len(subj_dirs),
              desc=("Comparing subjects from"
                    f"{adni_prot}: {tesla}" if tesla else f"{adni_prot}"),
              unit="subjects") as pbar:
        # Process directories in parallel using ThreadPoolExecutor
        with ThreadPoolExecutor() as executor:
            futures = {executor.submit(process_subject_dir, subj_dir):
                       subj_dir for subj_dir in subj_dirs}

            # Use as_complete to update the progress bar as futures complete
            for future in as_completed(futures):
                result = future.result()
                if result is not None:
                    current_list.append(result)
                # Update the progress bar as each directory is processed
                pbar.update(1)

# Filter out any None results and write them to CSV
outcsv = data_dir.parent / 'xfm_displacements.csv'
results = [(key, value)
           for key, values in results.items()
           for value in values if value is not None]

with outcsv.open('w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['PROTOCOL'] +
                    ['SUBJ'] +
                    ['VISIT'] +
                    [f'Corner{i}' for i in range(1,9)])
    for adni, (subj, visit, displacements) in tqdm(results,
                                                   desc="Writing CSV",
                                                   unit="rows",
                                                   total=len(results),
                                                   ):
        writer.writerow([adni, subj, visit] + displacements.tolist())

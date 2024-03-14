# Project-ActivMol

## Description

ActivMol is a software tool designed to visualize the active sites of structurally characterized Metalloproteins. It is built using Python and the PyOpenGL library, with protein structure information parsed from PDB files using the Biopython library. The user interface is implemented with PyQt5, allowing for interactive visualization.

## Video Demo

A video demonstration of Project-ActivMol can be found [here](https://drive.google.com/file/d/1yNi-xEyisEsdM72UPDX3zPqhwck-W/view?usp=share_link).

To install project:
pip install -r requirements.txt

## Usage

- Ensure you have a stable internet connection for downloading PDB files.
- Enter a valid PDB ID to visualize. Incorrect IDs or IDs without metal atoms may cause the software to crash.
- Press the 'Submit' button and wait for the PDB file to download in the background.
- The visualization window will appear, displaying the active site of the protein within a GLWidget and control buttons on the right side.

### Camera Controls

Activate the keyboard controls by clicking once inside the GLWidget (black screen), then use the following keys:

- `W`: Move the camera forward (zoom in)
- `S`: Move the camera backward (zoom out)
- `A`: Move the camera left
- `D`: Move the camera right
- `Q`: Move the camera up
- `E`: Move the camera down

### Additional Features

- "Save Image" button: Saves the current view of the GLWidget to the "pdb/{ID}" folder as "ID_{index}".
- "Back" button: Returns to the main window for new PDB file visualization.
- "Next Active Site" and "Previous Active Site" buttons: Navigate through multiple active sites within the protein.
- "Export Site" button: Exports information about the currently active site's coordination sphere.

## Active Site Information

After submitting a PDB ID, a corresponding folder is created within the "pdb" directory containing the PDB file and active site information. For example, submitting "1si4" creates a "pdb1si4" folder with files `Pdb1si4.ent` and `1si4`, the latter detailing the active site interactions[1].

![Visualization Window](https://cdn.mathpix.com/cropped/2024_03_14_850cce79f38d605d0f39g-2.jpg?height=848&width=1628&top_left_y=796&top_left_x=304)

![Active Site Information](https://cdn.mathpix.com/cropped/2024_03_14_850cce79f38d605d0f39g-4.jpg?height=872&width=1071&top_left_y=320&top_left_x=423)

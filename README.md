# Project-ActivMol

## Description

ActivMol is a software tool designed to visualize the active sites of structurally characterized Metalloproteins. It is built using Python and the PyOpenGL library, with protein structure information parsed from PDB files using the Biopython library. The user interface is implemented with PyQt5, allowing for interactive visualization.

## Video Demo

A video demonstration of Project-ActivMol can be found [here](https://drive.google.com/file/d/1yNi-x__EyisEsdM72UPDX3zPq-hwcK-W/view?usp=drive_link).

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

After submitting a PDB ID, a corresponding folder is created within the "pdb" directory containing the PDB file and active site information. For example, submitting "1si4" creates a "pdb1si4" folder with files `Pdb1si4.ent` and `1si4`, the latter detailing the active site interactions.

## Results

- The analysis of active sites in structurally characterized metalloproteins revealed significant insights into the coordination sphere of these proteins. 
- For the protein "1si4" (Human hemoglobin A2), four active sites corresponding to each chain were identified, showcasing an octahedral geometry with a coordination number of 6. This geometry allows for the visualization of multiple active sites within the software, enhancing the understanding of protein structure and function.

## Conclusion and Future Works

- **Project Completion**: The development of ActivMol has successfully enabled the visualization of multiple active sites in metalloproteins, along with providing detailed information about the bond distances from donor atoms.
- **Software Limitations**: 
  - The software currently crashes with incorrect PDB IDs or when a protein without active sites/metal atoms is submitted (e.g., "6joy", "6zhs"). 
  - The visualization inaccurately represents some donor atoms (e.g., should be C instead of N for Cyanide), due to the exclusion of C as a donor atom for generality, leading to distance inaccuracies (e.g., FE and N atom distance being 3.09A, more than 3A).
- **Recommendations for Improvement**:
  - Implement error handling for incorrect PDB IDs and proteins without active sites/metal atoms.
  - Revise the donor atom representation to include C for cases like Cyanide to improve accuracy.
  - Add mouse controls for camera movement for user convenience.
  - Improve camera focus on metals in the active site, especially in cases of multiple metal atoms.
  - Introduce labeling of residues for better information accessibility.
  - Include a feature for measuring bond angles by selecting any three atoms on the screen.
  - Allow users to change the background of saved images for customization.
  - Enhance the UI for a smoother user experience.

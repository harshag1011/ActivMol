import itertools

from PyQt5 import QtWidgets, QtGui
from PyQt5.QtCore import Qt, QPoint
from PyQt5.QtWidgets import QApplication, QDesktopWidget, QFileDialog
from OpenGL.GLUT import glutBitmapCharacter
from PyQt5.QtWidgets import QOpenGLWidget
from PyQt5.QtGui import QImage, QColor, QIcon
from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt5.QtCore import Qt
from utility import *
import os
import shutil

active_sites = []
HETATM_name_dic = []
ATOM_name_dic = []
file_name_id = None


class GLWidget(QOpenGLWidget):

    def __init__(self, parent=None):
        super().__init__(parent)
        self.i = 0
        # self.active_site_file_list = []
        self.active_site_file_list = set()

    def next_site(self):
        # self.current_site_index = (self.current_site_index + 1) % len(self.active_sites)
        self.i = (self.i + 1) % len(active_sites)
        #print("value of i is ", self.i)

    def previous_site(self):
        # self.current_site_index = (self.current_site_index - 1) % len(self.active_sites)
        self.i = self.i - 1
        #print("value of i is ", self.i)

    def initializeGL(self):
        glClearColor(0.0, 0.0, 0.0, 1.0)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHT0)
        glEnable(GL_LIGHTING)
        glEnable(GL_COLOR_MATERIAL)
        self.setFocusPolicy(Qt.ClickFocus)

        # initialize the camera position and orientation
        self.camera_pos = [0.0, 0.0, 5.0]
        self.camera_target = [-1, -1, -1]
        self.camera_up = [0, 1, 0]


    def keyPressEvent(self, event):
        super().keyPressEvent(event)
        key = event.key()
        if key == Qt.Key_W:
            # Move camera forward along its local z-axis
            self.camera_pos[2] -= 0.2
            self.update()
        elif event.key() == Qt.Key_S:
            # Move camera backward along its local z-axis
            self.camera_pos[2] += 0.2
            self.update()
        elif event.key() == Qt.Key_A:
            # Move camera left along its local x-axis
            self.camera_pos[0] -= 0.2
            self.update()
        elif event.key() == Qt.Key_D:
            # Move camera right along its local x-axis
            self.camera_pos[0] += 0.2
            self.update()
        elif event.key() == Qt.Key_Q:
            # Move camera up along its local y-axis
            self.camera_pos[1] += 0.2
            self.update()
        elif event.key() == Qt.Key_E:
            # Move camera down along its local y-axis
            self.camera_pos[1] -= 0.2
            self.update()

    def focusInEvent(self, event):
        super().focusInEvent(event)
        # print('focus in: ')

    def resizeGL(self, width, height):
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        aspect_ratio = width / height
        gluPerspective(45.0, aspect_ratio, 0.1, 100.0)
        glMatrixMode(GL_MODELVIEW)

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        gluLookAt(self.camera_pos[0], self.camera_pos[1], self.camera_pos[2],
                  self.camera_target[0], self.camera_target[1], self.camera_target[2],
                  self.camera_up[0], self.camera_up[1], self.camera_up[2])

        # print("For a particular site :\n", active_sites[0])
        active_site_atoms = []
        active_site_metal_atoms = []
        active_site_residues = set()
        # getting site's metal atoms, atoms(for normalizebox), residues
        for ele in active_sites[self.i]:  ################
            active_site_metal_atoms.append(ele[0])  # appending all metal atoms
            active_site_atoms.append(ele[0])
            for res in ele[1]:
                active_site_residues.add(res)

        metal_centers = active_site_metal_atoms  # (without oxygetn)

        # appending residue atoms to active_site_atoms which already contains site's metal atoms
        for res in active_site_residues:
            l1, l2 = get_atoms_of_residue(res)
            active_site_atoms.append(l1)
            if l2 != None:
                active_site_metal_atoms.append(l2)  # for O in water
                active_site_metal_atoms = flat(active_site_metal_atoms)

        active_site_atoms = flat(active_site_atoms)
        normalize_center, normalize_size = normalize_box(active_site_atoms)
        metal_coords = normalize_metal_atoms(active_site_metal_atoms, normalize_center, normalize_size, HETATM_name_dic)

        b_coords_all = []  # b_coords of all residues of an active site
        for res in active_site_residues:
            l1, _ = get_atoms_of_residue(res)
            if len(l1) != 1:
                atoms = flat(l1)
            else:
                atoms = []

            if len(atoms) != 0:
                bonded_atoms = find_bonded_atoms(atoms, HETATM_name_dic, ATOM_name_dic)
                b_coord = get_bonded_coords(bonded_atoms, HETATM_name_dic,
                                            ATOM_name_dic)  # [(x1,y1,z1),(x2,y2,z2),(ele1,ele2)]
                b_coords_all.append(normalize_residue(b_coord, normalize_center, normalize_size))

        atoms = []
        for res in active_site_residues:
            ls = Selection.unfold_entities(res, 'A')
            atoms.append(ls)
        atoms = flat(atoms)
        metal_centers = flat(metal_centers)

        # remove water from metal_centers
        temp = []
        for i in metal_centers:
            res_id = i.get_full_id()[3]
            if res_id[0][0] != "W":
                temp.append(i)
        #print(temp)
        metal_centers = temp
        if len(self.active_site_file_list) != 0:
            self.active_site_file_list.clear()

        for atom_metal in metal_centers:
            donor_atoms = find_interaction_atoms(atom_metal, atoms, HETATM_name_dic,
                                                 ATOM_name_dic)  # for active_site_metal[0]

            self.active_site_file_list.add((atom_metal,tuple(donor_atoms)))




            # print(donor_atoms)
            my_pairs = list(itertools.combinations(donor_atoms, 2))

            for donor in donor_atoms:
                start_point = normalize_atom(atom_metal, normalize_center, normalize_size)
                end_point = normalize_atom(donor, normalize_center, normalize_size)

                # Set up line properties
                glColor3f(3.0, 3.0, 3.0)
                glLineWidth(5)

                # Define dashed line pattern
                glEnable(GL_LINE_STIPPLE)
                glLineStipple(5, 0xAAAA)

                # Draw the line
                glBegin(GL_LINES)
                glVertex3f(start_point[0], start_point[1], start_point[2])
                glVertex3f(end_point[0], end_point[1], end_point[2])
                glEnd()

                # Disable dashed line mode
                glDisable(GL_LINE_STIPPLE)

            for pair in my_pairs:
                start_point = normalize_atom(pair[0], normalize_center, normalize_size)
                end_point = normalize_atom(pair[1], normalize_center, normalize_size)

                # Set up line properties
                glColor3f(2.0, 2.0, 0.0)
                glLineWidth(1)

                # Define dashed line pattern
                glEnable(GL_LINE_STIPPLE)
                glLineStipple(5, 0xAAAA)

                # Draw the line
                glBegin(GL_LINES)
                glVertex3f(start_point[0], start_point[1], start_point[2])
                glVertex3f(end_point[0], end_point[1], end_point[2])
                glEnd()

                # Disable dashed line mode
                glDisable(GL_LINE_STIPPLE)

        for m_coords in metal_coords:
            x, y, z = m_coords[0], m_coords[1], m_coords[2]
            name = m_coords[3]
            if name == 'O':
                glColor3f(1.0, 0.0, 0.0)
            else:
                glColor3f(1.0, 1.0, 1.0)  # white metal atoms
            glTranslatef(x, y, z)
            quad = gluNewQuadric()
            # gluSphere(quad, 0.3, 50, 50)
            gluSphere(quad, 0.10, 50, 50)
            gluDeleteQuadric(quad)
            glLoadIdentity()
            gluLookAt(self.camera_pos[0], self.camera_pos[1], self.camera_pos[2],
                      self.camera_target[0], self.camera_target[1], self.camera_target[2],
                      self.camera_up[0], self.camera_up[1], self.camera_up[2])

        for v in b_coords_all:
            for b_coord in v:
                x1, y1, z1 = b_coord[0][0], b_coord[0][1], b_coord[0][2]
                x2, y2, z2 = b_coord[1][0], b_coord[1][1], b_coord[1][2]
                ele_name1 = b_coord[2][0]
                ele_name2 = b_coord[2][1]

                if x2 < x1:
                    x1, x2 = x2, x1
                    y1, y2 = y2, y1
                    z1, z2 = z2, z1
                    ele_name1, ele_name2 = ele_name2, ele_name1

                # draw sphere
                if ele_name1 == 'C':
                    glColor3f(0.0, 1.0, 0.0)
                elif ele_name1 == 'O':
                    glColor3f(1.0, 0.0, 0.0)
                elif ele_name1 == 'N':
                    glColor3f(0.0, 0.0, 1.0)
                elif ele_name1 == 'P':
                    glColor3f(1.0, 0.49, 0.0)
                elif ele_name1 == 'S':
                    glColor3f(0.89, 0.77, 0.25)
                else:
                    glColor3f(1.0, 1.0, 1.0)

                glTranslatef(x1, y1, z1)  # position
                quad = gluNewQuadric()
                # gluSphere(quad, 0.3, 50, 50)
                gluSphere(quad, 0.10, 50, 50)
                gluDeleteQuadric(quad)
                glLoadIdentity()
                gluLookAt(self.camera_pos[0], self.camera_pos[1], self.camera_pos[2],
                          self.camera_target[0], self.camera_target[1], self.camera_target[2],
                          self.camera_up[0], self.camera_up[1], self.camera_up[2])  # draw sphere
                if ele_name2 == 'C':
                    glColor3f(0.0, 1.0, 0.0)
                elif ele_name2 == 'O':
                    glColor3f(1.0, 0.0, 0.0)
                elif ele_name2 == 'N':
                    glColor3f(0.0, 0.0, 1.0)
                elif ele_name2 == 'P':
                    glColor3f(1.0, 0.49, 0.0)
                elif ele_name2 == 'S':
                    glColor3f(0.89, 0.77, 0.25)
                else:
                    glColor3f(1.0, 1.0, 1.0)

                glTranslatef(x2, y2, z2)  # position
                quad = gluNewQuadric()
                # gluSphere(quad, 0.3, 50, 50)
                gluSphere(quad, 0.10, 50, 50)
                gluDeleteQuadric(quad)
                glLoadIdentity()
                gluLookAt(self.camera_pos[0], self.camera_pos[1], self.camera_pos[2],
                          self.camera_target[0], self.camera_target[1], self.camera_target[2],
                          self.camera_up[0], self.camera_up[1], self.camera_up[2])  # Set color to white
                glColor3f(1.0, 1.0, 1.0)
                # Draw cylinder
                d = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
                z_rot = -math.pi / 2 + math.atan((y2 - y1) / (x2 - x1))
                x_rot = math.asin((z2 - z1) / d)
                # cylinder_radius = 0.15
                cylinder_radius = 0.05
                cylinder_height = d
                glTranslatef(x1, y1, z1)
                glRotatef((z_rot * 180) / (math.pi), 0.0, 0.0, 1.0)
                glRotatef((x_rot * 180) / (math.pi), 1.0, 0.0, 0.0)

                num_steps = 20  # Number of steps to create a circle
                angle_step = 2 * math.pi / num_steps

                glBegin(GL_TRIANGLE_STRIP)
                for i in range(num_steps + 1):
                    angle = i * angle_step
                    x = cylinder_radius * math.cos(angle)
                    z = cylinder_radius * math.sin(angle)

                    glVertex3f(x, 0.0, z)
                    glVertex3f(x, cylinder_height, z)
                glEnd()
                glLoadIdentity()
                gluLookAt(self.camera_pos[0], self.camera_pos[1], self.camera_pos[2],
                          self.camera_target[0], self.camera_target[1], self.camera_target[2],
                          self.camera_up[0], self.camera_up[1], self.camera_up[2])

    def export_data_active_site(self):
        # Generate information to be saved in a file
        # active_site = active_sites[self.i]
        #print(self.active_site_file_list)
        #data = "Active site information"

        # Open a file dialog to get the file path
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getSaveFileName(self, "Save File", "", "Text Files (*.txt)", options=options)

        # Write data to file
        if file_path:
            with open(file_path, 'w') as f:
                # f.write(data)
                #print(self.active_site_file_list,file=f)
                j = 1
                for x in self.active_site_file_list:
                    #print(type(x[0]),file=f)
                    metal = x[0]
                    donors = x[1] # tuple
                    print("Coordination Sphere - {}".format(j),file=f)
                    print("Metal Center -",metal.get_name(),file=f)
                    print("Full id(structure id, model id, chain id, residue id, atom id)",file=f)
                    print(metal.get_full_id(),file=f)
                    print(file=f)
                    print("Donors({}) are -".format(len(donors)),file=f)
                    for donor in donors:
                        print(donor.get_name(),file=f,end='  ')
                        print(donor.get_full_id(),file=f)
                        print("res name_id {}_{}".format(donor.get_parent().get_resname(),donor.get_parent().get_id()[1]),file=f) # res name and id of a donor
                        print("distance from metal(angstrom) -",round(donor-metal,2),file=f)
                    j +=1
                    print(file=f)
                f.close()

    def capture_image(self):
        # capture the current OpenGL screen image
        img = glReadPixels(0, 0, self.width(), self.height(), GL_RGBA, GL_UNSIGNED_BYTE)
        img = np.frombuffer(img, np.uint8).reshape(self.height(), self.width(), 4)
        img = QImage(img.data, self.width(), self.height(), QImage.Format_RGBA8888)
        # img = img.transformed(QTransform().rotate(90))
        img = img.mirrored(vertical=True)  # Mirror the image about the x-axis
        id = file_name_id
        # folder_path = os.path.join(os.getcwd(), "images/{}".format(id))
        # os.makedirs(folder_path, exist_ok=True)
        # save the image to the current directory
        img.save("./pdb/pdb{}/{}_{}.png".format(id, id, self.i + 1))  # 1si4_1, 1 is active site 1

    def on_slider_value_changed(self, new_value):
        print(new_value)
        modify_variable(new_value)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("ActivMol")
        self.setWindowIcon(QIcon('icon.png'))
        # get the size of the desktop screen
        desktop = QApplication.desktop()
        screenRect = desktop.screenGeometry(desktop.primaryScreen())

        # calculate the center point of the screen
        center = screenRect.center()
        # print(center)
        # move the MainWindow to the center of the screen
        self.move(int(QPoint.x(center) - 1.2 * self.width()), int(QPoint.y(center) - self.height()))

        # Create the two screens as separate QWidget instances
        screen1 = QtWidgets.QWidget()
        screen2 = QtWidgets.QWidget()

        # Create the input widgets
        self.pdb_edit = QtWidgets.QLineEdit()
        self.pdb_edit.setPlaceholderText("Enter PDB id here...")
        self.pdb_edit.setMinimumHeight(30)
        font = QtGui.QFont('Consolas', 12)
        self.pdb_edit.setFont(font)
        self.pdb_edit.setAlignment(Qt.AlignCenter)
        self.pdb_edit.setMaxLength(4)
        # self.pdb_edit.setFixedWidth(100)
        self.submit_button = QtWidgets.QPushButton("Submit")
        self.submit_button.clicked.connect(lambda: self.stacked_widget.setCurrentIndex(1))
        # Connect the submit button to the load_pdb function
        self.submit_button.clicked.connect(self.load_pdb)
        self.label = QtWidgets.QLabel()
        self.label.setText("ActivMOL")
        font = QtGui.QFont()
        font.setFamily("Microsoft JhengHei UI Light")
        font.setPointSize(18)
        self.label.setFont(font)
        self.label.setAlignment(Qt.AlignCenter)

        # Create the GLWidget and other necessary widgets
        self.glWidget = GLWidget()
        # push_button2 = QtWidgets.QPushButton("Hello")
        # push_button3 = QtWidgets.QPushButton("World")
        self.save_button = QtWidgets.QPushButton("Save Image")
        self.back_button = QtWidgets.QPushButton("Back")
        self.export_data_active_site = QtWidgets.QPushButton('Export Site', self)
        self.back_button.clicked.connect(lambda: self.stacked_widget.setCurrentIndex(0))
        self.export_data_active_site.clicked.connect(self.glWidget.export_data_active_site)

        self.save_button.clicked.connect(self.glWidget.capture_image)
        self.next_button = QtWidgets.QPushButton("Next Active Site")
        self.previous_button = QtWidgets.QPushButton("Previous Active Site")
        self.next_button.clicked.connect(self.glWidget.next_site)
        self.previous_button.clicked.connect(self.glWidget.previous_site)
        # Add other widgets as necessary

        # Add the input widgets to a layout
        screen1_layout = QtWidgets.QVBoxLayout()
        screen1_layout.addWidget(self.label)
        screen1_layout.addSpacing(50)
        screen1_layout.addWidget(self.pdb_edit)
        screen1_layout.addSpacing(10)
        screen1_layout.addWidget(self.submit_button)
        screen1_layout.setContentsMargins(700, 350, 700, 450)  # left, top, right, bottom
        screen1_layout.setSpacing(2)
        screen1.setLayout(screen1_layout)

        layout2_2 = QtWidgets.QVBoxLayout()
        layout2_3 = QtWidgets.QHBoxLayout()
        screen2_layout = QtWidgets.QHBoxLayout()
        screen2_layout.addWidget(self.glWidget)
        self.glWidget.setFixedWidth(1200)
        layout2_3.addWidget(self.previous_button)
        layout2_3.addWidget(self.next_button)
        layout2_2.addLayout(layout2_3)
        layout2_2.addWidget(self.save_button)
        layout2_2.addWidget(self.export_data_active_site)
        layout2_2.addWidget(self.back_button)
        # screen2_layout.addWidget(push_button2)
        # screen2_layout.addWidget(push_button3)
        screen2_layout.addLayout(layout2_2)
        screen2.setLayout(screen2_layout)

        # Create the stacked widget and add the screens as widgets to it
        self.stacked_widget = QtWidgets.QStackedWidget()
        self.stacked_widget.addWidget(screen1)
        self.stacked_widget.addWidget(screen2)

        # Set the central widget of the MainWindow to be the stacked widget
        self.setCentralWidget(self.stacked_widget)



    def load_pdb(self):
        # Retrieve the PDB ID from the QLineEdit widget
        file_id = self.pdb_edit.text()
        # Load the PDB file and extract the necessary data
        structure = parse_pdb(file_id)
        os.makedirs("./pdb/pdb{}".format(file_id), exist_ok=True)
        shutil.copy('pdb{}.ent'.format(file_id), './pdb/pdb{}'.format(file_id))
        os.remove("./pdb{}.ent".format(file_id))
        ATOM_name, HETATM_name = atoms_name(file_id)
        metal_atoms, atom_all = find_metal_atoms(structure, file_id)
        pos_sites = find_pos_sites(metal_atoms, atom_all)
        sites_active = find_active_sites(pos_sites)
        print_active_sites(sites_active,file_id)
        # print(sites_active)
        global active_sites
        active_sites = sites_active
        global ATOM_name_dic
        ATOM_name_dic = ATOM_name
        global HETATM_name_dic
        HETATM_name_dic = HETATM_name
        global file_name_id
        file_name_id = file_id


if __name__ == '__main__':
    app = QApplication([])
    window = MainWindow()
    window.show()
    app.exec_()

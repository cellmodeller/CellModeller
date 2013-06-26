import sys, math
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *



class Viewer:

    def __init__(self):
        self.view_x = 0
        self.view_y = 0
        self.view_z = 400

        self.width = 1200
        self.height = 1200
        self.window_name = 'CLBacteriumViewer'

        self.bg_color = (0.7, 0.7, 0.7)
        self.grid_color = (0, 0, 0)
        self.cell_outline_color = (0, 0, 0)
        self.cell_color = (0.5, 0.5, 0.5)

        self.circ_pts = [(math.cos(math.radians(th)), math.sin(math.radians(th))) for th in range(-80,90,20)]

        self.n_cells = 0
        self.cells = []

        self.init_glut()
        self.init_gl()


    def init_glut(self):
        glutInit(sys.argv)
        glutInitWindowSize(self.width, self.height)
        glutInitWindowPosition(0, 0)
        glutCreateWindow(self.window_name)
        glutDisplayFunc(lambda: self.display())
        glutReshapeFunc(lambda w,h: self.reshape(w,h))
        glutKeyboardFunc(lambda *args: self.handle_keyboard(*args))
        glutIdleFunc(lambda: self.idle())


    def init_gl(self):
        glEnable(GL_LINE_SMOOTH)
        glEnable(GL_POLYGON_SMOOTH)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)


    # def cell_color(self, i):
    #     while i not in self.founders:
    #         i = self.parent[i]
    #     return self.founders[i]


    def place_camera(self):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(60.0, 1.0, 0.1, 1000.0)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glTranslatef(self.view_x, self.view_y, -self.view_z)


    def display(self):
        glClearColor(*(self.bg_color + (0,)))
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        self.place_camera()
        for i in range(self.n_cells):
            center,dir,len,rad = self.cells[i]
            self.display_cell(center, dir, len, rad)
        glFlush()
        glutSwapBuffers()



    def display_cell(self, p, d, l, r):
        glColor3f(*self.cell_color)
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
        self._display_cell(p, d, l, r)
        glColor3f(0.0, 0.0, 0.0)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
        glLineWidth(0.5)
        self._display_cell(p, d, l, r)


    def _display_cell(self, p, d, l, r):
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        ang = math.atan2(d[1], d[0]) * 360.0 / (2.0*3.141593)
        glTranslatef(p[0], p[1], 0.0)
        glRotatef(ang, 0.0, 0.0, 1.0)
        glBegin(GL_POLYGON)
        glVertex3f(-l/2.0, -r, 0)
        glVertex3f(l/2.0, -r, 0)
        for x,y in self.circ_pts:
            glVertex3f(l/2.0 + x*r, y*r, 0.0)
        glVertex3f(l/2.0, r, 0)
        glVertex3f(-l/2.0, r, 0)
        for x,y in self.circ_pts:
            glVertex3f(-l/2.0 -x*r, -y*r, 0.0)
        glEnd()
        glPopMatrix()


    def reshape(self, w, h):
        l = min(w, h)
        glViewport(0, 0, l, l)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glutPostRedisplay()


    def handle_keyboard(self, *args):
        if args[0] == 'j':
            self.view_x += 2
        elif args[0] == 'l':
            self.view_x -= 2
        elif args[0] == 'i':
            self.view_y -= 2
        elif args[0] == 'k':
            self.view_y += 2
        elif args[0] == 'e':
            self.view_z -= 2
        elif args[0] == 'd':
            self.view_z += 2
        elif args[0] == '\x1b':
            exit()


    def idle(self):
        #self.display()
        pass





if __name__ == '__main__':
    import sys
    import cPickle
    filename = sys.argv[1]
    n,ps,ds,rs,ls = cPickle.load(open(filename, 'rb'))
    viewer = Viewer()
    viewer.n_cells = n
    viewer.cells = zip(ps, ds, rs, ls)

    print 'n_cells:',n

    glutMainLoop()

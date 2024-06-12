from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
import pyqtgraph.opengl as gl
import pyqtgraph as pg
import pyqtgraph.parametertree.parameterTypes as pTypes
from pyqtgraph.parametertree import Parameter, ParameterTree
import OpenGL.GL as GL
import numpy as np
import sys

import Physics as ph
import Initialization as init

class ComplexParameter(pTypes.GroupParameter):
	def __init__(self, **opts):
		opts['type'] = 'bool'
		opts['value'] = True
		pTypes.GroupParameter.__init__(self, **opts)

		self.addChild({'name': 'Young Modules (E)', 'type': 'float', 'value': 2.93, 'siPrefix': True, 'step': 0.01})
		self.addChild({'name': 'Poisson ratio (Nu)', 'type': 'float', 'value': 0.425, 'siPrefix': True, 'step': 0.005})
		self.addChild({'name': 'Volume constant (kv)', 'type': 'float', 'value': 1.4, 'siPrefix': True, 'step': 0.1})
		self.addChild({'name': 'Density constant (ro_uni)', 'type': 'float', 'value': 2, 'siPrefix': True, 'step': 0.1})
		self.addChild({'name': 'Gamma 1', 'type': 'float', 'value': 0, 'siPrefix': True, 'step': 0.01})
		self.addChild({'name': 'Gamma 2', 'type': 'float', 'value': 0, 'siPrefix': True, 'step': 0.01})

		self.E = self.param('Young Modules (E)')
		self.Nu = self.param('Poisson ratio (Nu)')
		self.kv = self.param('Volume constant (kv)')
		self.ro_uni = self.param('Density constant (ro_uni)')
		self.Gamma1 = self.param('Gamma 1')
		self.Gamma2 = self.param('Gamma 2')

		self.E.sigValueChanged.connect(self.E_Changed)
		self.Nu.sigValueChanged.connect(self.Nu_Changed)
		self.kv.sigValueChanged.connect(self.kv_Changed)
		self.ro_uni.sigValueChanged.connect(self.ro_uni_Changed)
		self.Gamma1.sigValueChanged.connect(self.Gamma1_Changed)
		self.Gamma2.sigValueChanged.connect(self.Gamma2_Changed)

	def E_Changed(self):
		ph.E = self.E.value()

	def Nu_Changed(self):
		ph.Nu = self.Nu.value()

	def kv_Changed(self):
		ph.kv = self.kv.value()

	def ro_uni_Changed(self):
		ph.ro_uni = self.ro_uni.value()

	def Gamma1_Changed(self):
		ph.Gamma1 = self.Gamma1.value()

	def Gamma2_Changed(self):
		ph.Gamma2 = self.Gamma2.value()


class GLWidget(QtWidgets.QWidget):
	def __init__(self, parent=None):
		super(GLWidget, self).__init__(parent)
		self.w = pg.GraphicsLayoutWidget()
		self.w.setBackground("w")
		# self.w = gl.GLViewWidget()
		self.plot = gl.GLViewWidget()
		self.w.setGeometry(0, 110, 1200, 960)
		self.w.show()
		self.w.setWindowTitle('Physics-based simulation')
		self.plot.setCameraPosition(distance=10, elevation=0, azimuth=0)

		## Create some widgets to be placed inside
		Run = QtWidgets.QPushButton('Run')
		params =[ComplexParameter(name = 'Model Parameters')]

		## Create tree of Parameter objects
		p = Parameter.create(name='params', type='group', children=params)
		t = ParameterTree()
		t.setParameters(p, showTop=True)

		## CheckBox for external force
		self.checkBox1 = QtWidgets.QCheckBox('External Force on Sphere 1')
		self.checkBox2 = QtWidgets.QCheckBox('External Force on Sphere 2')
		self.checkBox3 = QtWidgets.QCheckBox('Force Vectors')


		## Create a grid layout to manage the widgets size and position
		layout = QtWidgets.QGridLayout()
		self.w.setLayout(layout)
		layout.addWidget(Run, 0, 0, 1, 1)  # button goes in upper-left
		layout.addWidget(t, 1, 0, 1, 1)
		layout.addWidget(self.plot, 1, 1, 4, 170)  # plot goes on right side, spanning 2 rows
		layout.addWidget(self.checkBox1, 2, 0, 1, 1)
		layout.addWidget(self.checkBox2, 3, 0, 1, 1)
		layout.addWidget(self.checkBox3, 4, 0, 1, 1)
		# layout.addWidget(t2, 1, 1, 1, 1)

		# Widgets changed
		self.run = False
		self.run_once = False
		Run.clicked.connect(self.run2true)

		self.initializeGL()
		self.paintGL()
		# self.w.setBackgroundColor(255, 255, 255, 255)


	def initializeGL(self):
		# glewInit()
		GL.glClearColor(0.0, 0.0, 0.0, 0.0)
		GL.glEnable(GL.GL_LIGHTING)
		GL.glEnable(GL.GL_LIGHT0)
		GL.glEnable(GL.GL_NORMALIZE)
		GL.glEnable(GL.GL_COLOR_MATERIAL)
		GL.glEnable(GL.GL_CULL_FACE)
		GL.glEnable(GL.GL_DEPTH_TEST)


	def paintGL(self):
		grid_x = gl.GLGridItem()
		grid_x.scale(2, 2, 2)
		grid_x.rotate(90, 0, 1, 0)
		grid_x.translate(-5, 0, 0)
		self.plot.addItem(grid_x)

		grid_y = gl.GLGridItem()
		grid_y.scale(2, 2, 2)
		grid_y.rotate(90, 1, 0, 0)
		grid_y.translate(0, -5, 0)
		self.plot.addItem(grid_y)

		grid_z = gl.GLGridItem()
		grid_z.scale(2, 2, 2)
		grid_z.translate(0, 0, -5)
		self.plot.addItem(grid_z)

		vertsOfSphere1 = init.sphere1
		vertsOfSphere2 = init.sphere2
		# create the faces and colorsOfSphere1 arrays
		colorsOfSphere1 = []
		for p in init.env.Sphere1:
			colorsOfSphere1.append(p.color)
		colorsOfSphere1 = np.array(colorsOfSphere1)

		colorsOfSphere2 = []
		for p in init.env.Sphere2:
			colorsOfSphere2.append(p.color)
		colorsOfSphere2 = np.array(colorsOfSphere2)

		self.sphere1 = gl.GLScatterPlotItem(pos = vertsOfSphere1, color = colorsOfSphere1)
		self.sphere2 = gl.GLScatterPlotItem(pos = vertsOfSphere2, color = colorsOfSphere2)

		self.plot.addItem(self.sphere1)
		self.plot.addItem(self.sphere2)

	def update(self):
		if (self.run):
			if self.checkBox1.isChecked():
				init.applyExtForce2Sphere1(init.env)
			# else:
			# 	for p in (init.env.Sphere1):
			# 		p.ExtForce = np.asarray([0, 0, 0], dtype=float)
			if self.checkBox2.isChecked():
				init.applyExtForce2Sphere2(init.env)
			# else:
			# 	for p in (init.env.Sphere2):
			# 		p.ExtForce = np.asarray([0, 0, 0], dtype=float)
			if self.checkBox3.isChecked():
				self.drawForceArrows()
				# self.updateForceArrows()

			init.env.update_physics()

			vertsOfSphere1 = np.array(init.env.extract_vertsOfSphere1())
			colorsOfSphere1 = np.array(init.env.extract_colorsOfSphere1())
			self.sphere1.setData(pos=vertsOfSphere1, color=colorsOfSphere1)

			vertsOfSphere2 = np.array(init.env.extract_vertsOfSphere2())
			# colorsOfSphere2 = np.array(init.env.extract_colorsOfSphere2())
			colorsOfSphere2 = []
			for p in init.env.Sphere2:
				colorsOfSphere2.append([1, 0, 0, 1])
			colorsOfSphere2 = np.array(colorsOfSphere2)
			self.sphere2.setData(pos=vertsOfSphere2, color=colorsOfSphere2)


	def start(self):
		if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
			QtWidgets.QApplication.instance().exec_()

	def run2true(self):
		self.run = True

	def drawForceArrows(self):
		self.Forcelines = []
		for p in init.env.particles:
			xx = p.X[0]
			yx = p.X[1]
			zx = p.X[2]

			xy = p.X[0] + (p.ExtForce[0]/np.linalg.norm(p.ExtForce))
			yy = p.X[1] + (p.ExtForce[1]/np.linalg.norm(p.ExtForce))
			zy = p.X[2] + (p.ExtForce[2]/np.linalg.norm(p.ExtForce))
			Xdot = np.array([xx, yx, zx])
			Ydot = np.array([xy, yy, zy])
			pts = np.array([Xdot, Ydot])
			self.Forcelines.append(gl.GLLinePlotItem(pos=pts, width=1, antialias=False, color=[0, 0, 1, 1]))
		for l in self.Forcelines:
			self.plot.addItem(l)


	def updateForceArrows(self):
		for l in self.Forcelines:
			self.plot.removeItem(l)


	def animation(self):
		timer = QtCore.QTimer()
		timer.timeout.connect(self.update)
		timer.start(10)
		self.start()
		# self.update()

if __name__ == '__main__':
	app = QtWidgets.QApplication(sys.argv)
	t = GLWidget()
	t.animation()
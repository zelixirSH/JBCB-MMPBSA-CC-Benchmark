"""
amd.py: Implements the aMD integration method.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Peter Eastman, Steffen Lindert
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from openmm import CustomIntegrator
from openmm.unit import *
import numpy as np


class AMDGroupLIntegrator(CustomIntegrator):

    def __init__(self, temperature, friction, dt, groups1, alphas, Emaxs, groups2=[]):

        # =============== prepare ================
        self._temperature = temperature
        if is_quantity(temperature):
            self._temperature = temperature / kelvin
        kT = 0.0083144598 * self._temperature
        self._friction = friction

        self._dt = dt
        if is_quantity(dt):
            self._dt = self._dt/picosecond

        vscale = np.exp(-self._friction*self._dt);
        fscale = self._dt if self._friction == 0 else (1.0-vscale)/self._friction
        noisescale = np.sqrt(kT*(1.0-vscale*vscale))

        numG1 = len(groups1)
        numG2 = len(groups2)

        # =============== setup ================
        CustomIntegrator.__init__(self, dt)

        self.addGlobalVariable("vscale", vscale)
        self.addGlobalVariable("fscale", fscale)
        self.addGlobalVariable("noisescale", noisescale)

        self.addGlobalVariable("alpha", 0)
        self.addGlobalVariable("Emax", 0)
        self.addGlobalVariable("Eg", 0)
        self.addPerDofVariable("fg", 0)
        self.addPerDofVariable("oldx", 0) # position before application of constraints
        #self.addPerDofVariable("im", 0)
        #self.addPerDofVariable("sim", 0)

        self.addUpdateContextState();

        #self.addComputePerDof("im", "1.0/m")
        #self.addComputePerDof("sim", "sqrt(im)")
        #self.addComputePerDof("v", "vscale*v + noisescale*gaussian*sim")
        self.addComputePerDof("v", "vscale*v + noisescale*gaussian/sqrt(m)")

        # calc the scaled force group by aMD
        for i in range(numG1):
            g = str(groups1[i])
            self.addComputePerDof("fg", "f"+g)
            self.addComputeGlobal("Eg", "energy"+g)
            self.addComputeGlobal("alpha", str(alphas[i]))
            self.addComputeGlobal("Emax", str(Emaxs[i]))
            self.addComputePerDof("v", """v+fscale*fg*scale/m;
                                          scale=1-s+s*(alpha/(alpha+dE))^2;
                                          s=step(dE); dE=Emax-Eg""")
        # calc the normal force
        if numG2 > 0 :
            for g in groups2:
                self.addComputePerDof("fg", "f"+str(g))
                self.addComputePerDof("v", "v+fscale*fg/m")

        # update pos and speed
        self.addComputePerDof("oldx", "x")
        self.addComputePerDof("x", "x+dt*v")
        self.addConstrainPositions()
        self.addComputePerDof("v", "(x-oldx)/dt")

    def setTemperature(self, temperature):

        if is_quantity(temperature):
            self._temperature = temperature / kelvin
        kT = 0.0083144598 * self._temperature

        vscale = np.exp(-self._friction*self._dt);
        noisescale = np.sqrt(kT*(1-vscale*vscale))
        self.setGlobalVariableByName("noisescale", noisescale)


class IaMDLIntegrator(CustomIntegrator):

    def __init__(self, temperature, friction, dt, group, a, E, M, group2=-1):
     
        # =============== prepare ================
        self._temperature = temperature
        if is_quantity(temperature):
            self._temperature = temperature / kelvin
        kT = 0.0083144598 * self._temperature
        self._friction = friction

        self._dt = dt

        if is_quantity(dt):
            self._dt = self._dt/picosecond

        vscale = np.exp(-self._friction*self._dt); 
        fscale = self._dt if self._friction == 0 else (1.0-vscale)/self._friction
        noisescale = np.sqrt(kT*(1.0-vscale*vscale))

        num = len(a)
        an = []
        En = []
        Mn = []
        for i in range(num):
            an.append("a"+str(i+1))
            En.append("E"+str(i+1))
            Mn.append("M"+str(i+1))

        CustomIntegrator.__init__(self, dt)

        self.addGlobalVariable("vscale", vscale)
        self.addGlobalVariable("fscale", fscale)
        self.addGlobalVariable("noisescale", noisescale)
        self.addGlobalVariable("nbeta", -1.0/kT)
        self.addGlobalVariable("dE", 0.0)
        self.addGlobalVariable("scale", 0.0)
        self.addGlobalVariable("a", 0.0)
        self.addGlobalVariable("E", 0.0)
        self.addGlobalVariable("M", 0.0)
        for i in range(num):
            self.addGlobalVariable(an[i], a[i])
            self.addGlobalVariable(En[i], E[i])
            self.addGlobalVariable(Mn[i], M[i])

        self.addGlobalVariable("Z", 0)
        self.addGlobalVariable("K", 0)
        self.addGlobalVariable("EXP", 0)

        self.addGlobalVariable("e", 0)
        self.addPerDofVariable("oldx", 0)
        self.addPerDofVariable("fg", 0)

        self.addUpdateContextState()

        self.addComputePerDof("v", "vscale*v + noisescale*gaussian/sqrt(m)")
        self.addComputePerDof("fg", "f"+str(group))
        self.addComputeGlobal("e", "energy"+str(group))

        self.addComputeGlobal("K", "0.0")
        self.addComputeGlobal("Z", "0.0")

        for i in range(num):
            self.addComputeGlobal("a", an[i])
            self.addComputeGlobal("E", En[i])
            self.addComputeGlobal("M", Mn[i])
            self.addComputeGlobal("dE", "step(E-e)*(E-e)^2/(a+E-e)")
            self.addComputeGlobal("scale", "(a/(a+step(E-e)*(E-e)))^2")

            self.addComputeGlobal("EXP", "exp(nbeta*(dE+M))")
            self.addComputeGlobal("K", "K+scale*EXP")
            self.addComputeGlobal("Z", "Z+EXP")
        self.addComputeGlobal("K", "K/Z")
        self.addComputePerDof("v", "v+fscale*fg*K/m")

        # calc normal force
        if group2 >= 0:
            self.addComputePerDof("fg", "f"+str(group2))
            self.addComputePerDof("v", "v+fscale*fg/m")

        # update pos and speed
        self.addComputePerDof("oldx", "x")
        self.addComputePerDof("x", "x+dt*v")
        self.addConstrainPositions()
        self.addComputePerDof("v", "(x-oldx)/dt")

    def setTemperature(self, temperature):

        if is_quantity(temperature):
            self._temperature = temperature / kelvin
        kT = 0.0083144598 * self._temperature

        vscale = np.exp(-self._friction*self._dt); 
        noisescale = np.sqrt(kT*(1-vscale*vscale))
        self.setGlobalVariableByName("noisescale", noisescale)

    def setM(self, M):
        for i in range(len(M)):
            self.setGlobalVariableByName("M"+str(i+1), M[i])


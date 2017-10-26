################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#               2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
################################################################################

import h5py
from numpy import *
import sys
import setups

# Generates a swift IC file for the Sedov blast test in a periodic cubic box

# Parameters
gamma = 5./3.      # Gas adiabatic index
rho0 = 1.          # Background density
P0 = 1.e-6         # Background pressure
E0= 1.             # Energy of the explosion
N_inject = 15      # Number of particles in which to inject energy
fileName = "sedov.hdf5"

if len(sys.argv) < 2:
  print "No number of particles provided!"
  exit()

numPart_request = int(sys.argv[1])

if not numPart_request in setups.setups:
  print "No setup found for requested number of particles!"
  exit()

num_copy = setups.setups[numPart_request]["num_copy"]
glass_name = "glassCube_{0}.hdf5".format(
               setups.setups[numPart_request]["glass_name"])

#---------------------------------------------------
glass = h5py.File(glass_name, "r")

# Read particle positions and h from the glass
pos = glass["/PartType0/Coordinates"][:,:]
h = glass["/PartType0/SmoothingLength"][:] * 0.3

# Make copies
pos /= num_copy
h /= num_copy

# x direction
pos_all = pos
h_all = h
for i in range(1, num_copy):
  pos_all = append(pos_all, pos + array([i * 1. / num_copy, 0., 0.]), axis=0)
  h_all = append(h_all, h)

# y direction
pos = pos_all
h = h_all
for i in range(1, num_copy):
  pos_all = append(pos_all, pos + array([0., i * 1. / num_copy, 0.]), axis=0)
  h_all = append(h_all, h)

# z direction
pos = pos_all
h = h_all
for i in range(1, num_copy):
  pos_all = append(pos_all, pos + array([0., 0., i * 1. / num_copy]), axis=0)
  h_all = append(h_all, h)

pos = pos_all
h = h_all

numPart = size(h)
vol = 1.

# Generate extra arrays
v = zeros((numPart, 3))
ids = linspace(1, numPart, numPart)
m = zeros(numPart)
u = zeros(numPart)
r = zeros(numPart)

r = sqrt((pos[:,0] - 0.5)**2 + (pos[:,1] - 0.5)**2 + (pos[:,2] - 0.5)**2)
m[:] = rho0 * vol / numPart    
u[:] = P0 / (rho0 * (gamma - 1))

# Make the central particles detonate
index = argsort(r)
u[index[0:N_inject]] = E0 / (N_inject * m[0])

#--------------------------------------------------

#File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [1., 1., 1.]
grp.attrs["NumPart_Total"] =  [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

#Runtime parameters
grp = file.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = 1

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.
grp.attrs["Unit mass in cgs (U_M)"] = 1.
grp.attrs["Unit time in cgs (U_t)"] = 1.
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

#Particle group
grp = file.create_group("/PartType0")
grp.create_dataset('Coordinates', data=pos, dtype='d')
grp.create_dataset('Velocities', data=v, dtype='f')
grp.create_dataset('Masses', data=m, dtype='f')
grp.create_dataset('SmoothingLength', data=h, dtype='f')
grp.create_dataset('InternalEnergy', data=u, dtype='f')
grp.create_dataset('ParticleIDs', data=ids, dtype='L')

file.close()

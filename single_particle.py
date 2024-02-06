import numpy as numpy
from particle_system import particle, hamiltonian
from constants import *

class single_particle(particle, hamiltonian):

    def single_particle(particle):
        particle_mass = particle.m
        particle_spin = particle.s
        particle_charge = particle.c

    def expectation(single_particle, hamiltonian):
        if hamiltonian.spacial_dim == 1:
            ...
        if hamiltonian.spacial_dim == 2:
            ...
        if hamiltonian.spacial.dim == 3:
            ...



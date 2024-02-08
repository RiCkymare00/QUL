import numpy as numpy
from particle_system import *
from constants import *

class single_particle:
    def __init__ (self, particle_type, x_values, problema):
        particle = wave_function(particle_type, x_values, problema)
        particle_mass = particle.m
        particle_spin = particle.s
        particle_charge = particle.c
        print(particle_mass,particle_spin,particle_charge)
        particle_evolution = particle.time_evolution.evolution
        #expectation_energy = particle.H
        #print(expectation_energy)

        



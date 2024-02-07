from single_particle import *
import numpy as np

particle_type = input("Inserisci il tipo di particella: ")
problema = input("Inserisci il tipo di problema: ")
x_values = np.linspace(0, 10,10)  
particella = single_particle(particle_type, x_values, problema)

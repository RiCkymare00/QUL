import numpy as np
from constants import *
import time
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.linalg import expm
from scipy.linalg import eigh
import sympy as sp

global Psi_nuovo 
Psi_nuovo = np.zeros(100)

def potencial(x_values, delta_t,t):
    V = np.zeros_like(x_values)
    velocity = len(x_values) / delta_t
    if t > 10:
        buca = len(x_values) - int(velocity*(t - 10))
        V[x_values > buca] = 0.1
        if t > 15 :
            buca = len(x_values) - int(velocity*(t - 15))
            V[x_values > buca] = 0
    #print(V)
    return V 

def solve_eigenproblem(K,V):
    H = np.diag(K) + np.diag(V)
    _, autovettori = sp.Matrix(H).diagonalize()
    autostati = []
    Psi_normalized = []
    num_righe, num_colonne = autovettori.shape
    for i in range(num_righe):
        for j in range(num_colonne):
            elemento = autovettori[i, j]
            if elemento != 0:
                autostati.append(elemento)
    autostati_array = np.array(autostati)
    Psi_normalized = autostati_array/((np.sum(autostati_array**2))**0.5)
    return Psi_normalized
class Time_evolution:
    @staticmethod

    def evolution(psi, K, x_values, frames=100, interval=100):
        initial_time = 0
        simulation_time = 120
        t_values = np.linspace(initial_time, simulation_time, frames)
        fig, ax = plt.subplots()
        line, = ax.plot([], [])
        if not np.array_equal(Psi_nuovo, np.zeros_like(Psi_nuovo)):
            psi = Psi_nuovo
        
        def update(frame,x_values, K, frames, psi):
            t = t_values[frame]
            operatore = []
            operatore = np.exp(-1j * ((K + potencial(x_values, frames, t)) * t / hbar))
            if t > 10:
                kinetic = K 
                poten = potencial(x_values, frames, t)
                for x in range(int(t)+1, int(t)+2):
                    for y in range(2):
                        psi = psi - solve_eigenproblem(kinetic,potencial(x_values, frames, t+x))*psi-solve_eigenproblem(kinetic,poten)*psi-solve_eigenproblem(kinetic,potencial(x_values, frames, t-y))*psi
                
            evolved_psi = psi * operatore
            Psi = np.abs(evolved_psi)**2
            Psi_array = np.array(Psi)  
            Psi_normalized = Psi_array/((np.sum(Psi_array**2))**0.5)
            line.set_data(x_values, Psi_normalized)
            Psi_nuovo = Psi_normalized
            ax.set_title(f'Tempo: {t}')
            ax.figure.canvas.draw()
            return line, 

        x_min = min(x_values)
        x_max = max(x_values)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(-0.1,0.5)

        ani = FuncAnimation(fig, update, frames=frames, fargs=(x_values, K , frames, psi), blit=True, interval=interval)
        #ani.save('evoluzione_temporale.mp4', writer='ffmpeg')
        plt.show()

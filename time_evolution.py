import numpy as np
from constants import *
import time
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def potencial(x_values, delta_t,t):
    V = np.zeros_like(x_values)
    velocity = len(x_values) / delta_t
    if t > 50 and t < 75 :
        buca = len(x_values) - int(velocity*(t - 50))
        V[x_values > buca] = - float ('inf')
        #print(V)
    if t > 75:
        buca = len(x_values) - int(velocity*(t - 50))
        V[x_values > buca] = 0
        #print(V)
    return V 

class Time_evolution:
    @staticmethod
    def evolution(psi, K, x_values, frames=100, interval=50):
        initial_time = 0
        simulation_time = 100
        t_values = np.linspace(initial_time, simulation_time, frames)
        fig, ax = plt.subplots()
        line, = ax.plot([], [])

        def update(frame,x_values, K, frames, psi):
            t = t_values[frame]
            operatore = np.exp(-1j * (K + potencial(x_values,frames,t)) * t / hbar)
            evolved_psi = psi * operatore
            Psi = np.abs(evolved_psi)**2
            norm_Psi = np.linalg.norm(Psi)
            Psi_normalized = Psi/norm_Psi  
            line.set_data(x_values, Psi_normalized)
            ax.set_title(f'Tempo: {t}')
            return line,

        x_min = min(x_values)
        x_max = max(x_values)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(0,0.3) 

        ani = FuncAnimation(fig, update, frames=frames, fargs=(x_values, K , frames, psi), blit=True, interval=interval)
        ani.save('evoluzione_temporale.mp4', writer='ffmpeg')
        #plt.show()

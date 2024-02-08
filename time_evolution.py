import numpy as np
from constants import *
import time
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Time_evolution:
    @staticmethod
    def evolution(psi, K, x_values, frames=100, interval=50):
        initial_time = 0
        simulation_time = 30
        t_values = np.linspace(initial_time, simulation_time, frames)

        fig, ax = plt.subplots()
        line, = ax.plot([], [])

        def update(frame, x_values, K):
            t = t_values[frame]
            if t < 10 :
                potencial = 0
            elif t > 10: 
                potencial = -20000
            operatore = np.exp(-1j * (K(x_values) + potencial) * t / hbar)
            evolved_psi = psi * operatore  
            line.set_data(x_values, np.abs(evolved_psi)**2)
            ax.set_title(f'Tempo: {t}')
            return line,

        x_min = min(x_values)
        x_max = max(x_values)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(0, 0.001) 

        ani = FuncAnimation(fig, update, frames=frames, fargs=(x_values,K), blit=True, interval=interval)
        ani.save('evoluzione_temporale.mp4', writer='ffmpeg')
        #plt.show()

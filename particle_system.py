from constants import *
import numpy as np 
import matplotlib.pyplot as plt
import scipy.special
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time_evolution import Time_evolution


def V (x_values):
     Energy = np.zeros_like(x_values)
     return Energy

def expectation_value(psi,hamiltonian,x_values):
    integral = 0
    psi_dagger = np.conj(psi).T
    Energy = np.dot(np.dot(psi_dagger, hamiltonian), psi)
    for i in range (1,len(x_values)):
        interval_width = x_values[i] - x_values[i-1]
        integral += Energy[i] * interval_width
    return integral

def norm(psi,x_values):
    psi_dagger = np.conj(psi).T
    val = np.dot(psi_dagger, psi)
    normal = val*len(x_values)
    return normal

def radial_wave_function(r, n, l):

    rho = 2 * r / n
    radial_factor = np.sqrt((2 / (n**3 * np.math.factorial(n))) * np.exp(-rho) * rho**l)
    laguerre = np.polyval(np.polyder(np.poly1d([(-1)**l, 2*l+1, -n+l, n**2 * (n+1)])), rho)
    R = radial_factor * laguerre * np.exp(-rho / 2)
    return R

def angular_wave_function(theta, phi, l, m):

    Y = scipy.special.sph_harm(m, l, phi, theta)
    return Y

class wave_function:
    def __init__(self, particle_type, x_values, problema):

        if particle_type == 'elettrone':
            self.m = m_e
            self.s = 1/2
            self.c = e 
        #elif particle_type == 'altra particella':
        #    ...
        if problema == 'buca di potenziale':
            n = int(input("Definire il primo numero quantico: "))
            psi = []
            function = []
            for x in x_values:
                function.append(np.sqrt(2/len(x_values))*np.sin(n*x*np.pi/len(x_values)))
                psi.append((np.abs(np.sqrt(2/len(x_values))*np.sin(n*x*np.pi/len(x_values))))**2)
            norma_psi = np.linalg.norm(psi)
            psi_normalizzato = psi / norma_psi
            #self.n = norm(psi,x_values)
            plt.style.use("dark_background")
            plt.plot(x_values,psi_normalizzato, linewidth=0.5, label="Wave function")
            plt.show()
          
        elif problema == 'atomo idrogeno':
            n = int(input("Definire il primo numero quantico: "))
            l = int(input("Definire il secondo numero quantico: "))
            m = int(input("Definire il terzo numero quantico: "))
            if l > (n-1) :
                l = input("Il secondo numero inserito non è ammesso, definire nuovamente il secondo numero quantico: ") 
            if m < -l or m > l:
                m = input("Il terzo numero inserito non è ammesso, definire nuovamente il terzo numero quantico: ")
            psi = []
            function = []
            X_values = []
            y_values = []
            z_values = []
            theta_values = np.linspace(0, np.pi, 100)
            phi_values = np.linspace(0, 2*np.pi, 100)
            for r in x_values:
                for theta in theta_values:
                    for phi in phi_values:
                        R = radial_wave_function(r, n, l)
                        Y = angular_wave_function(theta, phi, l, m)
                        function.append(R * Y)
                        psi.append(np.abs(R * Y)**2)
                        x = np.abs(R * Y)**2 * np.sin(theta) * np.cos(phi)
                        y = np.abs(R * Y)**2 * np.sin(theta) * np.sin(phi)
                        z = np.abs(R * Y)**2 * np.cos(theta)
                        X_values.append(x)
                        y_values.append(y)
                        z_values.append(z)         
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(X_values, y_values, z_values, c=psi, cmap='viridis', alpha=0.5)
            ax.set_title('$\psi$ atomo idrogeno; n: {}, l: {}, m: {}'.format(n, l, m))
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            plt.show()
    
        ''' Si richiede all'utente di definire il numero di dimensioni ed il tipo di sistema che si vuole studiare
        '''
        if problema == 'buca di potenziale':
            H_kinetic = -hbar/2 * m_e * np.gradient(np.gradient(function, x_values), x_values)
            H = H_kinetic + V(x_values)
            #self.H = expectation_value(psi,H,x_values)
            self.time_evolution = Time_evolution()
            self.time_evolution.evolution(psi_normalizzato,H_kinetic,x_values,frames=100, interval=100)
    
        elif problema == 'atomo idrogeno':
            H_kinetic = -hbar/2 * m_e * np.gradient(np.gradient(function, x_values), x_values)
            H = H_kinetic + V(x_values)
            self.H = expectation_value(psi,H,x_values)
            self.time_evolution = Time_evolution() 
            self.time_evolution.evolution(function,H_kinetic,x_values,frames=100, interval=100)
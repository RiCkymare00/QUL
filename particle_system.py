from constants import *
import numpy as np 
import matplotlib.pyplot as plt
import scipy.special
import plotly.graph_objects as go

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
    """Calcola la funzione d'onda radiale dell'atomo di idrogeno per un valore specifico di r."""
    rho = 2 * r / n
    radial_factor = np.sqrt((2 / (n**3 * np.math.factorial(n))) * np.exp(-rho) * rho**l)
    laguerre = np.polyval(np.polyder(np.poly1d([(-1)**l, 2*l+1, -n+l, n**2 * (n+1)])), rho)
    R = radial_factor * laguerre * np.exp(-rho / 2)
    return R

def angular_wave_function(theta, phi, l, m):
    """Calcola la funzione d'onda angolare dell'atomo di idrogeno per valori specifici di theta e phi."""
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
            for n in range (1,10):
                psi = []
                for x in x_values:
                    psi.append((np.abs(np.sqrt(2/len(x_values))*np.sin(n*x*np.pi/len(x_values))))**2)
                #self.n = norm(psi,x_values)
                plt.style.use("dark_background")
                plt.plot(x_values,psi, linewidth=0.5, label="Wave function")
            plt.show()
          
        elif problema == 'atomo idrogeno':
            for n in range (1,2):
                l = n-1
                for m in range(-l, l+1):
                    psi = []
                    X_values = []
                    y_values = []
                    z_values = []
                    theta_values = np.linspace(0, np.pi/2, 10)
                    phi_values = np.linspace(0, np.pi, 10)
                    for r in x_values:
                        for theta in theta_values:
                            for phi in phi_values:
                                R = radial_wave_function(r, n, l)
                                Y = angular_wave_function(theta, phi, l, m)
                                psi.append(np.abs(R * Y)**2)
                                x = r * np.sin(theta) * np.cos(phi)
                                y = r * np.sin(theta) * np.sin(phi)
                                z = r * np.cos(theta)
                                X_values.append(x)
                                y_values.append(y)
                                z_values.append(z)         
            fig = go.Figure(data=go.Scatter3d(x=X_values, y=y_values, z=z_values, mode='markers', marker=dict(size=5, color=psi, colorscale='Viridis', opacity=0.8)))
            fig.update_layout(scene=dict(xaxis=dict(title='x'), yaxis=dict(title='y'), zaxis=dict(title='z')), title='Funzione d\'onda $\psi$ in coordinate cartesiane', autosize=True)
            fig.show()
            fig.write_image("funzione_onda_atomo_idrogeno.png")
    
        ''' Si richiede all'utente di definire il numero di dimensioni ed il tipo di sistema che si vuole studiare
        '''
        if problema == 'buca di potenziale':
            H_kinetic = -hbar/2 * m_e * np.gradient(np.gradient(psi, x_values), x_values)
            H = H_kinetic + V(x_values)
            self.H = expectation_value(psi,H,x_values)
    
        elif problema == 'atomo idrogeno':
           ...

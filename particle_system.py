from constants import *
import numpy as np 

def wave_function(x_values):
    '''Si richiede all'utente di definire il tipo di particella presente nel sistema
    '''
    if('''condizioni per la buca di potenziale'''):
        for n in range (0,100): #numeri quantici
            for x in x_values:
                psi = np.sqrt(2/len(x_values))*np.sin(n*x*np.pi/len(x_values))

    if('''condizioni per atomo di idrogeno'''):
        ...
    

def hamiltonian():
    ''' Si richiede all'utente di definire il numero di dimensioni ed il tipo di sistema che si vuole studiare
    '''
    if('''condizioni per la buca di potenziale 1D'''):
        x_values = np.linspace(-1, 1, 100) 
        H = -hbar/2*m_e*(np.gradient(np.gradient(wave_function(x_values),x_values),x_values),x_values) + V(x_values)

    if('''condizioni per atomo di idrogeno'''):
        ...
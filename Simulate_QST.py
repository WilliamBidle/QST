''' Code written by William Bidle '''

##########################################################################################
######################################## Imports #########################################
##########################################################################################

''' For Array Manipulation '''
import numpy as np

''' For Mathematical Operations '''
from scipy import special
from scipy.special import factorial

''' For Visualization '''
import matplotlib.pyplot as plt

''' For Random Number Generation '''
import random

''' For Progress Bar '''
from tqdm.notebook import tnrange 

##########################################################################################
####################################### Functions ########################################
##########################################################################################

def Quad_Fock_Projection(n, x):
    
    ''' 
    ##################################################################
                 Quadrature Projection of a Fock State 
    ##################################################################
    
    INPUTS:
    - n (int) : the Fock state under consideration
    - x (1D NumPy Array) : the list of possible quadratures
    
    OUTPUTS:
    - result (1D NumPy Array) : the quadratures corresponding to n_list
    
    ##################### Example ######################
    
    from Quadrature_Simulation_Methods import *
    
    n = 1
    x = np.linspace(-5,5,100)
    
    result = Quad_Fock_Projection(n, x)
    
    plt.plot(x, result)
    plt.show()
    
    '''
    
    Hermite_n = special.hermite(n)(x)
    result = (1/np.pi)**(1/4)/(np.sqrt(2**n*factorial(n)))*np.exp(-x**2/2)*Hermite_n
    return result

##########################################################################################

def generate_psi(n_list, coefs, x):
    
    ''' 
    ##################################################################
                         Generate a Wavefunction 
    ##################################################################
    
    INPUTS:
    - n_list (1D NumPy array) : the Fock states under consideration
    - coefs (1D NumPy array) : the coefficients corresponding to the above Fock states
    - x (1D NumPy array) : the possible quadrature values
    
    OUTPUTS:
    - psi_tot (2D NumPy Array) : the possible quadrature values per Fock state
    
    ##################### Example ######################
    
    from Quadrature_Simulation_Methods import *
    
    n_list = np.linspace(0, 10, 11)
    coefs = np.random.rand(10)
    coefs = coefs/np.linalg.norm(coefs) 
    x = np.linspace(-5,5,100)
    
    psi_tot = generate_psi(n_list, coefs, x)
    
    '''
        
    # check that the coefs are properly normalized
    tot = 0
    for i in coefs:
        tot += i*np.conj(i)
    
    # if the coefficients have't been properly normalized, normalize them
    if np.round(tot, 3) != 1.: 
        coefs = coefs/np.linalg.norm(coefs) 
    
    psi_tot = []
    for i in range(len(n_list)):
        psi_tot.append(coefs[i]*Quad_Fock_Projection(n_list[i], x))
        
    return psi_tot

##########################################################################################

def generate_coherent_state(n_list, avg_n):
    
    ''' 
    ##################################################################
                         Generate a Coherent State 
    ##################################################################
    
    INPUTS:
    - n_list (list) : the Fock states under consideration
    - avg_n (float) : the average photon number
    
    OUTPUTS:
    - coherent_state_coefs (1D NumPy Array) : the coherent state coefficients corresponding to n_list
    
    ##################### Example ######################
    
    from Quadrature_Simulation_Methods import *
    
    n_list = np.linspace(0, 10, 11)
    avg_n = 2
    
    coherent_state_coefs = coherent_state(n_list, avg_n)
    
    plt.plot(n_list, coherent_state_coefs)
    plt.show()
    
    '''
    
    n_list = np.array(n_list, dtype = float)
    coherent_state_coefs = np.array(np.sqrt(avg_n**n_list*np.exp(-avg_n)/factorial(n_list)))

    return coherent_state_coefs

##########################################################################################

def generate_quadratures(n_list, coefs, theta, x, visualize, save_data, filename = None):
    
    ''' 
    ##################################################################
                    Generate a Quadrature Dataset 
    ##################################################################
    
    INPUTS:
    - n_list (1D NumPy array) : the Fock states under consideration
    - coefs (1D NumPy array) : the coefficients corresponding to the above Fock states
    - theta (1D NumPy array) : the list of phase values
    - x (1D NumPy array) : the list of possible quadratures
    - visualize (Boolean) : whether or not to visualize the generated quadratures
    - save_data (Boolean) : whether or not to save the data to a file
    - filename (string) : create and save the data to 'filename' if save_data = True
    
    OUTPUTS:
    - theta (1D NumPy Array) : the coherent state coefficients corresponding to n_list
    - quadratures (1D NumPy Array) : the coherent state coefficients corresponding to n_list
    
    ##################### Example ######################
    
    For an example of this in action, see the file "Generate Quadratures Simulation.ipynb" 
    
    '''
        
    psi = generate_psi(n_list, coefs, x) # create a wavefunction based on the declared parameters
    
    
    quadratures = []
    y_vis = 0
    for i in tnrange(len(theta), desc = 'Quadrature Generation'):
    # for i in range(len(theta)):
        y = 0
        for j in range(len(psi)):

            val = np.exp(n_list[j]*1j*(theta[i] - np.pi/2))*psi[j]
            y += val

        psi_norm = y*np.conj(y)

        # random_sample = [x[np.where(psi_norm == max(psi_norm))[0][0]]]
        random_sample = random.choices(x, weights = psi_norm) # randomly select a value of the domain weighted by psi_norm
        quadratures += random_sample
        
    if visualize == True:
        fig, ax = plt.subplots(figsize = (10, 6))

        ax.scatter(theta, quadratures, s = 2) # plot of quadrature sampling
        ax.set_xlabel('Phase', fontsize = 18)
        ax.set_ylabel('Quadrature (arb.)', fontsize = 18)

        plt.show()
        
    if save_data == True:
        
        try:
            np.savez("Data/" + filename, theta = theta, x = quadratures) # saves the data in the MLE folder to be directly used
            print('data saved at ' + filename + '.npz')
            
        except:
            os.mkdir('Data')
            np.savez("Data/" + filename, theta = theta, x = quadratures) # saves the data in the MLE folder to be directly used
            print('data saved at Data/' + filename + '.npz')
        
    return theta, quadratures

##########################################################################################
##########################################################################################
##########################################################################################

''' Code written by William Bidle '''

##########################################################################################
######################################## Imports #########################################
##########################################################################################

''' For Array Manipulation '''
import numpy as np

''' For Mathematical Operations '''
import scipy
import scipy.integrate as integrate
import scipy.special as special
from scipy.optimize import curve_fit

''' For Visualization '''
import matplotlib.pyplot as plt 

''' For Random Number Generation '''
import random

''' For Symbolic Manipulation '''
from sympy import * 

''' For Runtime Checking '''
from datetime import datetime

''' For Wigner Function Visualization '''
from qutip import *

##########################################################################################
####################################### Functions ########################################
##########################################################################################

def coherent_state(n_list, avg_n):
    
    ''' 
    ##################################################################
        Definition of a Coherent State (mostly used for fitting)
    ##################################################################
    
    INPUTS:
    - n_list (list) : the Fock states under consideration
    - avg_n (float) : the average photon number
    
    OUTPUTS:
    - coherent_state_coefs (1D NumPy Array) : the coherent state coefficients corresponding to n_list
    
    ##################### Example ######################
    
    from MLE_Functions import *
    
    n_list = np.linspace(0, 10, 11)
    avg_n = 2
    
    coherent_state_coefs = coherent_state(n_list, avg_n)
    
    plt.plot(n_list, coherent_state_coefs)
    plt.show()
    
    '''
    
    coherent_state_coefs = (np.exp(-avg_n/2)*np.sqrt(avg_n)**n_list / np.sqrt(special.factorial(n_list)))**2
    return coherent_state_coefs

##########################################################################################

def calc_rho(dims, num_iters, theta_data, quadrature_data, visualize_steps):
    
    ''' 
    ##################################################################
    Calculate most likely density matrix given some quadrature data
    ##################################################################
    
    INPUTS:
    - dims (int) : the number of dimensions/Fock states to consider
    - num_iters (int) : the number of iterations to perform
    - theta_data (1D NumPy Array) : a list of the phase measurement associated to a quadrature data point
    - quadrature_data (1D NumPy Array) : a list of the quadrature data associated to a phase measurement
    - visualize_steps (Boolean) : Whether or not to visualize the denisty matrix results after each iteration
    
    OUTPUTS:
    - rho (2D NumPy Array) : the most likely density matrix for the given dataset
    - likelihood_trend (1D NumPy Array) : the likelihood trend to see how the density matrix converges over time
    
    ##################### Example ######################
    
    from MLE_Functions import *
    
    dims = 10 
    num_iters = 20 
    theta_data = np.linspace(0, 2*np.pi, 100)
    quadrature_data = np.random.randn(100)
    visualize_steps = True

    rho, likelihood_trend = calc_rho(dims, num_iters, theta_data, quadrature_data, visualize_steps)
    
    '''
    
    likelihood_trend = [] # will be used to keep track of the likelihood of rho at each step
    
    n_list = np.linspace(0, dims - 1, dims, dtype = int) # creates an integer list of the values of n (Fock states) we are considering

    rho = np.zeros((dims, dims), complex) # initialize the correct dimension denisty matrix
    np.fill_diagonal(rho, 1) # fill the diagonals with ones
    rho = rho/sum(rho.diagonal()) # normalize the trace of the density matrix

    ''' The below matricies allow us to do the evaluation of the n wavefunctions at every data point (theta, x) in a vectorized fashion '''
    N = np.array(np.split(np.repeat(n_list, len(quadrature_data)), dims)) # A matrix that repeats the Fock state list in each column
    THETA = np.tile(theta_data, (dims, 1)) # A matrix that repeats the theta data list in each column
    X = np.tile(quadrature_data, (dims, 1)) # A matrix that repeats the quadrature data list in each row

    ################################################### Scipy Variable and Function Declarations ######################################################

    x = Symbol('x') # declare variable 'x' signifying the quadrature value
    theta = Symbol('\\theta') # declare variable 'theta' signifying the value of the phase
    n = Symbol('n', integer=True) # declare variable 'n' signifying the Fock states
    
    # Note: The coefs need to be done separately from the functions since it was giving me trouble at one point...
    coefs = (1/np.pi)**(1/4)/sqrt(2**n*factorial(n)) # the sympy function for the coefficients
    g = lambdify(n, coefs)(N) # evaluate the coefficents at the values of N

    func = exp(1j*n*(theta - np.pi/2))*exp(-x**2/2)*hermite(n, x) # the sympy function for the coefficients
    f = lambdify((n, theta, x), func)(N, THETA, X) # evaluate the coefficents at the values of N
    
    ''' Debug ''' 
    # print('The Wavefunction Equation:')
    # display(coefs*func)
    # print()
    
    ####################################################### Generate Useful Matricies ##########################################################
    
    psi_matrix = (f*g).transpose() # evaluation of the wavefunction at each data point (theta, x) for every consiered n - rows are (theta, x), columns are n 
    psi_squared_matrix = np.einsum('bi,bo->boi', psi_matrix, psi_matrix.conj()) # complex outer product
    likelihood = np.real(np.sum(np.log(np.sum(np.sum(psi_squared_matrix*np.real(rho), axis = 1), axis = 1)))) # calculate the likelihood of the current rho
    likelihood_trend.append(likelihood) # add likelihood to trend list
    
    ''' Debug '''
    # print('Psi Matrix is:\n', psi_matrix.transpose())
    # print('Psi Squared Matrix is:\n', psi_squared_matrix)
    
    # print('Initial Density Matrix:\n', rho, '\n')
    
    # plt.bar(n_list, np.diagonal(rho))
    # plt.show()

    # print('Likelihood For Current Rho:', likelihood)

    ################################################################ Perform the MLE ###################################################################
    
    for i in range(num_iters): # repeat the likelihood num_iters times

        R_matrix = np.real(np.diag(1/np.sum(np.sum(psi_squared_matrix*np.real(rho), axis = 1), axis = 1))/len(quadrature_data)) # R matrix for current rho

        rho = psi_matrix.conj().T@R_matrix@psi_matrix@rho@psi_matrix.conj().T@R_matrix@psi_matrix # calculate the new rho
        rho = rho/sum(np.diagonal(np.real(rho)))# need to normalize the density matrix such that the trace is 1
        
        likelihood = np.real(np.sum(np.log(np.sum(np.sum(psi_squared_matrix*np.real(rho), axis = 1), axis = 1)))) # calculate the likelihood of the new rho
        likelihood_trend.append(likelihood) # add the new likelihood to trend list
        
        ''' Debug '''
        # print('R Matrix is:\n', R_matrix)
        # print('Current Denisty Matrix is:\n', rho)
        # print('Diagonal Elements of Density Matrix are:\n', np.diagonal(np.real(rho)))
        # print('Likelihood For Current Rho:', likelihood)
        
        if visualize_steps == True: # visualize how the diagonal elements of the density matrix evolve over time
            
            fig, ax = plt.subplots(figsize = (10, 6))
            ax.bar(n_list, np.diag(np.real(rho)), color = 'lightcoral', edgecolor='black', label = 'iteration %d' % i + 1)
            ax.set_xticks(n_list)
            
            ax.grid(linestyle = '--')
            ax.set_axisbelow(True)

            ax.set_xlabel(r'$n$', fontsize = 16)
            ax.set_ylabel(r'$\rho_{nn}$', fontsize = 16)
            
            ax.legend(fontsize = 14)
        
            plt.show()
    
    return rho, likelihood_trend

##########################################################################################

def generate_quadratures(rho, num_data_points):

    ''' 
    ##################################################################
    Generate a set of quadratures that are associated with a certain phase for the input density matrix (for computing errors)
    ##################################################################
    
    INPUTS:
    - rho (2D NumPy Array) : the input density matrix
    - num_data_points (int) : the number of data points to simulate
    
    OUTPUTS:
    - theta_data (1D NumPy Array) : the simulated list of phases
    - quadrature_data (1D NumPy Array) : the simulated list of quadratures associated with the above phases
    
    ##################### Example ######################
    
    from MLE_Functions import *
    
    rho = np.random.rand(2,2)
    rho /= np.sum(np.diag(rho))
    num_data_points = 100
    
    theta_data, quadrature_data = generate_quadratures(rho, num_data_points)
    
    plt.scatter(theta_data, quadrature_data, s = 5)
    plt.show()
    
    '''
    
    dims = len(rho)
    
    n_list = np.linspace(0, dims - 1, dims, dtype = int)
    theta_data = np.linspace(0, 2*np.pi, num_data_points)
    x_data = np.linspace(-5, 5, num_data_points)

    N = np.array(num_data_points*list(np.array([np.tile(n_list, (num_data_points, 1))]))) # A matrix that repeats the Fock state list in each column
    THETA = np.array(np.split(np.repeat(theta_data, num_data_points*dims), num_data_points)).reshape(num_data_points, num_data_points, dims) # A matrix that repeats the theta data list in each column
    X = np.array(num_data_points*list(np.array([np.array(np.split(np.repeat(x_data, dims), num_data_points))]))) # A matrix that repeats the quadrature data list in each row

    ''' Debug '''
    # print(N)
    # print(X)
    # print(THETA)

    ##################### Scipy Function Declaration #####################

    x = Symbol('x') # declare variable 'x' signifying the quadrature value
    theta = Symbol('\\theta') # declare variable 'theta' signifying the value of the phase
    n = Symbol('n', integer=True) # declare variable 'n' signifying the Fock states

    ''' Note: The coefs need to be done separately since they need scipy in the labdify function to do factorial of lists '''
    coefs = (1/np.pi)**(1/4)/sqrt(2**n*factorial(n)) # the sympy function for the coefficients
    g = lambdify(n, coefs)(N) # evaluate the coefficents at the values of N

    func = exp(1j*n*(theta - np.pi/2))*exp(-x**2/2)*hermite(n, x) # the sympy function for the coefficients
    f = lambdify((n, theta, x), func)(N, THETA, X) # evaluate the coefficents at the values of N


    Chi_matrix = (f*g) # evaluation of each data point for every consiered n - rows are (x, theta) values and columns are n values
    psi_squared_matrix = np.einsum('...bi,...bo->...bio', Chi_matrix, Chi_matrix.conj())
    result = np.sum(np.sum(psi_squared_matrix*rho, axis = 3), axis = 2)

    ''' Debug '''
    # print('The Wavefunction Equation:')
    # display(coefs*func)
    # print()

    quadrature_data = []
    for i in range(num_data_points):
        quadrature_data.append(random.choices(x_data, weights = result[i])[0])

    return theta_data, quadrature_data

##########################################################################################

''' WORK IN PROGRESS '''
def calc_err(data_length, num_err_iters, rho):
    
    ''' The master MLE function to calculate the most likely density matrix based off of the quadrature data as well as compute errors
    
    INPUTS:
    - dims (int) : the number of dimensions/Fock states to consider
    - num_iters (int) : the number of iterations to perform
    - theta_data (1D NumPy Array) : a list of the phase measurement associated to a quadrature data point
    - quadrature_data (1D NumPy Array) : a list of the quadrature data associated to a phase measurement
    - visualize_steps (Boolean) : Whether or not to visualize the denisty matrix results after each iteration
    
    OUTPUTS:
    - rho (2D NumPy Array) : the most likely density matrix for the given dataset
    - likelihood_trend (1D NumPy Array) : the likelihood trend to see how the density matrix converges over time
    
    ##################### Example ######################
    
    WORK IN PROGRESS
    
    '''
    
    diags = np.array(np.diagonal(np.real(rho)))
    diags_prime = []
    
    for i in range(num_err_iters):
        theta_data, quadrature_data = generate_quadratures(rho, data_length) # generate a quadrature dataset from the calculated rho

        ''' Debug '''
        # plt.scatter(theta_data, quadrature_data, s = 2)
        # plt.show()

        rho_prime, likelihood_trend = calc_rho(len(rho), num_err_iters, theta_data, quadrature_data, visualize_steps = False)
        diags_prime.append(np.array(np.diagonal(np.real(rho_prime)), dtype = float)) # get the diagonal of this rho_prime


    err = np.sum(np.abs(diags - diags_prime), axis = 0)/num_err_iters
    return err

##########################################################################################

def perform_MLE(theta_data, quadrature_data, dims, num_iters, include_errors, num_err_iters, is_coherent_state, show_wigner, visualize_steps = False):
    
    ''' The master MLE function to calculate the most likely density matrix based off of the quadrature data as well as compute errors
    
    INPUTS:
    - theta_data (1D NumPy array) : the phase information corresponding to the quadrature data
    - quadrature_data (1D NumPy array) : the given quadrature data to analyze
    - dims (int) : the number of dimensions/Fock states to consider
    - num_iters (int) : the number of iterations to perform
    - include_errors (Boolean) : whether or not to calculate errors (WORK IN PROGRESS)
    - num_err_iters (int) : the number of iterations to perform error averaging
    - is_coherent_state (Boolean) : whether or not the state under consideration is a coherent state (fits the result)
    - show_wigner (Boolean) : whether or not to show the wigner distribution
    - visualize_steps (Boolean) : Whether or not to visualize the denisty matrix diagonals after each iteration
    
    OUTPUTS:
    - final_rho (2D NumPy Array) : the most likely density matrix for the given dataset
    
    ##################### Example ######################
    
    For an example of this in action, see the file "Better MLE.ipynb"
    
    '''
    
    ####################################################### Calculate the Most Likely Rho #########################################################

    start = datetime.now() # Record current timestamp

    final_rho, likelihood_trend = calc_rho(dims, num_iters, theta_data, quadrature_data, visualize_steps)
    diags = np.array(np.diagonal(np.real(final_rho)), dtype = float)

    end = datetime.now() # record loop end timestamp
    td = (end - start).total_seconds() * 10**3 # find difference loop start and end time and display
    print(f"The time to calculate the most likely density matrix was : {td:.03f}ms")

    ############################################################### Visualization #################################################################

    if include_errors == True:

        print()
        print("Calculating Errors...")
        print()

        start = datetime.now() # Record current timestamp

        err = calc_err(len(theta_data), num_err_iters, final_rho)

        end = datetime.now() # record loop end timestamp
        td = (end - start).total_seconds() * 10**3 # find difference loop start and end time and display
        print(f"The time to calculate the errors was : {td:.03f}ms")
        ############################################################## Plotting With Errors ###############################################################

        n_list = np.linspace(0, dims - 1, dims)
        avg_n = np.average(n_list, weights=diags)


        fig, ax = plt.subplots(figsize = (10, 6))
        ax.bar(n_list, diags, color = 'lightcoral', edgecolor='black')
        ax.errorbar(n_list, diags, yerr = err, capsize = 5, fmt = '.k')
        ax.set_xticks(n_list)

        ax.grid(linestyle = '--')
        ax.set_axisbelow(True)

        ax.set_xlabel(r'$n$', fontsize = 16)
        ax.set_ylabel(r'$\rho_{nn}$', fontsize = 16)

        plt.show()

    else:
        n_list = np.linspace(0, dims - 1, dims)
        avg_n = np.average(n_list, weights=diags)

        fig, ax = plt.subplots(figsize = (10, 6))
        ax.bar(n_list, diags, color = 'lightcoral', edgecolor='black')
        ax.set_xticks(n_list)

        ax.grid(linestyle = '--')
        ax.set_axisbelow(True)

        ax.set_xlabel(r'$n$', fontsize = 16)
        ax.set_ylabel(r'$\rho_{nn}$', fontsize = 16)

        if is_coherent_state == True: # include a fit of the coherent state to determine average photon number
            popt, pcov = curve_fit(coherent_state, n_list, diags, 
                                       p0 = [avg_n], 
                                       bounds=([0],[np.inf]))

            z = np.linspace(0, n_list[-1], 101)
            ax.plot(z, coherent_state(z, *popt), color = 'black', linestyle = '-.', label = r'Fit: $\alpha$ = %s' % np.round(np.sqrt(popt[0]), 3))

            ax.legend(fontsize = 14)

        plt.show()
        
    if show_wigner == True:

        print('The Wigner Distribution:\n')
        q_rho = Qobj(final_rho)

        xvec = np.linspace(-10, 10, 100)

        W = wigner(q_rho, xvec, xvec)
        cmap = wigner_cmap(W)
        X, Y = np.meshgrid(xvec, xvec)

        fig = plt.figure(figsize = (10, 10))
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_surface(X, Y, W, cmap="viridis", lw=0, rstride=1, cstride=1)
        ax.set_xlabel('B', fontsize = 22)
        ax.set_ylabel('E', fontsize = 22)
        ax.set_zlabel('W(E,B)', fontsize = 22)

        plt.show()
        
    return final_rho

##########################################################################################
##########################################################################################
##########################################################################################

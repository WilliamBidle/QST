# Quantum State Estimation through Maximum Likelihood Estimation

This work is a culmanation of my graduate school masters project at Stony Brook University in the [Quantum Information Science & Technology Group](http://qit.physics.sunysb.edu/wordpress/) under Dr. Eden Figueroa to perform quantum state tomography through maximum likelihood estimation. The main files to perform the simulation of quantum states ($\textit{Quadrature_Simulation_Methods.py}$) as well as to perform the maximum likelihood estimation of quantum states ($\textit{MLE.py}$) are described below. 

### Quadrature_Simulation_Methods.py (Quadrature Generation)

Serves to simulate real world examples of quantum states that can be measured through techniques such as Homodyne Detection.

It provides:

- An efficient way of simulating quadrature data for a given quantum state
- Can work with any arbitrary quantum state, including coherent states and (WORK IN PROGRESS) squeezed states

Usage:

The Jupyter Notebook $\textit{Generate Quadratures Simulation.ipynb}$ contains several working examples of how to easily generate different sets of quadrature data. Datasets from this notebook can be easily saved to the $\textit{Data}$ folder to be used for Maximum Likelihood Reconstruction. 

In order to use the functions in a different Python file, make sure to download the $\textit{Quadrature_Simulation_Methods.py}$ file to the same folder as the new code and use the command:

    from Quadrature_Simulation_Methods.py import *

The full list of functions, their usage, as well as some examples can be found within the above Python file.

### MLE.py (Maximum Likelihood Estimation)

Serves to reconstruct the density matrix of a given quantum state, and can be used to extract amplitude and phase information. Our group has shown it to work with real data gathered in the lab through Homodyne Detection of several different coherent states.

It provides:

- A vectorized method of finding the most likely density matrix corresponding to a given set of quadrature data
- Visualizations of the Wigner distribution of the quantum state through [QuTiP](https://qutip.org/docs/4.0.2/guide/guide-visualization.html)

Usage:

The Jupyter Notebook $\textit{MLE.ipynb}$ contains an illustritive working example of how to reconstruct a density matrix given a set of quadrature data. It can easily be used in conjuction with the $\textit{Generate Quadratures Simulation.ipynb}$. 

In order to use the functions in a different Python file, make sure to download the $\textit{MLE.py}$ file to the same folder as the new code and use the command:

    from MLE.py import *
    
The full list of functions, their usage, as well as some examples can be found within the above Python file.

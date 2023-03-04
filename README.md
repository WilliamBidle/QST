# QST

QST is a public Python package for performing quantum state tomography of quantum states through maximum likelihood estimation. This work is a culmanation of my graduate school masters project at Stony Brook University in the [Quantum Information Science & Technology Group](http://qit.physics.sunysb.edu/wordpress/) under Dr. Eden Figueroa. You can read more about the work in my included masters thesis, ***Quantum_State_Tomography.pdf***. 

### Installation

QST is avaliable on the Python Package Index [PyPI](https://pypi.org/project/QST/) and can be installed through pip:

    pip install QST

QST contains two major function libraries for performing quantum state tomography: QSim for simulations and MLE for maximum likelihood estimation.

### QST.QSim (Simulation of Quantum Statistics)

Serves to simulate real world examples of quantum states that can be measured through techniques such as Optical Homodyne Detection.

It provides:

- An efficient way of simulating quadrature data for a given quantum state
- Can work with any arbitrary quantum superposition, including coherent states and (WORK IN PROGRESS) squeezed states

Usage:

The Jupyter Notebook ***QST_Simulation.ipynb*** in the ***examples*** folder contains several working examples of how to easily generate different sets of quadrature data from different quantum states. Datasets from this notebook can be saved to the ***Data/*** subfolder to be used for maximum likelihood reconstruction. The full list of functions, their usage, as well as some examples can be found within the ***QSim.py*** file under ***src/QST/***.

### QST.MLE (Maximum Likelihood Estimation)

Serves to reconstruct the density matrix of a given quantum state, and can be used to extract amplitude and phase information. Our group has shown it to work with real data gathered in the lab through Homodyne Detection of several different coherent states.

It provides:

- A vectorized method of finding the most likely density matrix corresponding to a given set of quadrature data
- Visualizations of the Wigner distribution of the quantum state through [QuTiP](https://qutip.org/docs/4.0.2/guide/guide-visualization.html)
- Relative phase information between different quantum states

Usage:

The Jupyter Notebook ***MLE.ipynb*** in the ***examples*** folder contains an illustritive working example of how to reconstruct a density matrix given a set of quadrature data. It can easily be used in conjuction with the ***Generate_Quadratures_Simulation.ipynb***. The full list of functions, their usage, as well as some examples can be found within the ***MLE.py*** file under ***src/QST/***.

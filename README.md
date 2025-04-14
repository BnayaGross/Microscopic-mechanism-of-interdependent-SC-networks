# Microscopic-mechanism-of-interdependent-SC-networks
The files in this folder refer to the microscopic mechanism of interdependent superconducting networks studied in the article entitled "The microscopic origin of abrupt transitions in interdependent systems" and freely accessible at https://arxiv.org/abs/2403.03050



The folder "MATLAB code" contains the code needed to track the resistance changes over time of both networks with constant temperature and different currents. This allows us to observe an abrupt transition when the system reaches a steady state and to capture the plateau behavior during the transition. The main code is in the file "SC_dependent_networks_currents_simple.m" . In the absence of a Matlab license, one can use Octave as a compatible program to run the code.


The simulation process and function follow these steps:

Initialization: Characterize the system with its parameters and initialize the system in its initial state.

Loop: Modification of an external parameter in the system, dT or dI, in small steps. At each step, the following process is performed:

Step 1: Voltage_nodes (function) – Receives the system structure, resistances, and currents, computes the conductance matrix, and solves the Kirchhoff matrix equation. Returns the voltage at each node.

Step 2: System convergence check.

Step 3: find_currents  (function) – Receives the system structure, resistances, and voltages, and solves Kirchhoff’s equation for each resistor. Returns the current and heat for each resistor.

Step 4: heat_transfer_spectral  (function) – Receives the system structure and the heat dissipated by one network, and by using the diffusion equation, it computes the heat perceived by the other network.

Step 5: update_resistors_SC  (function) – Receives the system structure and its state at the current step, and determines the resistances and states of all resistors, using the Josephson I–V characteristic.

Repeat the steps until system convergence is achieved.

Measurement of the system’s total resistance.

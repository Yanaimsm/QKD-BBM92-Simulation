# QKD-BBM92-Simulation
This is our final Electrical Engineering Project.
It is a comprehensive MATLAB-based simulation framework for the BBM92 Quantum Key Distribution protocol. It encompasses the entire QKD process, including:

Parameter Initialization: Setting up simulation parameters such as detection efficiencies, dark count rates, and background noise levels.

Photon Pair Generation: Simulating entangled photon pairs and their polarization states.

Measurement Simulation: Modeling the measurement process at both Alice's and Bob's ends, including basis selection and detector responses.

Coincidence Detection: Determining valid detection events based on timing and basis alignment.

Key Sifting: Extracting raw keys from coincident detection events.

Error Correction: Implementing the Cascade algorithm to reconcile discrepancies between Alice's and Bob's keys.

Privacy Amplification: Reducing potential information leakage to eavesdroppers, ensuring the final key's security.

Additionally, the toolkit includes a parameter sweep module allowing users to study the impact of varying system parameters (e.g., background noise, detection efficiency) on the Quantum Bit Error Rate (QBER) and overall key generation performance.

This simulation environment serves as a valuable resource for researchers and students aiming to understand and analyze the dynamics of entanglement-based QKD systems.

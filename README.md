# Modeling the Effects of Climate on the Transmission Dynamics of High and Low Pathogenic Avian Influenza Using Computational Methods

## Introduction

A mathematical model(ODE system) is used to describe the transmission dynamics of avian influenza among a non-migratory waterfowl population. The manuscript under the same title was submitted to Mathematical Bioscience in Fall 2025 (currently still under review). This repository stores all the computational tools the authors used to produce the analytical and numerical results. 


## Avian Influenza Dynamics Solver 

To generate the transmission dynamics, the module named `AIVmodules.py` can be downloaded into your directory to solve for the intended system. The module includes functions: **Temp**, **Viral**, **Avian**, and **DynamicsSolver**. To generate the solutions and dynamics plot, users only need to call the **DynamicsSolver** function along with their initial conditions and required parameters for high pathogenic virus. It is recommended that users review the docstrings of each function carefully.








 

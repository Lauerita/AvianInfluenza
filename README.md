# Modeling the Effects of Climate on the Transmission Dynamics of High and Low Pathogenic Avian Influenza Using Computational Methods

## Introduction

A mathematical model(ODE system) is used to describe the transmission dynamics of avian influenza among a non-migratory waterfowl population. The manuscript under the same title was submitted to Mathematical Bioscience in Fall 2025 (currently still under review). This repository stores all the computational tools the authors used to produce the analytical and numerical results. 


## Avian Influenza Dynamics Solver 

To generate the transmission dynamics, the module named `AIVDynamics.py` can be downloaded into your directory to solve for the intended system. The module includes functions: **Temp**, **Viral**, **Avian**, and **DynamicsSolver**. To generate the solutions and dynamics plot, users only need to call the **DynamicsSolver** function along with their initial conditions and required parameters for high pathogenic virus. It is recommended that users review the docstrings of each function carefully.


## $R_0$ Calculator

The calculation for the time-invariant basic reproduction number is shown in the manuscript; one can simply compute $R_0$ by using the formula and desired parameters in the manuscript. However, a module for such calculation is provided in this repository under the file name `R0_calc.py`. Note that this file also plots the basic reproduction number as a function of temperature for $R_0^{LPAI}$, $R_0^{HPAI}$, as well as $R_0^{model}$. 





 

# pyranda
A Finite-Difference PDE solver framework.  

Pyranda is a simple Partial differential equation (PDE) solver framework for rapidly prototyping new algorithms and numerical methods.
Current capabilites are limited to a single structured mesh and equations are solved in strong conservative form.  The default integration 
scheme is a 5 stage RK4 method.  The finite difference method is used to approximate spatial derivatives and pyranda uses the FloATPy library 
for parallel, compact finite difference calculations.  

Dependencies include:
 - FloATPy
 - numpy, matplotlib, f2py

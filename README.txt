This code provides an algorithm for a steady state aerodynamic analysis of 3D wings using the Lifting Line Theory (LLT) with horseshore elements. It is based on the work done by Katz & Plotkin in Section 12.1, pg. 331 [1]
as well on the paper published by Phillips & Snyder [2]. 

The simulation parameters are defined in the steadyLLTrun.py file. In addition, to run the simulation, this file needs to be executed.

The output of the code depends on whether a range of angles of attack are considered or a single angle of attack is considered.

For a single angle of attack, the following will be outputed (see sampleResults/singleAOA):
	1. wing geometry plot.
	2. a plot showing the aerodyamic properties (such as aerodynamic coefficients) distribution along the wing span.
	3. a txt file with the values of the aerodynamic properties along the wing span.

For a range of angles of attack, the following will be outputed (see sampleResults/rangeOfAOA):
	1. a plot showing the polars of the aerodynamic coefficients (lift, drag, moment).
	2. a txt file with the values of the aerodynamic coefficients for each angle of attack.

Package requirements for running the simulations in a Conda virtual environment --> requirements.yml file.

Bibliography

[1] Katz, Joseph, and Allen Plotkin. Low-speed aerodynamics. Vol. 13. Cambridge university press, 2001.
[2] Phillips, Warren F., and D. O. Snyder. "Modern adaptation of Prandtl's classic lifting-line theory." Journal of Aircraft 37, no. 4 (2000): 662-670.
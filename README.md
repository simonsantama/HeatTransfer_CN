# HeatTransfer_CN
Crank-Nicholson solver for a 1-D heat transfer model.
Constant properties, homogenous and inert solid.

Not particularly efficient, but gets the job done.

1. Folder test_validation:

Heat diffusion in a semi-infinite solid with Dirichlet, Neunman and Robin boundary conditions validated against the analytical
solutions.

2. Folder implementation:
Uses the Crank-Nicolson scheme to predict the inert response of a one dimensional solid to different boundary conditions for a range of scenarios, including:
1) Different convective heat transfer coefficients
2) Different thermal conductivities
3) Different thermal diffusivities
4) Other


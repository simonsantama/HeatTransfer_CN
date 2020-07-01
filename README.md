# HeatTransfer_CN
Crank-Nicholson solver for a 1-D heat transfer model.
Constant properties, homogenous and inert solid.

0. test_validation/:

Heat diffusion in a semi-infinite solid with Dirichlet, Neunman and Robin boundary conditions validated against the analytical solutions.

1. BackSurface_Insulation/:

Uses CN scheme to predict the inert response of a one dimensional solid assuming an adiabatic back surface.
Several animations are created to visualize the evolution of the temperature profiles and absorbed energy for the different conditions studied.

a. Surface boundary conditions:

	Linear -> q_net = q_inc - hT(T_surf - T_gas). Where q_inc is a function of time and hT constant.

	Non-Linear -> q_net = q_inc - hc(T_surf - T_gas) - emissivity x stefan_boltzman x T_surf^4. Where q_inc is a function of time. hc and emissivity constant.



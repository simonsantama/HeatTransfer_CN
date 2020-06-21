# HeatTransfer_CN
Crank-Nicholson solver for a 1-D heat transfer model.
Constant properties, homogenous and inert solid.

Not particularly efficient, but gets the job done.

0. test_validation/:

Heat diffusion in a semi-infinite solid with Dirichlet, Neunman and Robin boundary conditions validated against the analytical solutions.

1. implementation_backinsulated/:

Uses CN scheme to predict the inert response of a one dimensional solid assuming an adiabatic back surface.

a. Surface boundary conditions:

	Linear -> q_net = q_inc - hT(T_surf - T_gas). q_inc is a function of time and hT constant.

	Non-Linear -> q_net = q_inc - hc(T_surf - T_gas) - emissivity x stefan_boltzman x T_surf^4. q_inc is a function of time. hc and emissivity constant.


2. implementation_backlosses/:

Uses CN scheme to predict the inert response of a one dimensional solid assuming conductive losses to an adiabatic aluminium block without a thermal contact resistance.

a. Surface boundary conditions:

	Linear -> q_net = q_inc - hT(T_surf - T_gas). q_inc is a function of time and hT constant.

	Non-Linear -> q_net = q_inc - hc(T_surf - T_gas) - emissivity x stefan_boltzman x T_surf^4. q_inc is a function of time. hc and emissivity constant.



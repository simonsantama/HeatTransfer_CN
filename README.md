# HeatTransfer_CN
Crank-Nicholson solver for a 1-D heat transfer model.
Constant properties, homogenous and inert solid.

Not particularly efficient, but gets the job done.

1. /test_validation/:

Heat diffusion in a semi-infinite solid with Dirichlet, Neunman and Robin boundary conditions validated against the analytical solutions.

2. /implementation_general/:

Uses CN scheme to predict the inert response of a one dimensional solid.
a. Surface boundary conditions:
	Linear -> q_net = q_inc - hT(T_surf - T_gas). q_inc is a function of time and hT constant.
	Non-Linear -> q_net = q_inc - hc(T_surf - T_gas) - emissivity x stefan_boltzman x T_surf^4. q_inc is a function of time. hc and emissivity constant.
b. Back boundary conditions: 
	Insulation -> q_conductive_back = 0
	Aluminium block -> q_conductive_back = q_net_alumium. Aluminium block at the back face. No thermal resistance. No heat losses. 20 mm depth.




! Predict along the energy parameter 
!   We observe that the curve of the control variables( x, z, vy, tp) w.r.t energy is almost linear 

! Take the second option!!!!
!   because the energy cj is not explicitly dependent on the period 
! 2 possible options
!   1:  use the derivative of energy w.r.t each control parameters   -- do not rely on the previous p.o. 
!   2:  use the difference bewteen the previous two p.o.s     -- rely on the previous p.o.s, but do not need to compute the derivative


! 	Input  Varaible 
!  curv 	the control variables in the initial conditions of the previous p.o.s, for general case, the full state vector 
!           	for generality, use as curv(5,*)
!  nds 		the number of current available rows in curv 
!  incrds 	flag to decide how to change the current arc step, 1: 2*ds, 0: remain, -1: ds/2


! 	Input  and Output Varaible 
!  yctr		the initial (and the updated) guess  of the prediction for the next p.o. 


 


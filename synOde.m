% Interactions between multiple sources of short term plasticity
% during evoked and spontaneous activity at the rat calyx of Held
% J Physiol, 2008
%
% Matthias H. Hennig, Michael Postlethwaite, Ian D. Forsythe, Bruce
% P. Graham
% MHH: mhhennig@gmail.com; BPG:  b.graham@cs.stir.ac.uk
%
% This function contains the set of ODEs of the synapse model:
% y(1) is the 1st inactivated state (i1)
% y(2) is the calcium transient amplitude (c1)
% y(3) represents effects of  mGluR/AMPAR activation (b)
% y(4) is the 2nd inactivated state (i2)
% y(5) is the vesicle pool occupancy (n)
% y(6) is the rate of activity-dependent vesicle recruitment (ke)

function dy=synOde(t,y)

global gcairel gcairel2 gfrel gglrel gkd gke kek ke

dy = zeros(6,1);

dy(1) = -y(1)/gcairel+y(4)/gcairel2;
dy(2) = -1/gfrel*(y(2) - 1+y(1)+y(3)+y(4) );
dy(3) = -y(3)/gglrel;
dy(4) = -y(4)/gcairel2;
dy(5) = (ke*y(6)+1/gkd)*(1-y(5));
dy(6) = -y(6)/kek;

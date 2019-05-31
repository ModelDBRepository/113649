% Interactions between multiple sources of short term plasticity
% during evoked and spontaneous activity at the rat calyx of Held
% J Physiol, 2008
%
% Matthias H. Hennig, Michael Postlethwaite, Ian D. Forsythe, Bruce
% P. Graham
% MHH: mhhennig@gmail.com; BPG:  b.graham@cs.stir.ac.uk
%
% This function is called to simulate a sequence of EPSCs in
% response to a series of presyanptic spikes.
%
% The funtion expects a series of inter-spike-intervals as argument
% (isi).
%
% It returns the following data:
% nresps - the normalised EPSC amplitude for each presyn. AP
% pprel - the release probability for each AP
% ns - the vesicle pool occypancy for each AP
% ppbase - the "basal" Ca++ calcium transient amplitude (c2) for each AP
% nrels - the number of released vesicles for each AP
% ppfac - the "release" Ca++ calcium transient amplitude (c1) for each AP
% rdess - level of postsyn. AMPAR desensitisation for each AP
% finalstate - the state of all variables after the last AP (allows
%              to continue simulations with different stimuli etc.)
% retrieved - amount of recruited vesicles for each AP

function [nresps, pprel, ns, ppbase, nrels, ppfac, ...
	  rdess, finalstate, retrieved] = releasef(isi)

% Below, all activation rates are given in units of 1/AP (per
% presyn. action potential), and decay time constants are given in
% units of seconds.

% variables required by synOde.m
global gcairel gcairel2 gfrel gglrel gkd gke kek ke

% globally control calcium channel modulation 
withca = 1;
% globally control slow calcium channel modulation
withslow = 1;

% Ca++ transient amplitude (release probability)
fac = 0.9894;
ca = 0.034*fac;

% facilitation activation
gfac = withca * 0.06;
% facilitation relaxation
gfrel = 0.04;

% activation of activity-dependent vesicle retrieval
kefrac = 1;
% decay time of activity-dependent vesicle retrieval
kek = 0.1;
ke = 6.0;
gke = 0.2373;

% passive vesicle recycling time constant
gkd   = 4.4;              

% fast Ca++ channel inactivation
gcai    = withslow * withca * 0.009;
% fast Ca++ channel relaxation
gcairel = 0.3;
% slow Ca++ channel inactivation
gcai2    = withslow * withca * 0.007;
% slow Ca++ channel relaxation
gcairel2 = 20;
% mGluR/AMPAR mediated Ca++ channel inhibition (presyn.)
ggli    = withslow * withca * 0.013;
% relaxation of mGluR/AMPAR effect
gglrel  = 10;    

% AMPAR desensitisation (postsyn.)
grr =  0.0230;
% AMPAR desensitisation recovery time
grd =  2.8955;

% factor to convert (ca*[0:1])^4 into release probability
k=138000*1.4;

% power law exponent for release probability
exponent = 4;

% ODE solver parameters
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare initial values

time = 1:length(isi)-1;

nrels = zeros(length(time),1);
rdess = zeros(length(time),1);
resps = zeros(length(time),1)';
ppbase = zeros(length(time),1);
ppfac = zeros(length(time),1);
pprel = zeros(length(time),1);

ns(1) = 1;
rdess(1) = 0;
pp = 1;
ppb = 1;
ppin = 0;
ppin2 = 0;
ppmglu  = 0;
prel(1)  = 1-exp(-k*(ca)^exponent);
pprel(1) = 1-exp(-k*(ca)^exponent);
ca_amp = pp;
ppbase(1) = ca_amp;
ppfac(1) = ca_amp;
enhrep = 0;
Y = [ 0 0 0 0 1 0];
T = [0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run simulation

for t=time,

  % next inter-spike interval
  dt = isi(t);

  % Desensitisation
  rr = 1-exp(-dt/grr);

  % Transmitter release
  nrels(t) = ns(t) * pprel(t);
  nt = ns(t) - nrels(t);
    
  % AMPAR desensitisation
  desstmp = 1-rdess(t);
  rdess(t+1) = rdess(t) + nrels(t) * grd * desstmp;
  rdess(t+1) = rdess(t+1) - rr * rdess(t+1);
   
  % postsynaptic response
  resps(t) = nrels(t) * desstmp;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % update all state variables in the model
  % mGluR activation
  pglut = ggli * ppb * nrels(t);
  % Ca++ channel inactivation  
  pcait = gcai *  ppb * ca_amp; 
  pcait2 = gcai2 * ppin * ca_amp;
  % Ca++ channel facilitation  
  facilt = gfac * ppb;
  % activity-dependent vesicle retrieval
  enhrep = enhrep + gke * ca_amp * (1-enhrep);
  
  % store number of vesicles
  retrieved(t) = nt;

  finalstate = [ ppin ppin2 ppmglu  desstmp (1-nt) enhrep];

  % solve the ODE system
  [T, Y] = ode45(@synOde, [0 dt], [ppin + pcait - pcait2, pp + facilt, ...
		    ppmglu + pglut, ppin2+pcait2, nt, enhrep],options);

  % find last element of calculated time course (used for backwards compatibility)
  tn = length(T);
  
  % release prob.
  pp = Y(tn,2);
  % inactivated Ca++ channels
  ppin = Y(tn,1);
  ppin2 = Y(tn,4);
  % mGluR effect on Ca++ channels
  ppmglu = Y(tn,3);
  % activity-dependent vesicle retrieval
  enhrep = kefrac*Y(tn,6);
  % vesicle pool size
  ns(t+1)  = Y(tn,5);
  
  % calculate number of retrieved vesicles
  retrieved(t) = ns(t+1)-retrieved(t);

  % base release probability
  ppb = 1-ppin-ppmglu-ppin2;
    
  % new release probability/Ca++ transient amplitude
  ca_amp = pp;
  ppbase(t+1) = ppb;
  ppfac(t+1) = ca_amp;
  pprel(t+1) = 1-exp(-k*(ca * ca_amp)^exponent);

end

% calculate normalised EPSCs
nresps = resps/resps(1);

% Interactions between multiple sources of short term plasticity
% during evoked and spontaneous activity at the rat calyx of Held
% J Physiol, 2008
%
% Matthias H. Hennig, Michael Postlethwaite, Ian D. Forsythe, Bruce
% P. Graham
% MHH: mhhennig@gmail.com; BPG:  b.graham@cs.stir.ac.uk
%
% This code produces a set of plots as shown in the paper in
% Figs. 2 and 3.

clear

% initialise graphics
h1 = figure(1);
clf
fs = [8.5 14]*1.5;
set(h1, 'PaperOrientation','portrait');
set(h1, 'PaperType','a4');
set(h1,'PaperUnits','centimeters');
set(h1,'Units','centimeters');
set(h1,'PaperPosition',[0 0 fs]);
set(h1,'Position',[0 19 fs]);
ms = 8;
colours = repmat([0 0 0],4,1)+repmat([0.4 0.8 0.6 0],3,1)';

%%%%%%%%%%%%%%%%%%%%%%%%%%

% load experimental data
load 'mean_epscs.mat'

frec = 0;
for freq = [10 20 50 100 ],

  frec = frec + 1;
  time = 1:freq;
  exptime = 1:freq;
  isi = repmat(1/freq, freq+1,1);
 

  % simulate...
  [resps, pprel, n, pb, nr, pf, rdes, final, retrieved ] = releasef(isi);

  % ...and plot the results

  subplot(3,2,1)
  p = plot(exptime./freq,epsc_exp{frec}(:,2),'x');
  set(p,'Color',colours(frec,:),'LineWidth',1,'MarkerSize',ms);
  hold on
  plot(time./freq,resps/resps(1),'k-','LineWidth',2,'Color',colours(frec,:))
  hold on

  subplot(3,2,2)
  loglog(exptime./freq,epsc_exp{frec}(:,2),'x','Color',colours(frec,:),'LineWidth',1,'MarkerSize',ms)
  hold on
  loglog(time./freq,resps/resps(1),'k-','LineWidth',2)
  
  subplot(3,2,3)
  plot(time./freq,pf(1:length(time)),'k.-','LineWidth',2,'Color',colours(frec,:))
  hold on
  
  subplot(3,2,4)
  plot(time./freq,pprel(1:length(time)),'k.-','LineWidth',2,'Color',colours(frec,:))
  hold on

  subplot(3,2,5)
  plot(time./freq,n(1:length(time)),'k.-','LineWidth',2,'Color',colours(frec,:))
  hold on

  subplot(3,2,6)
  plot(time./freq,retrieved,'k.-','LineWidth',1,'Color',colours(frec,:))
  hold on  

end

% set appropriate axis labels etc.

subplot(3,2,1)
set(gca,'XLimMode','manual');
set(gca,'XLim', [0.0 1.0]);
set(gca,'YLimMode','manual');
set(gca,'YLim', [0 1]);    
set(gca,'FontName','Helvetica-Narrow');
set(gca,'FontSize',8);
ylabel('Normalised EPSC')
xlabel('Time/s')

subplot(3,2,2)
set(gca,'XLimMode','manual');
set(gca,'XLim', [0.01 1]);
set(gca,'YLimMode','manual');
set(gca,'YLim', [0.1 1.0]);    
set(gca,'FontName','Helvetica-Narrow');
set(gca,'FontSize',8);
set(gca,'XTick',[0.01 0.1 1])
set(gca,'XTickLabel',[0.01 0.1 1])
set(gca,'YTick',[0.1  1])
set(gca,'YTickLabel',[ 0.1 1])
xlabel('Time/s')

subplot(3,2,3)
set(gca,'XLimMode','manual');
set(gca,'XLim', [0.0 1.0]);
set(gca,'YLimMode','manual');
set(gca,'YLim', [0.8 1.2]);    
set(gca,'FontName','Helvetica-Narrow');
set(gca,'FontSize',8);
ylabel('Calcium Transient')
xlabel('Time/s')

subplot(3,2,4)
set(gca,'XLimMode','manual');
set(gca,'XLim', [0.0 1.0]);
set(gca,'YLimMode','manual');
set(gca,'YLim', [0.1 0.4]);    
set(gca,'FontName','Helvetica-Narrow');
set(gca,'FontSize',8);
ylabel('Release Probability')
xlabel('Time/s')


subplot(3,2,5)
set(gca,'XLimMode','manual');
set(gca,'XLim', [0.0 1]);
set(gca,'YLimMode','manual');
set(gca,'YLim', [0 1.2]);    
set(gca,'FontName','Helvetica-Narrow');
set(gca,'FontSize',8);
ylabel('Release Pool Occupancy')
xlabel('Time/s')

subplot(3,2,6)
set(gca,'XLimMode','manual');
set(gca,'XLim', [0.0 1.0]);
%set(gca,'YLimMode','manual');
%set(gca,'YLim', [0.0 1.2]);    
set(gca,'FontName','Helvetica-Narrow');
set(gca,'FontSize',8);
xlabel('Time/s')
ylabel('Vesicles Retrieved')


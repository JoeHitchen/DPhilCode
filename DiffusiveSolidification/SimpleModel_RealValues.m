%% Salt Diffusion Investigation - RealValues
%
% This model investigates the effect of salt diffusion in the liquid phase
% on the growth rate for a variety of mushy diffusivities.
% The liquidus temperature at infinity can also be varied.
%
% It creates figure 3.5 in my pre-corrections thesis.
%
% (05/07/14)

clear
clc
addpath(genpath('../CoreFunctions'))


%% Properties

% Physical Parameters
LqdusInfVec = 0.01:0.01:0.99;
StefC_Ice = 30;
Lewis = 190;
StefC_Wettlaufer = 333.4e3/(4e3*0.085*70);
m1 = sqrt(1+StefC_Ice);
m2 = sqrt((1+StefC_Ice)/(1+StefC_Ice/Lewis));
mw = sqrt((1+StefC_Wettlaufer)/(1+StefC_Wettlaufer/Lewis));
s = sqrt(Lewis);


%% Solution Code

% Prepare
IntPos = NaN(4,length(LqdusInfVec));

% Loop over infinite liquidus and mushy diffusivity
tic
for iL = 1:length(LqdusInfVec)
    
    % Set parameters
    LqdusInf = LqdusInfVec(iL);
    
    % Find interface with infinite s
    IntPos(1,iL) = DiffGrow_Simp_Cond(m1, LqdusInf, inf, 0.1);
    IntPos(2,iL) = DiffGrow_Simp_Cond(m2, LqdusInf, inf, 0.1);
    IntPos(3,iL) = DiffGrow_Simp_Cond(m1, LqdusInf, sqrt(Lewis), 0.4);
    IntPos(4,iL) = DiffGrow_Simp_Cond(m2, LqdusInf, sqrt(Lewis), 0.5);
    IntPos(5,iL) = DiffGrow_Simp_Cond(mw, LqdusInf, sqrt(Lewis), 0.5);
end
toc

% Tidy
clear tmpM tmpLqdusInf FindInterfaceNoSalt FindInterfaceWiSalt FindInterfaceWettlaufer iL iM m LqdusInf

%% Output
FC = 0;
FS = 18;
LW = 3;

% Growth rates
if 1
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create graph (a)
    subplot(1,2,1)
    plot(LqdusInfVec,IntPos(1,:),'linewidth',LW)
    
    % Style (a)
    if 0
        title('(a) Ice Growth Rate Without Salt Diffusion','fontsize',FS)
        xlabel('Far field liquidus temperature, \theta_{L\infty}','fontsize',FS)
        ylabel('Growth Rate, \lambda','fontsize',FS)
        set(gca,'ylim',[0.00,0.8])
        set(gca,'xlim',[0.1,1])
        set(gca,'xtick',0.1:0.1:1,'ytick',0.0:0.2:0.8,'fontsize',FS)
        axis square
        box on
    else
        title('1TITLE','fontsize',FS)
        xlabel('1XLAB','fontsize',FS)
        ylabel('1YLAB','fontsize',FS)
        set(gca,'ylim',[0.00,0.8])
        set(gca,'xlim',[0.0,1])
        set(gca,'xtick',0.0:0.2:1,'ytick',0.0:0.1:0.8,'fontsize',FS)
        set(gca,'xticklabel',{'1X01','1X02','1X03','1X04','1X05','1X06','1X07','1X08','1X09','1X10'})
        set(gca,'yticklabel',{'1Y01','1Y02','1Y03','1Y04','1Y05','1Y06','1Y07','1Y08','1Y09'})
        axis square
        box on
    end
    
    
    
    % Create graph (b)
    subplot(1,2,2)
    tmpBase = [1;1;1]*IntPos(1,:);
    tmpAnalytical = 1/m1*sqrt(log(m1*LqdusInfVec./(1-LqdusInfVec)));
    hold on
    h = plot([1,1],[100,100],LqdusInfVec,(IntPos(2:4,:)-tmpBase)./tmpBase*100,'linewidth',LW);
    plot([0,1],[0,0],'k--')
    hold off
    
    % Style (b)
    if 0
        title('(b) Salt Diffusion Effects on Ice Growth Rate','fontsize',FS)
        xlabel('Far field liquidus temperature, \theta_{L\infty}','fontsize',FS)
        ylabel('Percentage Change in Growth Rate','fontsize',FS)
        legend(h(2:4),'Salt diffusion in mushy phase','Salt diffusion in liquid phase','Salt diffusion in both phases','fontsize',FS,'location','southeast')
        set(gca,'ylim',[-15,10])
        set(gca,'xlim',[0.1,1])
        set(gca,'xtick',0.1:0.1:1,'ytick',-15:5:10,'fontsize',FS)
        box on
        axis square
    else
        title('2TITLE','fontsize',FS)
        xlabel('2XLAB','fontsize',FS)
        ylabel('2YLAB','fontsize',FS)
        set(gca,'ylim',[-15,10])
        set(gca,'xlim',[0.0,1])
        set(gca,'xtick',0.0:0.2:1,'ytick',-15:5:10,'fontsize',FS)
        set(gca,'xticklabel',{'2X01','2X02','2X03','2X04','2X05','2X06','2X07','2X08','2X09','2X10'})
        set(gca,'yticklabel',{'2Y01','2Y02','2Y03','2Y04','2Y05','2Y06'})
        legend(h(2:4),'LEG1111111111111','LEG2','LEG3','fontsize',FS,'location','southeast')
        box on
        axis square
    end
    
    % Tidy
    clear tmpBase tmpAnalytical h
end

% Tidy
clear FC FS LW


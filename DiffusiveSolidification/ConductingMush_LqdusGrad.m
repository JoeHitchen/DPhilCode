%% Liquidus Gradient Investigation (Conducting, without salt)
% This script investigates the effects of the liquidus gradient, and to
% some extent the Stefan number, for a fixed initial liquidus temperature.
% The investigation is in the perfectly conducting limit and without salt
% diffusion.
%
% It produces figures 3.6, 3.7, 3.8 (a&b), and 3.9 of my pre-corrections
% thesis.
%
% (14/01/16)

clear
clc
addpath(genpath('../CoreFunctions'))


%% Setup Investigation

% Parameter space
LqdusInf     = 0.5;
StefanVec    = [0,10.^(-2:0.05:3)];
LqdusGradVec = [0,10.^(-3:0.2:1)];

% Computational Domain
Comp.DomainDepth = 15;
Comp.NGrid       = [189,89];
Comp.StefanIdx   = 54;


%% Run investigation

% Initialise result variables
GrowthRate_Simple = NaN(size(LqdusGradVec'*StefanVec));
GrowthRate_Mushy  = NaN(size(LqdusGradVec'*StefanVec));

Profs_Grids   = NaN(Comp.NGrid(1)+Comp.NGrid(2)+1,length(LqdusGradVec));
Profs_LqdFrac = NaN(Comp.NGrid(1)+Comp.NGrid(2)+1,length(LqdusGradVec));


% Loop over liquidus gr8adients
hTic = tic;
parfor iLg = 1:length(LqdusGradVec)
    
    % Initialise parallel variables
    PAR_GRM     = NaN(1,length(StefanVec));
    PAR_Grids   = NaN(Comp.NGrid(1)+Comp.NGrid(2)+1,1); %#ok<PFBNS>
    PAR_Temps   = NaN(Comp.NGrid(1)+Comp.NGrid(2)+1,1);
    PAR_LqdFrac = NaN(Comp.NGrid(1)+Comp.NGrid(2)+1,1);
    
    % Set initial guess
    IG  = [DiffGrow_Simp_Cond(1,LqdusInf),1/sqrt(pi),LqdusGradVec(iLg)/(LqdusGradVec(iLg)+LqdusInf)];
    
    % Loop over Stefan numbers
    for iSt = 1:length(StefanVec)
        
        % Finite liquidus gradient
        if LqdusGradVec(iLg) ~= 0
            
            % Get results
            [PAR_GRM(iSt),tmpProfs,tmpGrid,IG] = DiffGrow_Mushy_Cond_NS(StefanVec(iSt),LqdusGradVec(iLg),LqdusInf,Comp,IG);
            
        % Stefan limit
        else
            
            % Get results
            [PAR_GRM(iSt),tmpProfs,tmpGrid] = DiffGrow_Stefan(StefanVec(iSt),LqdusInf,Comp);
        end
        
        % Get profiles
        if iSt == Comp.StefanIdx
            PAR_Grids   = tmpGrid;
            PAR_Temps   = tmpProfs(:,2);
            PAR_LqdFrac = tmpProfs(:,5);
        end
    end
    
    % Parallel output
    GrowthRate_Mushy(iLg,:) = PAR_GRM;
    Profs_Grids(  :,iLg)    = PAR_Grids;
    Profs_Temps(  :,iLg)    = PAR_Temps;
    Profs_LqdFrac(:,iLg)    = PAR_LqdFrac;
    
    toc(hTic)
end

% Tidy
clear hTic iLg iSt PAR_GRS PAR_GRM PAR_Grids PAR_LqdFrac IG tmpProfs tmpTemps tmpGrids


%% Output
FC = 0;
LW = 3;
FS = 18;
FW = 'bold';
PSFrag = true;


% Self-similar growth rate (mushy) surface plot
if 1
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create graph
    CTicks = 0:0.1:1;
    Contours = 0.1:0.2:0.9;
    hold on
    contourf(LqdusGradVec(2:end),StefanVec(2:end),GrowthRate_Mushy(2:end,2:end)',50,'edgecolor','none')
    contour(LqdusGradVec(2:end),StefanVec(2:end),GrowthRate_Mushy(2:end,2:end)',Contours,'edgecolor','k','linewidth',LW)
    hold off
    if ~PSFrag
        ylabel('Stefan number,  S_t','fontsize',FS)
        set(gca,'yscale','log','ytick',10.^(-2:3))
        xlabel('Liquidus gradient,  \Lambda','fontsize',FS)
        set(gca,'xscale','log','xtick',10.^(-3:1))
        set(gca,'fontsize',FS)
        caxis(CTicks([1,end]))
        ylabel(...
            colorbar('ylim',CTicks([1,end]),'ytick',CTicks,'yticklabel',CTicks,'fontsize',FS),...
            'Mushy growth rate')
        box on
        axis square
        colormap(cbrewer('seq','YlGnBu',50))
    else
        ylabel('YLAB','fontsize',FS)
        set(gca,'yscale','log','ytick',10.^(-2:3),'yticklabel',{'Y1','Y2','Y3','Y4','Y5','Y6'})
        xlabel('XLAB','fontsize',FS)
        set(gca,'xscale','log','xtick',10.^(-3:1),'xticklabel',{'X1','X2','X3','X4','X5'})
        set(gca,'fontsize',FS)
        caxis(CTicks([1,end]))
        hC = colorbar('ylim',CTicks([1,end]),'ytick',CTicks,'yticklabel',{'C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11'},'fontsize',FS);
        ylabel(hC,'CLAB','fontsize',FS)
        for iC = Contours
            line('parent',hC,'xdata',[-5,5],'ydata',iC*[1,1],'linewidth',LW)
        end
        box on
        axis square
        colormap(cbrewer('seq','YlGnBu',50))
    end
    
    % Tidy
    clear CTicks
end


% Surface liquid fractions
if 1
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create graph
    LqdusGradVec2 = 10.^(-3:0.01:2);
    hold all
    plot(LqdusGradVec2,1./(1+0.5./LqdusGradVec2),'linewidth',LW);
    hold off
    if ~PSFrag
        xlabel('Liquidus gradient,  \Lambda','fontsize',FS)
        ylabel('Liquid fraction at surface,  \chi(0)','fontsize',FS)
        set(gca,'xscale','log','xlim',10.^[-3,2],'xtick',10.^(-3:2),'ytick',0:0.1:1,'fontsize',FS)
    else
        xlabel('XLAB','fontsize',FS)
        set(gca,'xscale','log','xlim',10.^[-3,2],'xtick',10.^(-3:2),'xticklabel',{'X1','X2','X3','X4','X5','X6'})
        ylabel('YLAB','fontsize',FS)
        set(gca,'ytick',0:0.1:1,'yticklabel',{'Y01','Y02','Y03','Y04','Y05','Y06','Y07','Y08','Y09','Y10','Y11'})
        set(gca,'fontsize',FS)
    end
    box on
    axis square
    
    % Tidy
    clear LqdusGradVec2
end


% Liquid fraction contour plot
if 1
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create graph
    CTicks = 0:0.1:1;
    hold on
    contourf(ones(Comp.NGrid(1)+1,1)*LqdusGradVec(2:end),Profs_Grids(1:Comp.NGrid(1)+1,2:end),Profs_LqdFrac(1:Comp.NGrid(1)+1,2:end),50,'edgecolor','none')
    contour( ones(Comp.NGrid(1)+1,1)*LqdusGradVec(2:end),Profs_Grids(1:Comp.NGrid(1)+1,2:end),Profs_LqdFrac(1:Comp.NGrid(1)+1,2:end),0.1:0.2:0.9,'edgecolor',0.5*[1,1,1],'linewidth',LW)
    contour( ones(Comp.NGrid(1)+1,1)*LqdusGradVec(2:end),Profs_Grids(1:Comp.NGrid(1)+1,2:end),Profs_LqdFrac(1:Comp.NGrid(1)+1,2:end),0.5,'edgecolor','r','linewidth',LW)
    plot(LqdusGradVec,GrowthRate_Mushy(:,Comp.StefanIdx),'k-','linewidth',LW)
    hold off
    if ~PSFrag
        ylabel('Depth,  \eta','fontsize',FS,'fontweight',FW)
        set(gca,'ydir','reverse')
        xlabel('Liquidus gradient,  \Lambda','fontsize',FS,'fontweight',FW)
        set(gca,'xscale','log','xtick',10.^(-3:1))
        set(gca,'fontsize',FS,'fontweight',FW)
        colormap(cbrewer('seq','Blues',50))
        caxis(CTicks([1,end]))
        hC = colorbar('ylim',CTicks([1,end]),'ydir','reverse','ytick',CTicks,'yticklabel',CTicks,'fontsize',FS,'fontweight',FW);
        ylabel(hC,'Liquid fraction,  \chi')
    else
        title('1TITLE','fontsize',FS)
        ylabel('1YLAB','fontsize',FS)
        set(gca,'ydir','reverse','yticklabel',{'1Y1','1Y2','1Y3','1Y4','1Y5'})
        xlabel('1XLAB','fontsize',FS)
        set(gca,'xscale','log','xtick',10.^(-3:1),'xticklabel',{'1X1','1X2','1X3','1X4','1X5'})
        set(gca,'fontsize',FS)
        colormap(cbrewer('seq','Blues',50))
        caxis(CTicks([1,end]))
        hC = colorbar('ylim',CTicks([1,end]),'ydir','reverse','ytick',CTicks,'yticklabel',{'1C01','1C02','1C03','1C04','1C05','1C06','1C07','1C08','1C09','1C10','1C11'},'fontsize',FS);
        ylabel(hC,'1CLAB')
    end
    box on
%    axis square
    line('parent',hC,'xdata',[-5,5],'ydata',[0.1,0.1],'color',0.5*[1,1,1],'linewidth',LW)
    line('parent',hC,'xdata',[-5,5],'ydata',[0.3,0.3],'color',0.5*[1,1,1],'linewidth',LW)
    line('parent',hC,'xdata',[-5,5],'ydata',[0.5,0.5],'color','r','linewidth',LW)
    line('parent',hC,'xdata',[-5,5],'ydata',[0.7,0.7],'color',0.5*[1,1,1],'linewidth',LW)
    line('parent',hC,'xdata',[-5,5],'ydata',[0.9,0.9],'color',0.5*[1,1,1],'linewidth',LW)
    
    % Tidy
    clear CTicks hC
end


% Scaled liquid fraction surface
if 1
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create data
    ScaledGrid = 1 - (...
        Profs_Grids(  1:Comp.NGrid(1)+1,:))...
        ./...
        (ones(Comp.NGrid(1)+1,1)*Profs_Grids(Comp.NGrid(1)+1,:)...
        );
    ScaledLqdFrac = ( 1-Profs_LqdFrac(1:Comp.NGrid(1)+1,:) )...
        ./...
        (ones(Comp.NGrid(1)+1,1)*( 1-Profs_LqdFrac(1,:) ));
    ScaledLqdFrac(ScaledLqdFrac>=1) = 1;
    
    % Create graph
    CTicks = 0:0.1:1;
    hold on
    contourf(ones(Comp.NGrid(1)+1,1)*LqdusGradVec(2:end),ScaledGrid(:,2:end),ScaledLqdFrac(:,2:end),50,'edgecolor','none')
    contour( ones(Comp.NGrid(1)+1,1)*LqdusGradVec,ScaledGrid,ScaledLqdFrac,0.1:0.2:0.9,'edgecolor',0.5*[1,1,1],'linewidth',LW)
    contour( ones(Comp.NGrid(1)+1,1)*LqdusGradVec,ScaledGrid,ScaledLqdFrac,0.5,'edgecolor',[1,0.5,0],'linewidth',LW)
    plot(LqdusGradVec,4.0*LqdusGradVec,'--','color',[1,0.5,0],'linewidth',LW)
    contour( ones(Comp.NGrid(1)+1,1)*LqdusGradVec,ScaledGrid,Profs_LqdFrac(1:Comp.NGrid(1)+1,:),0.5,'edgecolor','r','linewidth',LW)
    hold off
    cmap = flipud(cbrewer('seq','Blues',50));
    colormap(cmap)
    patch(10.^[-3,1,1,0,-1,-3],10.^[-4,-4,-1.5,-1.5,-1.7,-3.5],cmap(1,:),'edgecolor','none')
    if ~PSFrag
        xlabel('Liquidus gradient,  \Lambda','fontsize',FS,'fontweight',FW)
        set(gca,'xscale','log','xtick',10.^(-3:1))
        ylabel('Scaled distance from interface,  1-\eta/\lambda','fontsize',FS,'fontweight',FW)
        set(gca,'ydir','normal','yscale','log')
        set(gca,'fontsize',FS,'fontweight',FW)
        caxis(CTicks([1,end]))
        hC = colorbar('ylim',CTicks([1,end]),'ydir','reverse','ytick',CTicks,'yticklabel',1-CTicks,'fontsize',FS,'fontweight',FW);
        ylabel(hC,'Scaled liquid fraction,  1-(1-\chi)/(1-\chi_0)')
    else
        title('2TITLE','fontsize',FS)
        xlabel('2XLAB','fontsize',FS)
        set(gca,'xscale','log','xtick',10.^(-3:1),'xticklabel',{'2X1','2X2','2X3','2X4','2X5'})
        ylabel('2YLAB','fontsize',FS)
        set(gca,'ydir','normal','yscale','log','ylim',10.^[-3,0],'ytick',10.^(-3:0),'yticklabel',{'2Y1','2Y2','2Y3','2Y4','2Y5'})
        set(gca,'fontsize',FS)
        caxis(CTicks([1,end]))
        hC = colorbar('ylim',CTicks([1,end]),'ydir','normal','ytick',CTicks,'yticklabel',{'2C01','2C02','2C03','2C04','2C05','2C06','2C07','2C08','2C09','2C10','2C11'},'fontsize',FS);
        ylabel(hC,'2CLAB')
    end
    box on
%    axis square
    line('parent',hC,'xdata',[-5,5],'ydata',0.1*[1,1],'color',0.5*[1,1,1],'linewidth',LW)
    line('parent',hC,'xdata',[-5,5],'ydata',0.3*[1,1],'color',0.5*[1,1,1],'linewidth',LW)
    line('parent',hC,'xdata',[-5,5],'ydata',0.5*[1,1],'color',[1,0.5,0],'linewidth',LW)
    line('parent',hC,'xdata',[-5,5],'ydata',0.7*[1,1],'color',0.5*[1,1,1],'linewidth',LW)
    line('parent',hC,'xdata',[-5,5],'ydata',0.9*[1,1],'color',0.5*[1,1,1],'linewidth',LW)
    
    % Tidy
    clear ScaledGrid ScaledLqdFrac CTicks hC
end


% Difference from Stefan Limit
if 1
   
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create graph
    CTicks = [0.001,0.003,0.01,0.03,0.1,0.3,1,3,10];
    Data = RelErr(GrowthRate_Mushy(2:end,2:end)',(ones(length(LqdusGradVec)-1,1)*GrowthRate_Mushy(1,2:end))').*(ones(length(StefanVec)-1,1)*LqdusGradVec(2:end).^(0));
    hold on
    contourf(LqdusGradVec(2:end),StefanVec(2:end),log10(Data),50,'edgecolor','none')
    contour( LqdusGradVec(2:end),StefanVec(2:end),log10(Data),log10(CTicks),'edgecolor','k','linewidth',1)
    contour( LqdusGradVec(2:end),StefanVec(2:end),log10(Data),[-1,10],'edgecolor','r','linewidth',LW)
    hold off
    if ~PSFrag
        ylabel('Stefan number,  S_t','fontsize',FS)
        set(gca,'yscale','log','ytick',10.^(-2:3))
        xlabel('Liquidus gradient,  \Lambda','fontsize',FS)
        set(gca,'xscale','log','xtick',10.^(-3:1))
        set(gca,'fontsize',FS)
        colorbar
        caxis(log10(CTicks([1,end])))
        hC = colorbar('ylim',log10(CTicks([1,end])),'ytick',log10(CTicks),'yticklabel',CTicks,'fontsize',FS);
        ylabel(hC,'Fractional growth rate increase,  (\lambda-\lambda_S)/\lambda_S')
        box on
        axis square
        colormap(cbrewer('seq','YlGnBu',50))
        line('parent',hC,'xdata',[-5,5],'ydata',[-1,-1],'color','k','linewidth',LW)
    else
        ylabel('YLAB','fontsize',FS)
        set(gca,'yscale','log','ytick',10.^(-2:3),'yticklabel',{'Y1','Y2','Y3','Y4','Y5','Y6'})
        xlabel('XLAB','fontsize',FS)
        set(gca,'xscale','log','xtick',10.^(-3:1),'xticklabel',{'X1','X2','X3','X4','X5'})
        set(gca,'fontsize',FS)
        colorbar
        caxis(log10(CTicks([1,end])))
        hC = colorbar('ylim',log10(CTicks([1,end])),'ytick',log10(CTicks),'yticklabel',{'C1','C2','C3','C4','C5','C6','C7','C8','C9'},'fontsize',FS);
        ylabel(hC,'CLAB')
        box on
        axis square
        colormap(cbrewer('seq','YlGnBu',50))
        for iC = CTicks
            line('parent',hC,'xdata',[-5,5],'ydata',log10(iC)*[1,1],'color','k','linewidth',1)
        end
        line('parent',hC,'xdata',[-5,5],'ydata',[-1,-1],'color','r','linewidth',LW)
    end
    
    % Tidy
    clear CTicks hC Data
end   

%%% Plots that possibly belong in other investigations


% Relative difference surface
if 0
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create graph
    CTicks = [0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1];
    hold on
    contourf(LqdusGradVec,StefanVec(2:end),log10(RelErr(GrowthRate_Simple(:,2:end),GrowthRate_Mushy(:,2:end),1))',50,'edgecolor','none')
    contour( LqdusGradVec,StefanVec(2:end),log10(RelErr(GrowthRate_Simple(:,2:end),GrowthRate_Mushy(:,2:end),1))',[-2,-1],'edgecolor','k','linewidth',LW)
    hold off
    ylabel('Stefan number,  S_t','fontsize',FS)
    set(gca,'yscale','log','ytick',10.^(-2:3))
    xlabel('Liquidus gradient,  \Lambda','fontsize',FS)
    set(gca,'xscale','log','xtick',10.^(-3:1),'fontsize',FS)
    colormap(cbrewer('seq','Reds',50))
    caxis(log10(CTicks([1,end])))
    hC = colorbar('ylim',log10(CTicks([1,end])),'ytick',log10(CTicks),'yticklabel',CTicks,'fontsize',FS);
    ylabel(hC,'Fractional error between growth rates')
    line('parent',hC,'xdata',[-5,5],'ydata',[-1,-1],'color','k','linewidth',LW)
    line('parent',hC,'xdata',[-5,5],'ydata',[-2,-2],'color','k','linewidth',LW)
    box on
    
    % Tidy
    clear CTicks hC
end
















%



% Scaled liquid fraction surface
if 0
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create data
    ScaledGrid = (1 - (...
        Profs_Grids(  1:Comp.NGrid(1)+1,:))...
        ./...
        (ones(Comp.NGrid(1)+1,1)*Profs_Grids(Comp.NGrid(1)+1,:)...
        )) .* (ones(Comp.NGrid(1)+1,1)*LqdusGradVec.^-1);
    %ScaledLqdFrac = ( 1-Profs_LqdFrac(1:Comp.NGrid(1)+1,:) )...
    %    ./...
    %    (ones(Comp.NGrid(1)+1,1)*( 1-Profs_LqdFrac(1,:) ));
    ScaledLqdFrac = 1-Profs_LqdFrac(1:Comp.NGrid(1)+1,:);
    
    % Create graph
    CTicks = 0:0.1:1;
    hold on
    contourf(ones(Comp.NGrid(1)+1,1)*LqdusGradVec(2:end),ScaledGrid(:,2:end),ScaledLqdFrac(:,2:end),50,'edgecolor','none')
    contour( ones(Comp.NGrid(1)+1,1)*LqdusGradVec(2:end),ScaledGrid(:,2:end),ScaledLqdFrac(:,2:end),0.5,'edgecolor',[1,0.5,0],'linewidth',LW)
    contour( ones(Comp.NGrid(1)+1,1)*LqdusGradVec(2:end),ScaledGrid(:,2:end),ScaledLqdFrac(:,2:end),20,'edgecolor','k')
    plot(10.^[-3,-1],ScaledGrid(:,5)*[1,1])
    plot(LqdusGradVec,LqdusGradVec.^-1,'-','color',[0,0,0],'linewidth',LW)
    contour( ones(Comp.NGrid(1)+1,1)*LqdusGradVec,ScaledGrid,Profs_LqdFrac(1:Comp.NGrid(1)+1,:),0.5,'edgecolor','r','linewidth',LW)
    patch(10.^[-3,1,1],10.^[3,-1,3],[0,0,0])
    hold off
    xlabel('Liquidus gradient,  \Lambda','fontsize',FS,'fontweight',FW)
    set(gca,'xscale','log','xlim',10.^[-3,-1],'xtick',10.^(-3:1))
    ylabel('Scaled distance from interface,  \Lambda(1-\eta/\lambda)','fontsize',FS,'fontweight',FW)
    set(gca,'ydir','normal','yscale','log','ylim',10.^[-1,2])
    set(gca,'fontsize',FS,'fontweight',FW)
    caxis(CTicks([1,end]))
    hC = colorbar('ylim',CTicks([1,end]),'ydir','reverse','ytick',CTicks,'yticklabel',1-CTicks,'fontsize',FS,'fontweight',FW);
    ylabel(hC,'Solid fraction,  1-\chi')
    box on
    colormap(cbrewer('seq','Blues',50))
    line('parent',hC,'xdata',[-5,5],'ydata',[0.5,0.5],'color',[1,0.5,0],'linewidth',LW)
    
    % Tidy
    %clear ScaledGrid ScaledLqdFrac CTicks hC
end


% Scaled liquid fraction surface
if 0
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create data
    ScaledTemp = Profs_Temps(1:Comp.NGrid(1)+1,:)/LqdusInf;
    ScaledGrid = ( Profs_Grids(1:Comp.NGrid(1)+1,:) ./ ...
        (ones(Comp.NGrid(1)+1,1)*Profs_Grids(Comp.NGrid(1)+1,:)) ) .* ...
        (ones(Comp.NGrid(1)+1,1)*LqdusGradVec.^-1);
    
    % Create graph
    hold on
    contourf(ones(Comp.NGrid(1)+1,1)*LqdusGradVec(2:end),ScaledGrid(:,2:end),log10(ScaledTemp(:,2:end)),-4:0.1:0)
    patch(10.^[-3,1,1],10.^[3,-1,3],[0,0,0])
    hold off
    xlabel('Liquidus gradient,  \Lambda','fontsize',FS,'fontweight',FW)
    set(gca,'xscale','log','xlim',10.^[-3,-1],'xtick',10.^(-3:1))
    ylabel('Scaled distance from interface,  \Lambda(1-\eta/\lambda)','fontsize',FS,'fontweight',FW)
    set(gca,'yscale','log','ylim',10.^[-1,+2])
    set(gca,'fontsize',FS,'fontweight',FW)
    caxis(log10([0.001,1]))
    colorbar('ytick',log10(10.^(-3:0.5:0)),'yticklabels',[0.001,0.003,0.01,0.03,0.1,0.3,1],'fontsize',FS,'fontweight',FW)
    colormap(cbrewer('seq','Blues',50))
end






%%% Finish up

% Tidy
clear FC LW FS



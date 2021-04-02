% CA Current Analysis File
% This file calls in results from prior runs and graphs them
clear all
close all

% SensCase= 1 % Base Case


SensCase=1;

STR='C:\Users\jsanchir\Box'; 

PATH=[STR '\Matlab 2018\CACurrentCode_GitHub\'];

% Used in base case 
PelicanAVG=[95.2287;75.3698;88.332;80.415];

SAVEFIG=1;

load (strcat(PATH,'Solution_Tim')) 
load (strcat(PATH,'openaccess_dynamics'))
load (strcat(PATH,'openaccess_pelicanT'))
load (strcat(PATH,'optimal_pelicanT'))
load (strcat(PATH,'SolutionSavesAll'))
load (strcat(PATH,'Parameters'),'data_params')
load (strcat(PATH,'Solution_NPV'),'NPV')

PELICAN=1; i=1; OA_CASE=7;
% Unpack the data 
Param = data_params(12:13);
 flag= data_params(1);
 capT=data_params(2); % indicator variable that was for whether imposing the T conditions
 cc_area = data_params(3);
 INITCOND = data_params(5:11); 
  BETA=data_params(14:16); 
  dis=data_params(4);

  
for II=1:4
    for JJ=1:5
        for KK=1:7
        Timenop{II,JJ,KK}=DATA{II,JJ,KK};
        Xanop{II,JJ,KK}=DATA{II,JJ+5,KK};
        Nanop{II,JJ,KK}=DATA{II,JJ+10,KK};
        Xsnop{II,JJ,KK}=DATA{II,JJ+15,KK};
        Nsnop{II,JJ,KK}=DATA{II,JJ+20,KK};
        Xhnop{II,JJ,KK}=DATA{II,JJ+25,KK};
        Nhnop{II,JJ,KK}=DATA{II,JJ+30,KK};
        Npnop{II,JJ,KK}=DATA{II,JJ+35,KK};
        E1nop{II,JJ,KK}=DATA{II,JJ+40,KK};
        E2nop{II,JJ,KK}=DATA{II,JJ+45,KK};
        E3nop{II,JJ,KK}=DATA{II,JJ+50,KK};
        end
    end
end
  
Time=Timenop{1,1,1};
time=Time(Time<35);
L=length(time); trun=1:L;  

for VarFlag=1:4
     
    
     [xazero, xszero, a,~,Mpj,tjuv, theta, vbk,kappa,d,wr,Y, alphaa, alphas, Cmax, Mp, amp_a, amp_s, per_a, per_s, start_a, start_s,...
            recbara, wra, Ma, winfa, kappaa, recbars, wrs, Ms, winfs, kappas,...
            z, phi, Kbp, theta_1, theta_2, theta_3, Gamma_a, Gamma_s, Gamma_other, Mbp, ...
            theta_1s, theta_2s, theta_3s, P,q,c0,c1]=CC_Param_File(flag, VarFlag-1, capT,Param,OA_CASE);

    
    for CASE=1:7
     
     E1=E1nop{VarFlag,PELICAN,CASE,i}(trun);
     E2=E2nop{VarFlag,PELICAN,CASE,i}(trun);
     E3=E3nop{VarFlag,PELICAN,CASE,i}(trun);
     Xa=Xanop{VarFlag,PELICAN,CASE,i}(trun);
     Xs=Xsnop{VarFlag,PELICAN,CASE,i}(trun);
     Xh=Xhnop{VarFlag,PELICAN,CASE,i}(trun);   
        Ra=P(1)*q(1)*E1.*(Xa.^BETA(1))-c0(1)*(E1)-c1(1)*(E1.^2); 
        Rs=P(2)*q(2)*E2.*(Xs.^BETA(2))-c0(2)*(E2)-c1(2)*(E2.^2);
        Rh=P(3)*q(3)*E3.*(Xh.^BETA(3))-c0(3)*(E3)-c1(3)*(E3.^2);
        
        f = fittype('pchipinterp');
                     g1a=fit(time,exp(-dis.*time).*Ra,f);
                     int1a=integrate(g1a,0:1:time(end),0);
                     g2s=fit(time,exp(-dis.*time).*Rs,f);
                     int2s=integrate(g2s,0:1:time(end),0);
                     g3h=fit(time,exp(-dis.*time).*Rh,f);
                     int3h=integrate(g3h,0:1:time(end),0);
                     
                     NpvA(VarFlag,CASE)=int1a(end);
                     NpvS(VarFlag,CASE)=int2s(end);
                     NpvH(VarFlag,CASE)=int3h(end);
        
                     
                 
    end
end

% Calculating per-period rent (at collocation nodes) for the different
% cases
for Varflag=1:4
    for PELICAN=1:5
        for CASE=1:7
            E1= E1nop{Varflag,PELICAN,CASE,i}(trun);
            E2=E2nop{Varflag,PELICAN,CASE,i}(trun);
            E3=E3nop{Varflag,PELICAN,CASE,i}(trun);
            Xa=Xanop{Varflag,PELICAN,CASE,i}(trun);
            Xs=Xsnop{Varflag,PELICAN,CASE,i}(trun);
            Xh=Xhnop{Varflag,PELICAN,CASE,i}(trun);            
            Rent_a{Varflag, PELICAN,CASE,i}=P(1)*q(1)*E1.*(Xa.^BETA(1))-c0(1)*(E1)-c1(1)*(E1.^2);
            Rent_s{Varflag, PELICAN, CASE,i}=P(2)*q(2)*E2.*(Xs.^BETA(2))-c0(2)*(E2)-c1(2)*(E2.^2);
            Rent_h{Varflag, PELICAN, CASE,i}=P(3)*q(3)*E3.*(Xh.^BETA(3))-c0(3)*(E3)-c1(3)*(E3.^2);
        end
    end
end


%% Figures in the main paper

T5550=0;
close all
C1=[0.0000	0.1262	0.3015]; % A alone[1 0 0]; % Red
C2=[0.2728	0.3144	0.4206]; %S Alone[0 0 0]; % Black
C3=rgb('Olive'); %[0.5970	0.5709	0.4696]; %H Alone%[0 1 0]; % Green
C4=rgb('SaddleBrown');[0.7560	0.6979	0.4290]; %A+S
C5=[0.941176474094391 0.941176474094391 0.941176474094391]; %[0.211	0.211	0.211];
C6=[0.7560	0.6979	0.4290]; %A+S alone
C7=rgb('Wheat'); %S+H
C8=rgb('Tomato'); %EBFM


time=Time(Time<35);
L=length(time); trun=1:L;

ddetime(1)=0;
for tt=2:length(OA_DYNAMICS{1,8}(1,:)')
    ddetime(tt)=.1+ddetime(tt-1);
end


CC=[C1; C2; C3;C4;C6;C7; C8];
CC1=[C8; C7; C6;C4;C3;C2; C1];


figure1 = figure;

% Create subplot
for Varflag=1:4
    W(Varflag)=OA_pelicanT{Varflag,7};
end

for Varflag=1:3
    subplot1 = subplot(1,3,Varflag,'Parent',figure1);
    w1 = 0.825;   w2 = .65;
    for i=1:7
        X(i)=100*(OA_pelicanT{Varflag,i} - W(Varflag))./W(Varflag);
        Y(i)=100*(Opt_pelicanT{Varflag,i}- W(Varflag))./W(Varflag); %-(OA_pelicanT{Varflag,i}-Opt_pelicanT{Varflag,i});
    end
    
    
    
    X=flip(X); Y=flip(Y);
    if Varflag<3
        
        for i=1:7
            bar1(i)=bar(i,X(i),w1,'FaceColor',CC1(i,:),'EdgeColor',CC1(i,:),'Parent',subplot1)
            hold on
        end
        hold off
        
        hold on
        bar(Y,w2,'FaceColor',C5,'EdgeColor',C1,'Parent',subplot1)
        if Varflag==1
            plot(0:8,100*ones(1,9)*(95.2287- W(Varflag))./W(Varflag),'k','LineWidth',2,'Parent',subplot1)
        elseif Varflag==2
            plot(0:8,100*ones(1,9)*(75.3698 - W(Varflag))./W(Varflag),'k','LineWidth',2,'Parent',subplot1)
        elseif Varflag==3
            plot(0:8,100*ones(1,9)*(88.332- W(Varflag))./W(Varflag),'k','LineWidth',2,'Parent',subplot1)
        else
            plot(0:8,100*ones(1,9)*(80.415- W(Varflag))./W(Varflag),'k','LineWidth',2,'Parent',subplot1)
        end
        
        hold off
        grid on
        
        
        set(subplot1,'XTick',[1 2 3 4 5 6 7],'XTickLabel',...
            {' ',' ',' ',' ',' ',' ',' '});
        
    end
    
    
    for i=1:7
        bar1=bar(i,X(i),w1,'FaceColor',CC1(i,:),'EdgeColor',CC1(i,:),'Parent',subplot1)
        hold on
    end
    hold off
    
    hold on
    bar(Y,w2,'FaceColor',C5,'EdgeColor',C1,'Parent',subplot1)
    if Varflag==1
        plot(0:8,100*ones(1,9)*(95.2287- W(Varflag))./W(Varflag),'k','LineWidth',2,'Parent',subplot1)
    elseif Varflag==2
        plot(0:8,100*ones(1,9)*(75.3698 - W(Varflag))./W(Varflag),'k','LineWidth',2,'Parent',subplot1)
    elseif Varflag==3
        plot(0:8,100*ones(1,9)*(88.332- W(Varflag))./W(Varflag),'k','LineWidth',2,'Parent',subplot1)
    else
        plot(0:8,100*ones(1,9)*(80.415- W(Varflag))./W(Varflag),'k','LineWidth',2,'Parent',subplot1)
    end
    hold off
    grid on
    
    
    box(subplot1,'on');
    % Set the remaining axes properties
    ax = subplot1;
    ax.XTick = [1 2 3 4 5 6 7];
    ax.XTickLabels = {'A+S+H','S+H','A+H','S+A','H','S','A'};
   
    ax.XTickLabelRotation = 45;
    ylim(subplot1,[0 140]);
  
    if Varflag==1
        ylabel({'Percent Increase in Pelican at T from Open Access '})
    end
    if Varflag==1
        title('No variability')
    elseif  Varflag==2
        title('Synchrony')
    elseif Varflag==3
        title({'Anti-synchrony'})
    else
        title({'Different period','same starting phase'})
    end
    
end


LineTYPE={'-b','-r','-k','-.b','-.r','-.k',':g'};

PELICAN=1; i=1;

figure1=figure;
clear Z
clear plot1
for j=1:9
    subplot1 = subplot(3,3,j,'Parent',figure1);
    
    if j<4
        VARFLAG=1;
    elseif (3<j) && (j<7)
        VARFLAG=2;
    elseif (6<j) && (j<9)
        VARFLAG=3;
    end
   
    if j==1||j==4||j==7||j==13
        for CASE=1:7
            Z(:,CASE)=Xanop{VARFLAG,PELICAN,CASE,i}(trun)./INITCOND(1);
        end
        if j==1
            title('Anchovy', 'Parent',subplot1)
            ylabel('Biomass', 'Parent',subplot1); %set(subplot1,'ylabel','Fishing Mortality Rate');
        end
        set(subplot1,'XTickLabel',{});
        
    elseif j==2||j==5||j==8||j==14
        for CASE=1:7
            Z(:,CASE)=Xsnop{VARFLAG,PELICAN,CASE,i}(trun)./INITCOND(3);
        end
        if j==2
            title( 'Sardine ', 'Parent',subplot1)
        end
        set(subplot1,'XTickLabel',{});
    elseif j==3||j==6||j==9||j==15
        for CASE=1:7
            Z(:,CASE)=Xhnop{VARFLAG,PELICAN,CASE,i}(trun)./INITCOND(5);
        end
        if j==3
            title( 'Halibut ', 'Parent',subplot1)
        end
        set(subplot1,'XTickLabel',{});
    
    end
    
    
   
        
    
    plot1=plot(time,Z(:,1),time,Z(:,2),time,Z(:,3),time,Z(:,4),time,Z(:,5),...
        time,Z(:,6),time,Z(:,7),'Parent',subplot1);
    
    set(plot1(1),'DisplayName','A','Color',CC(1,:),'LineWidth',2,'LineStyle','--'); %,'MarkerSize',5,'Marker','o');
    set(plot1(2),'DisplayName','S','Color',CC(2,:),'LineWidth',2,'LineStyle',':'); %,'MarkerSize',5,'Marker','+');
    set(plot1(3),'DisplayName','H','Color',CC(3,:),'LineWidth',2,'LineStyle','-.'); %,'MarkerSize',5,'Marker','s');
   
    set(plot1(4),'DisplayName','A+S','Color',CC(4,:),'LineWidth',2,'LineStyle','-');
    set(plot1(5),'DisplayName','A+H','Color',CC(5,:),'LineWidth',2,'LineStyle','-');
    set(plot1(6),'DisplayName','S+H','Color',CC(6,:),'LineWidth',2,'LineStyle','-');
    set(plot1(7),'DisplayName','A+S+H','Color',CC(7,:),'LineWidth',2,'LineStyle','-');
    
    axis(subplot1,[0 round(time(end)) 0 2])
    
    
     if j==1||j==4||j==7
            hold on;
            OAD=interp1(ddetime,OA_DYNAMICS{VARFLAG,8}(1,:)', time);
            plot2=plot(time,OAD./INITCOND(1),'b--','Parent',subplot1);
            set(plot2,'DisplayName','','Color',rgb('Black'),'LineWidth',2,'LineStyle','--')
            hold off;
            ylim(subplot1,[0 3]);
        end
        if j==2||j==5||j==8
            hold on;
            OAD=interp1(ddetime,OA_DYNAMICS{VARFLAG,8}(3,:)', time);
            plot2=plot(time,OAD./INITCOND(3),'b--','Parent',subplot1);
            set(plot2,'DisplayName','','Color',rgb('Black'),'LineWidth',2,'LineStyle','--')
            hold off;
            ylim(subplot1,[0 3]);
        end
        if j==3||j==6||j==9
            hold on;
            OAD=interp1(ddetime,OA_DYNAMICS{VARFLAG,8}(5,:)', time);
            plot2=plot(time,OAD./INITCOND(5),'b--','Parent',subplot1);
            set(plot2,'DisplayName','','Color',rgb('Black'),'LineWidth',2,'LineStyle','--')
            hold off;
            ylim(subplot1,[0 4]);
        end
    
    if j==9
        legend1 = legend([plot1(1) plot1(2) plot1(3) plot1(4) plot1(5) plot1(6) plot1(7)],...
            {'A','S', 'H', 'A+S','A+H','S+H','A+S+H'});
        set(legend1,...
            'Position',[0.188649534480679 0.354341630261094 0.648979586241196 0.0286713280044235],...
            'Orientation','horizontal');
    end
    
    box(subplot1,'on');
    grid(subplot1,'on');
    
    if j==1
        ylabel({'No Variability','Biomass'}, 'Parent',subplot1);
    elseif  j==4
        ylabel({'Synchrony','Biomass'}, 'Parent',subplot1);
    elseif j==7
        ylabel({'Anti-synchrony','Biomass'}, 'Parent',subplot1);
    elseif j==13
        ylabel({'Different Period','Biomass'}, 'Parent',subplot1);
    end
    
    if j<7
        set(subplot1,'XTickLabel',{});
    end
    
   
    
    if j>6
       
        xlabel('Time','Parent',subplot1);
    end
    
    if j==1
        title('Anchovy')
    elseif  j==2
        title('Sardine')
    elseif j==3
        title('Halibut')
    
    end
    
    ylim(subplot1,[0 2.5]);
    
    
end


%% Figure in main paper with HCR and Hockey Stick
PELICAN=1; clear Z
figure1=figure
ddetime(1)=0;
for tt=2:length(OA_DYNAMICS{VARFLAG,8}(1,:)')
    ddetime(tt)=.1+ddetime(tt-1);
end
for j=1:9
    subplot1 = subplot(3,3,j,'Parent',figure1);
    
    if j<4
        VARFLAG=1;
    elseif (3<j) && (j<7)
        VARFLAG=2;
    elseif (7<j) && (j<10)
        VARFLAG=3;
    elseif j>12
        VARFLAG=4;
    end
    
    if j==1||j==4||j==7||j==10
        for CASE=[1 2 4 6 7]
            Z(:,CASE)=q(2)*E2nop{VARFLAG,PELICAN,CASE,i}(trun);
        end
        Z(:,1)=interp1(ddetime,q(2)*OA_DYNAMICS{VARFLAG,7}(9,:)', time);
        Z(:,3)=interp1(ddetime,q(2)*OA_DYNAMICS{VARFLAG,9}(9,:)', time);
            Z(:,5)=interp1(ddetime,q(2)*OA_DYNAMICS{VARFLAG,10}(9,:)', time);
        set(subplot1,'XTickLabel',{});
        
    elseif j==2||j==5||j==8||j==11
        for CASE=[1 2 4 6 7]
            Z(:,CASE)=Xsnop{VARFLAG,PELICAN,CASE,i}(trun)./INITCOND(3);
        end
            Z(:,1)=interp1(ddetime,OA_DYNAMICS{VARFLAG,7}(3,:)', time)./INITCOND(3);
            Z(:,3)=interp1(ddetime,OA_DYNAMICS{VARFLAG,9}(3,:)', time)./INITCOND(3);
            Z(:,5)=interp1(ddetime,OA_DYNAMICS{VARFLAG,10}(3,:)', time)./INITCOND(3);
        if j==2
            title( 'Sardine ', 'Parent',subplot1)
        end
        set(subplot1,'XTickLabel',{});
    else j==3||j==6||j==9||j==12
        for CASE=[1 2 4 6 7]
            Z(:,CASE)=Npnop{VARFLAG,PELICAN,CASE,i}(trun);
        end
            Z(:,1)=interp1(ddetime,OA_DYNAMICS{VARFLAG,7}(7,:)', time);
            Z(:,3)=interp1(ddetime,OA_DYNAMICS{VARFLAG,9}(7,:)', time);
            Z(:,5)=interp1(ddetime,OA_DYNAMICS{VARFLAG,10}(7,:)', time);
            
        set(subplot1,'XTickLabel',{});
        
    end
    
    plot1=plot(time,Z(:,1),time,Z(:,2),time,Z(:,3),time,Z(:,5),...
        time,Z(:,7),'Parent',subplot1);
    
    set(plot1(1),'DisplayName','A','Color',CC(1,:),'LineWidth',2,'LineStyle','--'); %,'MarkerSize',5,'Marker','o');
    set(plot1(2),'DisplayName','S','Color',CC(2,:),'LineWidth',2,'LineStyle',':'); %,'MarkerSize',5,'Marker','+');
    set(plot1(3),'DisplayName','S_H_S','Color',rgb('Blue'),'LineWidth',2,'LineStyle','-.'); %,'MarkerSize',5,'Marker','s');
   
    set(plot1(4),'DisplayName','S_H_C_R','Color',rgb('Green'),'LineWidth',2,'LineStyle','-');
    
    set(plot1(5),'DisplayName','A+S+H','Color',CC(7,:),'LineWidth',2,'LineStyle','-');
    
    axis(subplot1,[0 round(time(end)) 0 max(max(Z))])
    if j==3||j==6||j==9
        
       axis(subplot1,[0 round(time(end)) 50 100])
    elseif j==1||j==4||j==7
        
        axis(subplot1,[0 round(time(end)) 0 1])
   elseif j==2||j==5||j==8
        
     
    end
    
    if j==9
        legend1 = legend([plot1(1) plot1(2) plot1(3)  plot1(4) plot1(5)],...
            {'OA','S', 'S_H_S','S_H_C_R','A+S+H'});
        set(legend1,...
            'Position',[0.188649534480679 0.354341630261094 0.648979586241196 0.0286713280044235],...
            'Orientation','horizontal');
    end
    
    box(subplot1,'on');
    grid(subplot1,'on');
    
    if j==1
        ylabel({'No Variability'}, 'Parent',subplot1);
    elseif  j==4
        ylabel({'Synchrony'}, 'Parent',subplot1);
    elseif j==7
        ylabel({'Anti-synchrony'}, 'Parent',subplot1);
    elseif j==10
        ylabel({'Different Period'}, 'Parent',subplot1);
    end
    
    if j==1
        title('Fishing Mortality','Parent',subplot1)
    elseif j==2
        title('Sardine Biomass','Parent',subplot1)
    elseif j==3
        title('Pelican','Parent',subplot1)
    end
    
    if j<7
        set(subplot1,'XTickLabel',{});
    end
    
   
    if j>6
        set(subplot1,'XTickLabel',{0, 10, 20, 30})
        xlabel('Time','Parent',subplot1);
    end
    
    
end

% Looking at the hockey stick analysis
for VarFlag=1:4
    EsHS=interp1(ddetime,OA_DYNAMICS{VarFlag,9}(9,:)', time);
    XsHS=interp1(ddetime,OA_DYNAMICS{VarFlag,9}(3,:)', time);
    NpHS(:,VarFlag)=interp1(ddetime,OA_DYNAMICS{VarFlag,9}(7,:)', time);
    RentSHS(VarFlag,:)=P(2)*q(2)*EsHS.*(XsHS.^BETA(2))-c0(2)*(EsHS)-c1(2)*(EsHS.^2);
    f = fittype('pchipinterp');
    g1s=fit(time,exp(-dis.*time).*RentSHS(VarFlag,:)',f);
    int1s=integrate(g1s,0:1:time(end),0);
    NPV_SHS(VarFlag)=int1s(end);
    SumNpSHS(VarFlag)=sum(NpHS(:,VarFlag));
end

% Catch control rule
for VarFlag=1:4
    EsHCR=interp1(ddetime,OA_DYNAMICS{VarFlag,10}(9,:)', time);
    XsHCR=interp1(ddetime,OA_DYNAMICS{VarFlag,10}(3,:)', time);
    NpHCR(:,VarFlag)=interp1(ddetime,OA_DYNAMICS{VarFlag,10}(7,:)', time);
    RentSCR(VarFlag,:)=P(2)*q(2)*EsHCR.*(XsHCR.^BETA(2))-c0(2)*(EsHCR)-c1(2)*(EsHCR.^2);
    f = fittype('pchipinterp');
    g1s=fit(time,exp(-dis.*time).*RentSCR(VarFlag,:)',f);
    int1s=integrate(g1s,0:1:time(end),0);
    NPV_SHCR(VarFlag)=int1s(end);
    SumNpHCR(VarFlag)=sum(NpHCR(:,VarFlag));
    
end

%% How binding is the pelican constraint?
clear plot1
for PELICAN=5:5
    figure1=figure
    clear Z
    for VARFLAG=1:3
        
        subplot2 = subplot(1,3,VARFLAG,'Parent',figure1);
        for CASE=1:8
            if CASE<8
                Z(:,CASE)=Npnop{VARFLAG,1,CASE,i}(trun)./INITCOND(7);
            else
                Z(:,CASE)=NpHS(:,VARFLAG)./INITCOND(7);
            end
            
        end
        
        plot1=plot(time,Z(:,1),time,Z(:,2),time,Z(:,3),time,Z(:,4),time,Z(:,5),...
            time,Z(:,6),time,Z(:,7),time,Z(:,8),'Parent',subplot2);
        
        
        set(plot1(1),'DisplayName','A','Color',CC(1,:),'LineWidth',2,'LineStyle','--'); %,'MarkerSize',5,'Marker','o');
        set(plot1(2),'DisplayName','S','Color',CC(2,:),'LineWidth',2,'LineStyle',':'); %,'MarkerSize',5,'Marker','+');
        set(plot1(3),'DisplayName','H','Color',CC(3,:),'LineWidth',2,'LineStyle','-.'); %,'MarkerSize',5,'Marker','s');
        %set(plot1(1),'DisplayName','A','Color',C3,'LineWidth',2,'LineStyle','--');
        %set(plot1(2),'DisplayName','S','Color',C4,'LineWidth',2,'LineStyle',':');
        %set(plot1(3),'DisplayName','H','Color',C1,'LineWidth',2,'LineStyle','-.');
        set(plot1(4),'DisplayName','A+S','Color',CC(4,:),'LineWidth',2,'LineStyle','-');
        set(plot1(5),'DisplayName','A+H','Color',CC(5,:),'LineWidth',2,'LineStyle','-');
        set(plot1(6),'DisplayName','S+H','Color',CC(6,:),'LineWidth',2,'LineStyle','-');
        set(plot1(7),'DisplayName','A+S+H','Color',CC(7,:),'LineWidth',2,'LineStyle','-');
         set(plot1(8),'DisplayName','S_HS','Color',rgb('Black'),'LineWidth',2,'LineStyle','-');
        axis(subplot2,[0 30 min(min(Z)) max(max(Z))])
        
        
        %             legend1 = legend([plot1(1) plot1(2) plot1(3) plot1(4) plot1(5) plot1(6) plot1(7)],...
        %                 {'A','S', 'H', 'A+S','A+H','S+H','A+S+H'});
        %             set(legend1,...
        %                 'Position',[0.36 0.013 0.326 0.036],...
        %                 'Orientation','Vertical');
        %title( 'Pelican ', 'Parent',subplot2)
        
        if VARFLAG==1
            ylabel('Percent increase in P(t)','Parent',subplot2);
        end
        %if VARFLAG>2
        xlabel('Time','Parent',subplot2);
        %end
        %if VARFLAG<3
        %    set(subplot2,'XTickLabel',{});
        %end
        
        if VARFLAG==3
            legend1 = legend(subplot2,'show');
            set(legend1,...
                'Position',[0.164841619727545 0.0116708478226733 0.72492400590772 0.032958198590294],...
                'Orientation','horizontal');
        end
        if VARFLAG==1
            title('No variability')
        elseif  VARFLAG==2
            title('Synchrony')
        elseif VARFLAG==3
            title({'Anti-synchrony'})
        else
            title({'Different period','same starting phase'})
        end
    end
    
end
%% Tables for CC analysis

 clear Tcross Tchange Delta_NPV Delta_Ptotal NPV_P

 
for Varflag=1:4
        
   
        for PELICAN=5
        for CASE=[1 2 4 5 6 7]
            YY(CASE)=(PelicanAVG(Varflag)-Opt_pelicanT{Varflag,CASE});
                if YY(CASE)>0
                    if PELICAN==2
                        ConsP=.05*YY(CASE)+Opt_pelicanT{Varflag,CASE};
                    elseif PELICAN==3
                        ConsP=.1*YY(CASE)+Opt_pelicanT{Varflag,CASE};
                    elseif PELICAN==4
                        ConsP=.15*YY(CASE)+Opt_pelicanT{Varflag,CASE};
                    elseif PELICAN==5
                        ConsP=.2*YY(CASE)+Opt_pelicanT{Varflag,CASE};
                    end
 
                end
            Delta_NPV(CASE,Varflag)=100*(NPV{Varflag,1,CASE,i}-NPV{Varflag,PELICAN,CASE,i})./NPV{Varflag,1,CASE,i};
            NPV_P(CASE,Varflag)=(NPV{Varflag,1,CASE,i}-NPV{Varflag,PELICAN,CASE,i})./(sum(Npnop{Varflag,PELICAN,CASE,i})-sum(Npnop{Varflag,1,CASE,i}));
            Delta_Ptotal(CASE,Varflag)=100*(sum(Npnop{Varflag,PELICAN,CASE,i})-sum(Npnop{Varflag,1,CASE,i}))./sum(Npnop{Varflag,1,CASE,i});
            
            if isempty(min(Time(find(ConsP - Npnop{Varflag,1,CASE,i}(:)>0))))==1
                Tcross(CASE,Varflag)=-99;
            else
                Tcross(CASE,Varflag)=min(Time(find(ConsP - Npnop{Varflag,1,CASE,i}(:)>0)));
            end
            
            temp=100*(Npnop{Varflag,PELICAN,CASE,i}(:)-Npnop{Varflag,1,CASE,i}(:))./Npnop{Varflag,1,CASE,i}(:);
            if isempty(min(Time(find(temp>.1))))==1
                Tchange(CASE,Varflag)=-99
            else
            Tchange(CASE,Varflag)=(min(Time(find(temp>.1))));
            end
            
        end
        end
        
        % Calculating similar measures for the catch control rules for
        % Sardine optimally and
        for CASE=[2 7]
        Delta_Ptotal_hcr(CASE,Varflag) = 100*(SumNpHCR(VarFlag)-sum(Npnop{Varflag,1,CASE,i}(trun)))./sum(Npnop{Varflag,1,CASE,i}(trun));
        Delta_Ptotal_hs(CASE,Varflag) = 100*(SumNpSHS(VarFlag)-sum(Npnop{Varflag,1,CASE,i}(trun)))./sum(Npnop{Varflag,1,CASE,i}(trun));
         dNPV_P_hcr(CASE,Varflag)=  100*( NPV{Varflag,1,CASE,i}-NPV_SHS(VarFlag))./NPV{Varflag,1,CASE,i};
         dNPV_P_hs(CASE,Varflag) =  100*( NPV{Varflag,1,CASE,i}-NPV_SHCR(VarFlag))./NPV{Varflag,1,CASE,i};
         Elas_hcr(CASE,Varflag) = dNPV_P_hcr(CASE,Varflag)./Delta_Ptotal_hcr(CASE,Varflag); 
         Elas_hs(CASE,Varflag) = dNPV_P_hs(CASE,Varflag)./Delta_Ptotal_hs(CASE,Varflag); 
        end
        
end


mgmtregime={'A';'S';'H';'A+H';'A+S';'S+H';'A+S+H'};
Table_results=[Tcross(:,1), Tchange(:,1), Delta_NPV(:,1)./Delta_Ptotal(:,1)];
Latex_Table=latex(vpa(sym(Table_results),2));
clear tempcell
for JJ=1:7
tempcell{JJ}=[mgmtregime(JJ),  Tcross(JJ,1), Tchange(JJ,1), Delta_NPV(JJ,1)./Delta_Ptotal(JJ,1),Tcross(JJ,2), Tchange(JJ,2), Delta_NPV(JJ,2)./Delta_Ptotal(JJ,2)...
    Tcross(JJ,3), Tchange(JJ,3), Delta_NPV(JJ,3)./Delta_Ptotal(JJ,3)];
end
for JJ=1:7
tbl2latex(tempcell{JJ})
end

vregime={'No Variability';'Synchronous';'Anti-Synchronous'};
Table_results=[dNPV_P_hcr([2 7],1:3)', Elas_hcr([2 7],1:3)' dNPV_P_hs([2 7],1:3)', Elas_hs([2 7],1:3)'];
Latex_Table=latex(vpa(sym(Table_results),2));
clear tempcell
for JJ=1:3
tempcell{JJ}=[vregime(JJ),  Table_results(JJ,1) Table_results(JJ,2) Table_results(JJ,3) Table_results(JJ,4) Table_results(JJ,5) Table_results(JJ,6) Table_results(JJ,7) Table_results(JJ,8)];
end
for JJ=1:3
tbl2latex(tempcell{JJ})
end


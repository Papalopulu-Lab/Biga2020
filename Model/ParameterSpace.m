% HesNotchModel.m Time simulation of a multicellular Hes 5 via a (not explicitly modelled) Notch/Delta
% interaction. Number of cells can be defined by rows and cols parameters,
% and will produce a hexagonal grid of cells that each interact with the
% six surrounding cells, with the option of periodic or non-periodic 
% boundary conditions.
function ParameterSpace()
clear;clc;
%==========================================================================
%%             Define simulation conditions and grid size
%==========================================================================
rows=26;         % Defines how many rows of cells to simulate
cols=10;         % Defines how many columns of cells to simulate
cells=rows*cols; % Total number of cells

Stochastic            = 1;    % 1 = stochastic, 0 = deterministic   
Boundary              = 1;    % 0 = 'hard' boundary, 1 = periodic boundary, 2 = 'soft' boundary 
TurnOffAutorepression = 0;    % 0 = full model with Hes autorepression, 1 = model without Hes autorepression
ImplementSwapping     = 0;
CrudeDiff             = 0;
Replace               = 0;
AvM                   = 0;

if cells==1
    Boundary=0; % Other boundaries don't make sense for 1 cell!
end
%==========================================================================
%%                      Defining simulation time
%==========================================================================

t0=0;           % Start time
tf_hours=100;   % Final time (hours)
tf=tf_hours*60; % Final time (min)
dt=3;           % Time step size

if Stochastic==1
    dt=1;       % Stochastic uses Euler-Maruyama method rather than Runge-Kutta fourth order, so this requires a smaller step size      
end

Nt=(tf-t0)/dt;  % Number of time elements
T=t0:dt:tf;     % Time vector

%==========================================================================
%%                          Model parameters
%==========================================================================
HH=4;   %Number of points to plot in parameter space - Hill Coefficient (Z-axis)
II=20;  %Number of points to plot in parameter space - Repression threshold(y-axis)
JJ=20;  %Number of points to plot in parameter space - Time delay (x-axis)
KK=1;   %Number of times to repeat a simulation per parameter point with random initial conditions

% hh_arr=linspace(1,5,HH);     %Hill coefficient
hh_arr=[1,2,4,6];
% hh_arr=4;                      % Hill coefficient
ii_arr=linspace(500,60000,II); %Repression threshold for NICD
jj_arr=linspace(0,200,JJ);    %Time delay between cells
% jj_arr=600;
% hh_arr=4.3;
% ii_arr=500;
% jj_arr=128.6;

wl_mins=200;   %Window length in mins to take the moving mean from
wl=wl_mins/dt; %Window length converted to number of vector elements

load accepted_parameters
load summary_stats
AP=accepted_parameters;
SS=summary_stats;
s=3700;
Prmt=[AP(s,1),  AP(s,2),  AP(s,3), AP(s,5), AP(s,4),   500, 4,   150,  log(2)/30, log(2)/90]; %Coher=0.96, very weakly coupled, clusters of antiphase
% Prmt=[1, 30, 25000, 3, 30, 15000, 4, 128,  log(2)/30, log(2)/90]; %'Regular' param set

parfor_progress(HH*II);
Pk2PkPer=zeros(1,cells);

%==========================================================================
%%                            Main loop(s)
%==========================================================================

tic
for hh=1:HH
    for ii=1:II    
        clc
        parfor_progress;
        parfor jj=1:JJ
            
            VertProtrusions=0;
            %% Produce neighbour matrix
            %The neighbour function produces a 'neighbour matrix' (NM) which is used to
            %multiply out correct Hes values from the protein matrix (P) for each cell.

            [NM,NumNeigh]=neighbours(rows,cols,Boundary,VertProtrusions); 
            NM=sparse(NM);
            eps=(1./NumNeigh);
%             eps=0; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    %         if eps==0 && ii.*jj==1
    %             fprintf('HEADS UP! THE CELLS ARE CURRENTLY UNCOUPLED (eps=0), IS THIS REALLY WHAT YOUR HEART DESIRES?')
    %             pause(5)
    %         end
            %% Parameters
    %         Prmt=[58.5,  1.18,  27378, 6,    40,   ii_arr(ii), 1.5,   jj_arr(jj),  log(2)/30, log(2)/90]; %Coher=0.96, very weakly coupled, clusters of antiphase

            a_m   = Prmt(1); a_p   = Prmt(2); P_H0  = Prmt(3); n_H   = Prmt(4);     
            TauH  = Prmt(5); P_ND0 = ii_arr(ii); n_ND = hh_arr(hh); TauND = jj_arr(jj);
            u_m   = Prmt(9); u_p   = Prmt(10);

            TauH_step=round(TauH/dt);     % Conversion to simulation time steps
            TauND_step=round(TauND/dt); % Conversion to simulation time steps
            gamma=1;

%             a_m=1;
%             a_p=30;
%             u_m=log(2)/30;
%             u_p=log(2)/90; % Protein degradation rate (1/min)
%             P_H0=25000;
%             P_ND0=ii_arr(ii); % Repression threshold for NICD effects/neighbouring cells (units of concentration) (=100 is interesting)
%             n_H=3;         % Hill coefficient (intracellular)
%             n_ND=hh_arr(hh);        % Hill coefficient (intercellular)
%             TauH=30;
%             TauND=jj_arr(jj);    



            DiffTime_hours=50;         %Time at which differention can start to occur (hours)
            DiffTime=DiffTime_hours*60;
            S=0.0000003; % Multiplier for crude diff to tune the rate of differentiation.

            if CrudeDiff==0
                DiffTime=Nt*2; %ensures crude diff doesnt occur
            end

            %Cell swap/migration parameters
            Pm=0.5; % Total probability of movement either left or right
            SwapThresh=1-Pm/2; % Treshold to be used in probability of swapping function within diffSolver.m
            Ts=1; % Time (in mins) between each swapping event

            %% DDEs (using anonymous functions (@ symbol) rather than standard functions)
            dm=@(m,p,p_delay,Psum_delay,gamma) a_m*1./(1+(p_delay./P_H0).^n_H).*gamma./(1+(eps.*Psum_delay./P_ND0).^n_ND) - u_m.*m; % Describe mRNA change in time


            if TurnOffAutorepression==1
                dm=@(m,p,p_delay,Psum_delay,gamma) a_m.*gamma./(1+(eps.*Psum_delay./P_ND0).^n_ND) - u_m.*m; % No self repression version for simple Notch Delta modelling
            end

            dp=@(m,p,p_delay,Psum_delay,gamma) a_p.*m - u_p.*p; % Describes protein change in time

            %% SDDEs
            dm1=dm;

            dm2=@(m,p,p_delay,Psum_delay,gamma) sqrt(a_m.*gamma./(1+(p_delay./P_H0).^n_H)*1./(1+(eps.*Psum_delay./P_ND0).^n_ND) + u_m.*m);

            if TurnOffAutorepression==1
               dm2=@(m,p,p_delay,Psum_delay,gamma) sqrt(a_m.*gamma./(1+(eps.*Psum_delay./P_ND0).^n_ND) + u_m.*m);
            end

            dp1=dp;
            dp2=@(m,p,p_delay,Psum_delay,gamma) sqrt(a_p.*m + u_p.*p);
            
%--------------------------------------------------------------------------
            autoCorr_kk       = zeros(1,KK);
            avgPXX2_period_kk = zeros(1,KK);
            clusterMean_kk    = zeros(1,KK);
            clusterOcc_kk     = zeros(1,KK);
            Coherence_kk  = zeros(1,KK);
            fourierPer_kk = zeros(1,KK);
            fourierPerFilt_kk = zeros(1,KK);
            Sync_val_kk   = zeros(1,KK);
            BSf_kk        = zeros(1,KK); 
            LSf_kk        = zeros(1,KK);
            Ff_kk         = zeros(1,KK);
            bsPerM_kk     = zeros(1,KK);
            bsPerSD_kk    = zeros(1,KK);
            lsPerM_kk     = zeros(1,KK);
            lsPerSD_kk    = zeros(1,KK);
            fPerM_kk      = zeros(1,KK);
            fPerSD_kk     = zeros(1,KK);
            meanP_kk      = zeros(1,KK);
            rangePop_kk   =zeros(1,KK);
            rangeSingle_kk=zeros(1,KK);
            
            for kk=1:KK %This loop is to generate sevral simulations to take an average over

%==========================================================================
%%                       Main loop (solves DDEs)
%==========================================================================
                if KK==1
                    rng(3);
                end
                
                P=[rndrng(cells,1,0.2*P_H0,1.8*P_H0),zeros(cells,Nt-1),zeros(cells,1)]; %Vector that will store protein values from the simulation
                M=[rndrng(cells,1,0,20),zeros(cells,Nt-1),zeros(cells,1)];

                PercentCounter=0;
                [P, ~, t_step, DiffYN, DiffYNflash, CT, DiffThresh, MM]=diffSolver(P, M, Nt, TauND_step, TauH_step, NM, gamma, dm, dp, dm1, dm2, dp1, dp2, dt, Stochastic, rows, cols, DiffTime, S, ImplementSwapping, SwapThresh, Ts, PercentCounter, AvM, Replace, wl);

                add=meshgrid(0:cells:cells*(Nt+1)-cells, 1:cells);
                [~,sortMatrix]=sort(CT);
                recoverCellPos=sortMatrix+add;

                P=P(recoverCellPos);

%==========================================================================
%%                     Temporal Fourier transform 
%==========================================================================
                
                
                Osc= max(P(:,0.8*Nt:Nt),[],2)>5000; %Useful in salt and pepper for removing non oscillating cells
            %         Osc= min(P(:,0.5*Nt:Nt),[],2)>=0;

            
%                 figure(15)
%                 clf
%                 subplot(211)
%                 plot(T./60,P)
%                 subplot(212)
%                 plot(T./60,P(5,:));hold on
%      
%                 drawnow
                
                
                if Stochastic==1
                    frac=0.6;
                else   
                    frac=0.01;
                end

                Y_raw=P(Osc,:);

                %Detrending
                polyOrder=3;                  %Order of detrending polynomial in detrend.m function
                frameTime=30;                 %Frame length in hours for detrending window
                frameLength=frameTime*60./dt; %Conversion to window length in elements  

                [Ydetrend, t, Ysmooth,f,P1,coherence,f_C1,P_C1,avgFourier,ind_per,I]=tempFourier(T,Y_raw,polyOrder,frac,frameLength,Nt,dt);

%                 figure(16)
%                 clf
%                 cell=3;
%                 subplot(311)
%                 plot(t./60, Ysmooth(cell,:)+Ydetrend(cell,:));hold on
%                 plot(t./60,Ysmooth(cell,:))
%                 subplot(312)
% %                 plot(60*f,mean(P1,1))
%                 plot(60*f,P1(:,:))
%                 xlim([0 0.4])
%                 subplot(313)
% %                 plot(60*f,mean(P1,1))
%                 plot(60*f,P1(cell,:))
%                 xlim([0 0.4])
%                 
                
                Coherence_kk(kk)=coherence; 
                fourierPer_kk(kk)=1./(60*f(I));
                
                fourierPerFiltPre=ind_per;
                fourierPerFiltPre(ind_per>10)=[];
                if numel(fourierPerFiltPre)>1
                    fourierPerFilt_kk(kk)=mean(fourierPerFiltPre);
                else
                    fourierPerFilt_kk(kk)=fourierPer_kk(kk);
                end


%==========================================================================
%%                      Kuramoto order parameter
%==========================================================================            
                
                [PH,PH_unwrap,COP]=phase(t,Ydetrend); % Uses Hilbert transform to extract phase from oscillatory signals
                Sync_val_kk(kk)=mean(abs(COP(round(0.5*size(Ydetrend,2)):round(0.9*size(Ydetrend,2)))));
%                 figure(1)
%                 clf
%                 subplot(211)
%                 plot(t,PH_unwrap)
%                 subplot(212)
%                 plot(t,abs(COP))
%                 drawnow
                
                
% %==========================================================================                
% %%                          Auto-correlation
% %==========================================================================
% 
                startTime=60; %Removes signal before this time for the processing
                time_frac=startTime/tf_hours;
                t=T(time_frac*Nt:end)./60;
                ti=time_frac*Nt;
% 
                Y_RAW=(P((cols/2-1)*rows+1:cols/2*rows,ti:end)+P(cols/2*rows+1:(cols/2+1)*rows,ti:end))./2;
%             %     Y_RAW=(P((cols/2-1)*rows+1:cols/2*rows,time_frac*Nt:end)+P(cols/2*rows+1:(cols/2+1)*rows,time_frac*Nt:end)+P((cols/2+1)*rows+1:(cols/2+2)*rows,time_frac*Nt:end))./3; %Three Cols
%             %     Y_RAW=P((cols/2-1)*rows+1:cols/2*rows,time_frac*Nt:end); %single column
% 
% %                 YDETREND=Y_RAW-repmat(mean(Y_RAW,1),[rows 1]); %Detrending in spatial direction by subrtracing poulation mean at each time point
%                 YDETREND=Y_RAW-mean(Y_RAW(:));
%                 
%                 NN=numel(Y_RAW(1,:));
%                 AVGPeriod=zeros(NN,1);
%                 for nn=1:NN
%                     y=YDETREND(:,nn)';
%                     Ny=length(y);
%                     S=1:Ny;
%                     [CorrSignal,AvgPeriod]=autocorrelation2(y,1);
%                     AVGPeriod(nn)=AvgPeriod;
%                 end
% 
%                 autoCorr_kk(kk)=nanmean(AVGPeriod);
                

%==========================================================================
%%                 Spatial power spectrum significance
%==========================================================================  
                kymoWidth=2; %Number of cells to use for kymograph selection width
                K=cols/kymoWidth;

                if mod(cols,2)==1
                    error('Make cols even')
                end
                
                Yraw=zeros(size(Y_RAW,1),size(Y_RAW,2),K);
                Y_RAW=zeros(size(Y_RAW,1),size(Y_RAW,2)*K);
                
                for k=1:K

                    Yraw(:,:,k)=(P(2*(k-1)*rows+1:(2*k-1)*rows,ti:end)+P((2*k-1)*rows+1:2*k*rows,ti:end))./2;

                    if k==1
                        Y_RAW=Yraw(:,:,1);
                    else
                        Y_RAW=[Y_RAW Yraw(:,:,k)];
                    end

                end
                
%                 figure(1)
%                 clf
%                 
%                 subplot(121)
%                 Y=1:rows;
%                 j=surf(Y_RAW);
%                 title('Kymograph')
% %                 ylim([Y(1) Y(end)])
% %                 xlim([t(1) t(end)])
%                 xlabel('Time (hrs)')
%                 ylabel('Cell index/row number')
%                 set(j, 'LineStyle','none')
%                 set(gca,'YTickLabel',[]);
%                 set(gca,'FontSize',15)
%                 colormap(gca,'jet')
%                 colorbar
%                 view(0,90)
                
%                 subplot(122)
%                 plot(T,P(1:60,:))
%                 drawnow
%                 
                
                
                %Split Y_RAW into x hour windows of expression and then average
                t_elem=(1-time_frac)*Nt*K;
                split_time=0;
                if split_time==0
                    split_elem=1;
                else
                    split_elem=round(split_time*60/dt);
                end

                M=1:split_elem:t_elem;
                M=round(M);
                I=length(M);

                split_t=linspace(0,tf*K,I);

                Y_split=Y_RAW(:,M);            
                YDETREND2=Y_split-repmat(mean(Y_split,1),[rows 1]); %Detrending in spatial direction by subtracing poulation mean at that time point
%                 YDETREND2=Y_split-mean(Y_RAW(:));
                
                Y=fft(YDETREND2,[],1);
                L=length(Y(:,1));
                P2 = abs(Y/L).^2;
                P1 = P2(1:L/2+1,:);
                P1(2:end-1,:) = 2*P1(2:end-1,:);

                Fs=1;
                f = Fs*(0:(L/2))/L;

                KymLong=YDETREND2;
                zs=zscore(KymLong); %z-scored data
            %     [PXX2,F2]=periodogram(zs,[],[],1);         %Run periodogram and boostrap intervals
                PXX2=P1(1:end,:);
                F2=f(1:end);
                
                avgPXX2=mean(PXX2,2);
                [~,J]=max(avgPXX2);
                avgPXX2_period_kk(kk)=1/F2(J);
    
          
% %--------------------------------------------------------------------------
% %                    Frequency analysis and stats
% %--------------------------------------------------------------------------
% 
                zs=zscore(KymLong); %z-scored data
                [PXX1,F1,Pth]=plomb(zs,1,'Pd',0.95,'psd'); %Run Lomb Scargle
%                 [PXX2,F2]=periodogram(zs,[],[],1);         %Run periodogram and boostrap intervals
                [PXXbs,Fbs,BSth]=GetBootstrapPS(zs);
%                 
%                 
                I=size(KymLong,2);

                %Preallocation for storing stats
                LSPperiod=zeros(I,1);
                fisherG=zeros(I,1);
                pval=zeros(I,1);
                peakLoc=zeros(I,1);
                peakHeight=zeros(I,1);
                period=zeros(I,1);

                for i=1:I
                    % find if peak is significant for Lomb Scargle
                    pxx1=PXX1(:,i);
                    idx1=find(pxx1==max(pxx1));
                    if pxx1(idx1)>Pth(i) % get period if peak is significant
                        LSPperiod(i) = 1/F1(idx1);
                    else
                        LSPperiod(i)=nan;
                    end

                    pxx=PXX2(:,i);

                    [fisherG(i),pval(i),idx]=GetFisherG(pxx); % Find the peak and collect Fisher G-statistic
                    peakLoc(i)=idx;
                    peakHeight(i)=pxx(idx);
                    if pval(idx)<0.05
                        period(i)=1/F2(idx);
                    else
                        period(i)=nan;
                    end
                    
                end
                Avg_Sig_per=nanmean(period);
                
%                 figure(102)
%                 clf
%                 big=find(period>30);
%                 subplot(211)
%                 plot(KymLong(:,big(2)))
%                 subplot(212)
%                 plot(F2,PXX2(:,big(2)))

%                 remove=find(period>20);
%                 period(remove)=NaN;
%                
%                 pval(remove)=1;

%     figure(31)
%     clf
%     
%     subplot(221)
%     X=t;
%     Y=1:rows;
%     j=surf(X,Y,Y_RAW);
%     title('Kymograph')
%     ylim([Y(1) Y(end)])
%     xlim([t(1) t(end)])
%     xlabel('Time (hrs)')
%     ylabel('Cell index/row number')
%     set(j, 'LineStyle','none')
%     set(gca,'YTickLabel',[]);
%     set(gca,'FontSize',15)
%     colormap(gca,'jet')
%     colorbar
%     view(0,90)
%     pbaspect([2 1 1])
%     
%     subplot(222)
%     X=split_t;
%     Y=1:rows;
%     j=surf(X,Y,YDETREND2);
%     title('Detrended kymograph')
%     ylim([Y(1) Y(end)])
%     xlim([t(1) t(end)])
%     xlabel('Time (hrs)')
%     ylabel('Cells/row number')
%     set(j, 'LineStyle','none')
%     set(gca,'YTickLabel',[]);
%     set(gca,'FontSize',15)
%     colormap(gca,'jet')
%     colorbar
%     view(0,90)
%     pbaspect([2 1 1])
%     
%     
%     subplot(223)
%     load cmap
%     colormap(gca,cmap)
% %     [map2]=colourMap2(0.4);
% %     colormap(gca,map2)
%     t_plot=linspace(time_frac*tf/60,tf/60,numel(P1(1,:)));
%     h=surf(t_plot,F2,PXX2);
%     xlabel('Time (hours)')
%     ylabel('Frequency (1/cell)')
%     title('Spatial Fourier transform over time')
%     xlim([t(1) t(end)])
%     ylim([f(1) f(end)])
%     view(0,90)
%     colorbar
%     set(h,'LineStyle','none')
%     set(gca,'FontSize',15)
%     pbaspect([2 1 1]);
%     
%     subplot(224)
%     plot(t_plot,period,'color','w','linewidth',2)
%     title(sprintf('Dominant spatial period: %.1f', nanmean(period)))
% %     title('Significant period (Fisher G)')
%     ylabel('Period (cells)')
%     xlabel('Time (hours)')
%     set(gca,'FontSize',15)
%     xlim([t(1) t(end)])
% %     ylim([0 10])
%     drawnow

%% Cluster detection
    thresh=0.8;
    l=zeros(I,1);
    for i=1:I
        vect=KymLong(:,i);
        [l(i),~]=microClustDetect(vect,thresh);
    end
    clusterMean_kk(kk)=nanmean(l);
    clusterOcc_kk(kk)=1-sum(isnan(l))/numel(l);
%     figure(13),histogram(l,'normalization','probability'),ylabel('Frequency/occurence'),xlabel('Cluster size (radius)');

%--------------------------------------------------------------------------
%                     Storing stats for param space
%--------------------------------------------------------------------------
                
                BSf_kk(kk)=sum(peakHeight>BSth)/I;   %Fraction of power spectrums peaks above bootstrap threshold
                LSf_kk(kk)=1-sum(isnan(LSPperiod))/I;
                Ff_kk(kk)=sum(pval<0.05)/I;

                bsPerM_kk(kk)    = nanmean(period);
                bsPerSD_kk(kk)   = nanstd(period);
                lsPerM_kk(kk)  = nanmean(LSPperiod);
                lsPerSD_kk(kk) = nanstd(LSPperiod);
                fPerM_kk(kk)=nanmean(period(pval<0.05));
                fPerSD_kk(kk)=nanstd(period(pval<0.05));
                meanP_kk(kk)=mean(Y_RAW(:));
                rangePop_kk(kk)=range(mean(P(:,0.5*Nt:end),2));
                rangeSingle_kk(kk)=mean(range(P(:,0.5*Nt:end),2));
                
                
            end
            autoCorr(hh,ii,jj)   = nanmean(autoCorr_kk);
            avgPXX2_period(hh,ii,jj) = nanmean(avgPXX2_period_kk);
            clusterMean(hh,ii,jj) = nanmean(clusterMean_kk);
            clusterOcc(hh,ii,jj)  =nanmean(clusterOcc_kk);
            Coherence(hh,ii,jj)  = nanmean(Coherence_kk); 
            fourierPer(hh,ii,jj) = nanmean(fourierPer_kk);
            fourierPerFilt(hh,ii,jj)=nanmean(fourierPerFilt_kk);
            Sync_val(hh,ii,jj)   = nanmean(Sync_val_kk);
            BSf(hh,ii,jj)        = nanmean(BSf_kk);   %Fraction of power spectrums peaks above bootstrap threshold
            LSf(hh,ii,jj)        = nanmean(LSf_kk);
            Ff(hh,ii,jj)         = nanmean(Ff_kk);
            bsPerM(hh,ii,jj)     = nanmean(bsPerM_kk); %Bootstrap period mean
            bsPerSD(hh,ii,jj)    = nanmean(bsPerSD_kk);
            lsPerM(hh,ii,jj)     = nanmean(lsPerM_kk); %Lomb-Scargle period mean
            lsPerSD(hh,ii,jj)    = nanmean(lsPerSD_kk);
            fPerM(hh,ii,jj)      = nanmean(fPerM_kk);  % Fisher period mean
            fPerSD(hh,ii,jj)     = nanmean(fPerSD_kk); 
            meanP(hh,ii,jj)      = mean(meanP_kk);
            rangePop(hh,ii,jj)   = mean(rangePop_kk);
            rangeSingle(hh,ii,jj)= mean(rangeSingle_kk);
            
        end
    end
end
toc

%==========================================================================
%%                               Plots
%==========================================================================
GraphAppearance=0;
% close all;
if GraphAppearance==1
    
    BackgroundColour=[64/255 64/255 64/255];
    TextColour=[1 1 1];
    get(0,'Factory');
    set(0,'defaultfigurecolor',BackgroundColour)
    set(0,'DefaultAxesFontSize', 10)
    set(0,'defaultAxesColor',BackgroundColour)
    set(0,'defaultAxesZColor',TextColour)
    set(0,'defaultAxesXColor',TextColour)
    set(0,'defaultAxesYColor',TextColour)
    set(0,'defaultLegendTextColor',TextColour)
    set(0,'defaultTextColor',TextColour)
    
elseif GraphAppearance==0

    get(0,'Factory');
    set(0,'defaultfigurecolor',[0.94 0.94 0.94])
    set(0,'DefaultAxesFontSize', 10)
    set(0,'defaultAxesColor','w')
    set(0,'defaultAxesXColor','k')
    set(0,'defaultAxesYColor','k')
    set(0,'defaultLegendTextColor','k')
    set(0,'defaultTextColor','k')
end

textSize=10;
colbar=1;
%==========================================================================
%                        Stats parameter space
%=========================================================================
% figure(1)
% clf
% set(gcf,'renderer','Painters') %For EPS file export
% for hh=1:HH
%     
%     subplot(4,HH,hh)
% %     DATA=fourierPer;
%     DATA=fourierPerFilt;
%     data=squeeze(DATA(hh,:,:));
%     surf(jj_arr,ii_arr,data);
%     title({sprintf('n_{ND}=%.1f',hh_arr(hh)) ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (min)'); 
%     if hh==1
%         ylabel({'Repression','threshold'});
%     end
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     if colbar==1
%         colorbar
%     end
%     caxis([min(DATA(:)) max(DATA(:))])
%     colormap(magma(1000))
%     
%     subplot(4,HH,hh+1*4)    
%     DATA=Sync_val;
%     data=squeeze(DATA(hh,:,:));
%     surf(jj_arr,ii_arr,data);
% %     title({sprintf('n_{ND}= %.1f',hh_arr(hh)) ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (min)');
%     if hh==1
%         ylabel({'Repression','threshold'});
%     end
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     if colbar==1
%         colorbar
%     end
%     caxis([0 1])
%     colormap(magma(1000))
%      
%     
%     
%     subplot(4,HH,hh+3*4)    
%     DATA=clusterOcc;
%     data=squeeze(DATA(hh,:,:));
%     surf(jj_arr,ii_arr,data);
% %     title({sprintf('n_{ND}= %.1f)',hh_arr(hh)) ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (min)');
%     if hh==1
%         ylabel({'Repression','threshold'});
%     end
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     caxis([0 0.5])
%     if colbar==1
%         colorbar
%     end
%     colormap(magma(1000))
%     
% end
% 
% 
% remove=0.15;
% 
% bsPerM2 = bsPerM;
% lsPerM2 = lsPerM;
% FPerM2  = fPerM;
% bsPerM2(BSf<remove)=nan;
% lsPerM2(LSf<remove)=nan;
% FPerM2(Ff<remove)=nan;
% avgPXX2_period_NaN=avgPXX2_period;
% avgPXX2_period_NaN(Ff<remove)=nan;
% avgPXX2_period_NaN(avgPXX2_period>3)=nan;
% 
% figure(1) 
% set(gcf,'renderer','Painters') %For EPS file export
% fig = gcf;
% if GraphAppearance==1
%     fig.InvertHardcopy = 'off';
% elseif GraphAppearance==0
%     fig.InvertHardcopy = 'on';
% end
% 
% for hh=1:HH
%     
%     data=squeeze(avgPXX2_period_NaN(hh,:,:));
%     data=[data(1,:); data(1:end-1,:)];
%     data=[data(:,1), data(:,1:end-1)];
%     
%     subplot(4,HH,hh+2*4)
%     surf(jj_arr,ii_arr,data);
% %     title({sprintf('n_{ND}= %.1f',hh_arr(hh)) ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (min)');    
%     if hh==1
%         ylabel({'Repression','threshold'});
%     end
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     if colbar==1
%         colorbar
%     end
%     caxis([2 4])
%     colormap(magma(1000))
% 
% end
% drawnow
%   
% if HH==1
%     %Main plot
%     figure(112)
%     clf
%     
%     subplot(231) % Temporal period
%     surf(jj_arr,ii_arr,squeeze(fourierPer(hh,:,:)));
%     title('Temporal period');
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (min)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     colorbar
%     h2=colorbar;
%    
%     subplot(232)
%     surf(jj_arr,ii_arr,squeeze(avgPXX2_period(hh,:,:)));
%     title('Spatial period');
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     colorbar
%     
%     subplot(233) %Sync val
%     surf(jj_arr,ii_arr,squeeze(Sync_val(hh,:,:)));
%     title('Sync value');
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     colorbar
%     
%     subplot(234) % Temporal period
%     surf(jj_arr,ii_arr,squeeze(meanP(hh,:,:)));
%     title('Mean expression');
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     colorbar
%     h2=colorbar;
%     
%     subplot(236) % Temporal period
%     surf(jj_arr,ii_arr,squeeze(Coherence(hh,:,:)));
%     title('Coherence');
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     colorbar
%     h2=colorbar;
% end





remove=0.15;

bsPerM2 = bsPerM;
lsPerM2 = lsPerM;
FPerM2  = fPerM;
bsPerM2(BSf<remove)=nan;
lsPerM2(LSf<remove)=nan;
FPerM2(Ff<remove)=nan;
avgPXX2_period_NaN=avgPXX2_period;
avgPXX2_period_NaN(Ff<remove)=nan;
avgPXX2_period_NaN(avgPXX2_period>3)=nan;

figure(10)
clf
set(gcf,'renderer','Painters') %For EPS file export
for hh=1:HH
    
    data=squeeze(clusterOcc(hh,:,:));
    data=[data(1,:); data(1:end-1,:)];
    data=[data(:,1), data(:,1:end-1)];
    
    subplot(2,HH,hh)
    surf(jj_arr,ii_arr,data);
    title({sprintf('n_{ND}= %.1f',hh_arr(hh)) ''});
    xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
    xlabel('Time delay (min)');    
    if hh==1
        ylabel({'Repression','threshold'});
    end
    axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
    if colbar==1
        colorbar
    end
    caxis([0.1 0.5])
    
    data=squeeze(avgPXX2_period_NaN(hh,:,:));
    data=[data(1,:); data(1:end-1,:)];
    data=[data(:,1), data(:,1:end-1)];
    
    subplot(2,HH,hh+1*4)
    surf(jj_arr,ii_arr,data);
%     title({sprintf('n_{ND}= %.1f',hh_arr(hh)) ''});
    xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
    xlabel('Time delay (min)');    
    if hh==1
        ylabel({'Repression','threshold'});
    end
    axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
    if colbar==1
        colorbar
    end
    caxis([2 4])
    colormap(magma(1000))

end





















































%--------------------------------------------------------------------------
% % Temp fourier
% figure(3)
% clf
% set(gcf,'renderer','Painters') %For EPS file export
% for hh=1:HH
%     
%     DATA=fourierPer;
% %     DATA=fourierPerFilt;
%     data=squeeze(DATA(hh,:,:));
%     
%     subplot(1,HH,hh)
% %     surf(jj_arr,ii_arr,squeeze(fourierPer(hh,:,:)));
%     surf(jj_arr,ii_arr,data);
%     title({sprintf('n_{ND}=%.1f',hh_arr(hh)) ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     colorbar
% %     h2=colorbar;
% %     collim=get(h2,'ylim');
% %     caxis([min(fourierPer(:)) max(fourierPer(:))])
%     caxis([min(DATA(:)) max(DATA(:))])
%     colormap(magma(1000))
%     
% end


% figure(4) 
% set(gcf,'renderer','Painters') %For EPS file export
% fig = gcf;
% if GraphAppearance==1
%     fig.InvertHardcopy = 'off';
% elseif GraphAppearance==0
%     fig.InvertHardcopy = 'on';
% end
% clf
% textSize=9;

% %Plot fisher G for all n-values

% for hh=1:HH
%     
%     data=squeeze(Sync_val(hh,:,:));
% %     data=[data(1,:); data(1:end-1,:)];
% %     data=[data(:,1), data(:,1:end-1)];
%     
%     subplot(1,HH,hh)
%     surf(jj_arr,ii_arr,data);
%     title({sprintf('n_{ND}= %.1f',hh_arr(hh)) ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% %     colorbar
%     caxis([0 1])
%     colormap(magma(1000))
% 
%     
% end


% figure(5) 
% set(gcf,'renderer','Painters') %For EPS file export
% fig = gcf;
% if GraphAppearance==1
%     fig.InvertHardcopy = 'off';
% elseif GraphAppearance==0
%     fig.InvertHardcopy = 'on';
% end
% clf
% textSize=9;

% %Coherence
% 
% for hh=1:HH
%     
%     data=squeeze(Coherence(hh,:,:));
% %     data=[data(1,:); data(1:end-1,:)];
% %     data=[data(:,1), data(:,1:end-1)];
%     
%     subplot(1,HH,hh)
%     surf(jj_arr,ii_arr,data);
%     title({sprintf('n_{ND}= %.1f',hh_arr(hh)) ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     colorbar
%     caxis([min(Coherence(:)) max(Coherence(:))])
% 
%     
% end

% figure(6) 
% set(gcf,'renderer','Painters') %For EPS file export
% fig = gcf;
% if GraphAppearance==1
%     fig.InvertHardcopy = 'off';
% elseif GraphAppearance==0
%     fig.InvertHardcopy = 'on';
% end
% clf
% textSize=9;
% 
% for hh=1:HH
%     
%     data=squeeze(meanP(hh,:,:));
% %     data=[data(1,:); data(1:end-1,:)];
% %     data=[data(:,1), data(:,1:end-1)];
%     
%     subplot(1,HH,hh)
%     surf(jj_arr,ii_arr,data);
%     title({sprintf('n_{ND}= %.1f',hh_arr(hh)) ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     colorbar
%     caxis([min(meanP(:)) max(meanP(:))])
% 
%     
% end





% %% Autocorrelation plot
% figure(12)
% clf
% surf(jj_arr,ii_arr,squeeze(autoCorr(hh,:,:)));
% title('Autocorrelation period');
% xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% xlabel('Time delay (mins)');    ylabel('Repression threshold');
% axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% colorbar
% 
% %% Fourier average peak
% figure(13)
% clf
% fig = gcf;
% if GraphAppearance==1
% fig.InvertHardcopy = 'off';
% elseif GraphAppearance==0
% fig.InvertHardcopy = 'on';
% end

% subplot(121)
% surf(jj_arr,ii_arr,squeeze(avgPXX2_period(hh,:,:)));
% title('Average spatial Fourier peak');
% xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% xlabel('Time delay (mins)');    ylabel('Repression threshold');
% axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% colorbar
% 
% % subplot(122)
% surf(jj_arr,ii_arr,squeeze(avgPXX2_period_NaN(hh,:,:)));
% title('Average spatial Fourier peak (non-sig removed)');
% xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% xlabel('Time delay (mins)');    ylabel('Repression threshold');
% axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% colorbar
% caxis([2 6])


% %% Clusters
% figure(14)
% max(clusterMean(:))
% for hh=1:HH
%     
% %     subplot(2,HH,hh)
% %     surf(jj_arr,ii_arr,squeeze(clusterMean(hh,:,:)));
% %     title({sprintf('Cluster size (n_{ND}=%.1f)',hh_arr(hh)) ''});
% %     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% %     xlabel('Time delay (mins)');    ylabel('Repression threshold');
% %     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% %     h1=colorbar;
% %     collim=get(h1,'ylim');
% %     caxis([3 7])
%     
%     
%     subplot(1,HH,hh)
%     surf(jj_arr,ii_arr,squeeze(clusterOcc(hh,:,:)));
%     title({sprintf('n_{ND}= %.1f)',hh_arr(hh)) ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     caxis([0 0.5])
%     colorbar
%     colormap(magma(1000))
% end
% drawnow

% figure(15)
% textSize=9;
% for hh=1:HH
%     subplot(1,HH,hh)
%     surf(jj_arr,ii_arr,squeeze(clusterOcc(hh,:,:)));
%     title({sprintf('n_{ND}=%.1f',hh_arr(hh)) ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     caxis([0 0.5])
%     colorbar
% end
% drawnow
% max(clusterOcc(:))
% 
% 
remove=0.2;
clusterMean_NaN=clusterMean;
clusterMean_NaN(clusterOcc<remove)=nan;
figure(16)
set(gcf,'renderer','Painters')
for hh=1:HH
    subplot(1,HH,hh)
    surf(jj_arr,ii_arr,squeeze(clusterMean_NaN(hh,:,:)));
    title({sprintf('n_{ND}=%.1f',hh_arr(hh)) ''});
    xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
    xlabel('Time delay (mins)');    ylabel('Repression threshold');
    axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
    caxis([3 7])
    colorbar
    colormap(magma(1000))
end
drawnow

remove=0.12;
clusterOcc_NaN=clusterOcc;
clusterOcc_NaN(clusterOcc<remove)=nan;
figure(17)
set(gcf,'renderer','Painters')
for hh=1:HH
    subplot(1,HH,hh)
    surf(jj_arr,ii_arr,squeeze(clusterOcc_NaN(hh,:,:)));
    title({sprintf('n_{ND}=%.1f',hh_arr(hh)) ''});
    xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
    xlabel('Time delay (mins)');    ylabel('Repression threshold');
    axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     caxis([3 7])
    colorbar
    colormap(magma(1000))
end
drawnow




% figure(14)
%     subplot(2,1,1)
%     surf(jj_arr,ii_arr,squeeze(clusterMean(hh,:,:)));
%     title({sprintf('Cluster size (n_{ND}=%.1f)',hh_arr(hh)) ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     h1=colorbar;
% %     collim=get(h1,'ylim');
% %     caxis([2 collim(2)])
%     
%     
%     subplot(2,1,2)
%     surf(jj_arr,ii_arr,squeeze(clusterOcc(hh,:,:)));
%     title({sprintf('Cluster occurence (n_{ND}=%.1f)',hh_arr(hh)) ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% %     caxis([0 1])
%     colorbar
    
    
%     %% Coherence
%     figure(102)
%     clf
%     surf(jj_arr,ii_arr,squeeze(Coherence(hh,:,:)));
%     title({sprintf('Coherence)') ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     h1=colorbar;

    
    
    
    
% figure(11)
% clf
% subplot(121)
% surf(jj_arr,ii_arr,squeeze(Sync_val(hh,:,:)));
% title('Sync value');
% xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% xlabel('Time delay (mins)');    ylabel('Repression threshold');
% axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% colorbar
% 
% subplot(122)
% surf(jj_arr,ii_arr,1-squeeze(Sync_val(hh,:,:)));
% title('Error value');
% xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% xlabel('Time delay (mins)');    ylabel('Repression threshold');
% axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% colorbar
% 
% 
% 
% 
% 
% figure(24)
% clf
% surf(jj_arr,ii_arr,squeeze(fourierPerFilt(hh,:,:)));
% title('Filtered fourier');
% xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% xlabel('Time delay (mins)');    ylabel('Repression threshold');
% axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% colorbar
% colormap(magma(1000))    
    
    
%     figure(1) 
% fig = gcf;
% set(gcf,'renderer','Painters') %For EPS file export
% if GraphAppearance==1
%     fig.InvertHardcopy = 'off';
% elseif GraphAppearance==0
%     fig.InvertHardcopy = 'on';
% end
% clf
% textSize=8;
% 
% %Plot fisher G for all n-values
% 
% for hh=1:HH
%     
%     subplot(2,HH,hh)
% %     surf(jj_arr,ii_arr,squeeze(fPerM(hh,:,:)));
%     surf(jj_arr,ii_arr,squeeze(avgPXX2_period(hh,:,:)));
%     title({sprintf('Spatial period (n_{ND}=%.1f)',hh_arr(hh)) ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     h1=colorbar;
%     collim=get(h1,'ylim');
%     caxis([2 collim(2)])
%     
%     
%     subplot(2,HH,hh+HH)
%     surf(jj_arr,ii_arr,squeeze(Ff(hh,:,:)));
%     title({sprintf('Sig occurence (n_{ND}=%.1f)',hh_arr(hh)) ''});
%     xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
%     xlabel('Time delay (mins)');    ylabel('Repression threshold');
%     axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
%     caxis([0 1])
%     colorbar
% end
% drawnow






















% subplot(234) % Average period
% surf(jj_arr,ii_arr,BS_thresh_frac);
% title({'Fraction that passed BS' ''});
% xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% xlabel('Time delay (mins)');    ylabel('Repression threshold');
% axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% caxis([0 1])
% colorbar
% %--------------------------------------------------------------------------
% 
% subplot(232) % Average period
% surf(jj_arr,ii_arr,lsPerM);
% title({'LS Period mean (hours)' ''});
% xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% xlabel('Time delay (mins)');    ylabel('Repression threshold');
% axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% colorbar
% 
% subplot(235) % Average period
% surf(jj_arr,ii_arr,LS_thresh_frac);
% title({'Fraction that passed LS' ''});
% xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% xlabel('Time delay (mins)');    ylabel('Repression threshold');
% axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% caxis([0 1])
% colorbar
% %--------------------------------------------------------------------------
% subplot(233) % Average period
% surf(jj_arr,ii_arr,FPerM);
% title({'Fisher G period' ''});
% xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% xlabel('Time delay (mins)');    ylabel('Repression threshold');
% axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% colorbar
% 
% subplot(236) % Average period
% surf(jj_arr,ii_arr,F_sig_frac);
% title({'Fisher G above 0.05' ''});
% xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% xlabel('Time delay (mins)');    ylabel('Repression threshold');
% axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% colorbar
% 
% %==========================================================================
% remove=0.15;
% 
% bsPerM2 = bsPerM;
% lsPerM2 = lsPerM;
% FPerM2  = FPerM;
% bsPerM2(BS_thresh_frac<remove)=nan;
% lsPerM2(LS_thresh_frac<remove)=nan;
% FPerM2(F_sig_frac<remove)=nan;
% 
% figure(3) 
% fig = gcf;
% if GraphAppearance==1
%     fig.InvertHardcopy = 'off';
% elseif GraphAppearance==0
%     fig.InvertHardcopy = 'on';
% end
% clf
% textSize=12;
% 
% subplot(131) % Average period
% surf(jj_arr,ii_arr,bsPerM2);
% title({'BS period mean (hours)' ''});
% xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% xlabel('Time delay (mins)');    ylabel('Repression threshold');
% axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% colorbar
% 
% %--------------------------------------------------------------------------
% 
% subplot(132) % Average period
% surf(jj_arr,ii_arr,lsPerM2);
% title({'LS Period mean (hours)' ''});
% xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% xlabel('Time delay (mins)');    ylabel('Repression threshold');
% axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% colorbar
% 
% %--------------------------------------------------------------------------
% subplot(133) % Average period
% surf(jj_arr,ii_arr,FPerM2);
% title({'Fisher G period' ''});
% xlim([jj_arr(1) jj_arr(end)]); ylim([ii_arr(1) ii_arr(end)]);
% xlabel('Time delay (mins)');    ylabel('Repression threshold');
% axis square; view(0,90); set(gca,'fontsize',textSize); colormap(gca,parula);
% colorbar
% 
% 
% 
% if eps==0
%     figure(4)
%     fig = gcf;
%     if GraphAppearance==1
%     fig.InvertHardcopy = 'off';
%     elseif GraphAppearance==0
%     fig.InvertHardcopy = 'on';
%     end
%     clf
%     subplot(131)
%     h1=histogram(BS_thresh_frac);
%     set(h1,'edgecolor','w','facecolor',[0.15 1 0.5])
%     title('Bootstrap')
%     ylabel('Counts')
%     xlabel('Fraction passing significance')
%     subplot(132)
%     h2=histogram(LS_thresh_frac);
%     set(h2,'edgecolor','w','facecolor',[0.15 1 0.5])
%     title('Lomb Scargle')
%     ylabel('Counts')
%     xlabel('Fraction passing significance')
%     subplot(133)
%     h3=histogram(F_sig_frac);
%     set(h3,'edgecolor','w','facecolor',[0.15 1 0.5])
%     title('Fisher G')
%     ylabel('Counts')
%     xlabel('Fraction passing significance')
% end
% 
% 

end

function R=rndrng(m,n,low,high)
% rng(1);
R=(high-low)*rand(m,n)+low;
end

function [PXX,F,Thresh,RandMat]=GetBootstrapPS(mat)
    s=size(mat); 
    RandMat=[];
    BootMat=[];
    for i=1:s(2)
        vect=mat(:,i);%-mean(mat(:,i));
        kidx=randperm(numel(vect)); 
        vect_rand=vect(kidx);
        RandMat=[RandMat,vect_rand(:)];
    end
    % run periodogram on randomized data
    [PXX,F]=periodogram(RandMat,[],[],1);
    for i=1:size(PXX,2)
        M(i)=prctile(PXX(:,i),95);
    end
    Thresh=max(M);
end

function [fisher_g,pval,idx]=GetFisherG(Pxx)
    idx=find(Pxx==max(Pxx),1,'first');
    fisher_g=Pxx(idx)/sum(Pxx);
    N = length(Pxx);
        upper  = floor(1/fisher_g);
        for nn = 1:3
            I(nn) = (-1)^(nn-1)*nchoosek(N,nn)*(1-nn*fisher_g)^(N-1);
        end
    pval = sum(I);
end

function [y,pos]=microClustDetect(vect,thresh)
    % vect
    p=polyfit([1:numel(vect)]',vect,3); 
    f=polyval(p,1:numel(vect));
    idx=find(f==max(f),1,'first');
    if idx>1 && idx<numel(vect)
        % position at center
        kymo=vect;
        cent=kymo(idx);
        % fold around the max
        v1=kymo(1:idx);
        % flip around center
        v1=v1(end:-1:1);
        v2=kymo(idx:end);
        % fit a function to the average and estimate the band size
        n1=numel(v1);
        n2=numel(v2);
        if n1>n2
            vmean(1:n2)=(v1(1:n2)+v2(1:n2))/2;
            vmean(n2+1:n1)=v1(n2+1:n1);
        else
            vmean(1:n1)=(v1(1:n1)+v2(1:n1))/2;
            vmean(n1+1:n2)=v2(n1+1:n2);
        end
        lidx=find(vmean<thresh*max(vmean),1,'first');
        % size is double lidx since pattern is observed from center
        y=2*lidx-1;
        if ~isempty(lidx)& y>=2
            y=min(2*lidx-1,numel(vect));
            pos=idx;
        else
            y=NaN;
            pos=NaN;
        end
    else % if center of pattern found at edge discard
        y=NaN;
        pos=NaN;
    end
end


        

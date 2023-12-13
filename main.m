clc;
clear;
close all;

% statistics: just a timer to check the simulation time
fprintf('Script start time: %s [HH:MM:SS]\n',datetime("now"));
tic

saveData=1; % save data for the further processing (e.g. figures.m)
format default

fc = 2e9;%carrier frequency
c = physconst('LightSpeed'); %light speed
lambda=c/fc; %wavelength
dy = lambda/2; % Spacing between elements on each row (m) --> (subarray)
dz = lambda/2; % Spacing between elements on each column (m) --> (subarray)
Dy= lambda; % Spacing between elements on each row (m) --> (array)
Dz= lambda; % Spacing between elements on each column (m) --> (array)
aGranularity=1; % degree granularity of the beam pattern
angs=-90:aGranularity:90;
ang = [0;0]; %steering angles 
antenna = phased.NRAntennaElement('FrequencyRange',[1.97e9 2.3e9],'Beamwidth',[80,80]);

%% Step1: Beam pattern calculation

% Array creation (fix distance between array elements)
nR=[2 4]; % nR x nR is the number of rad. elements inside the single platform 
nE=1024; %total number of rad. elements of the entire system
nE_min=100; % minimum number of rad. elements
nPmax=6; % max number of simulated points
dP=[1.5/lambda 3/lambda]; %distance between the platform
nG=3; % number of simulated geometries

BW_mat=zeros(nG,length(nR),length(dP),nPmax); % beamwidth matrix
D_mat=zeros(nG,length(nR),length(dP),nPmax); % directivity matrix
Pat_mat=zeros(nG,length(nR),length(dP),nPmax,length(angs)); % beampattern matrix

for a=1:nG % for each geometry
    
    for k=1:length(nR) % for each configuration of nR

        % Creation of the single platform array (sub-array)
        sura = phased.URA([nR(k) nR(k)],[dz dy],'Element',antenna);
        N = ceil(sqrt(nE/nR(k)^2));   % Number of elements on each row
        M = N;   % Number of elements on each column
        nP(k,:)=floor(linspace(ceil(sqrt(nE_min/nR(k)^2)),N,nPmax)); % number of platforms
        ne=((nP(k,:).^2)*(nR(k)^2)); % number of radiating elements
        
        for i=1:length(dP) % for each configuration of dP

            for j=1:length(nP) % for each configuration of nP

                switch a % geometry selection

                    case 1 % c-URA
                        array=phased.URA('Size',[sqrt(ne(j)),sqrt(ne(j))],'ElementSpacing',[dy,dz],'Element',antenna);                    
                    case 2 % d-URA
                        array = phased.ReplicatedSubarray('Subarray', sura, ...
                                           'Layout','Rectangular',...
                                           'GridSize',[nP(k,j) nP(k,j)],'GridSpacing',[dP(i)*Dy dP(i)*Dz]);
                    case 3 % d-ELSA
                        numP=nP(k,j)^2;
                        [y,z]=elsaGeometry(numP,dP(i)*Dy,1);
                        x=y*0;
                        sPos=[x;y;z];
                        sNor=zeros(2,numP);
                        array = phased.ReplicatedSubarray('Subarray', sura, ...
                                           'Layout','Custom',...
                                           'SubarrayPosition',sPos,...
                                           'SubarrayNormal',sNor);
                end
                
                % BW, Beam pattern, Maximum directivity calculation

                BW_mat(a,k,i,j)= beamwidth(array,fc,'Cut','Elevation','dBDown',3); % beamwith computing
%                 D_mat(a,k,i,j)= directivity(array,fc,ang); %Directivity computing
%                 fprintf('D calculation step [%d,%d,%d,%d] - elapsed time: %s [HH:MM:SS:FFF]\n',a,k,i,j,datestr(seconds(toc),'HH:MM:SS:FFF'));
                p=pattern(array,fc,angs,0);
                Pat_mat(a,k,i,j,:)=p;
                D_mat(a,k,i,j)=max(p);
                fprintf('Step1 [%d/%d,%d/%d,%d/%d,%d/%d] - elapsed time: %s [HH:MM:SS] - Core part\n',a,nG,k,length(nR),i,length(dP),j,length(nP),duration(0,0,toc, 'Format', 'hh:mm:ss'));
            end
        end
    end
end
%% Step2: Single Beam Channel Capacity Without Interference

pE=0.35; % watt per single radiating element

SNR=zeros(nG,length(nR),length(dP),length(nP)); % Carrier-to-Noise ratio matrix
PRX=zeros(nG,length(nR),length(dP),length(nP)); % Received power matrix
EIRP=zeros(nG,length(nR),length(dP),length(nP)); % EIRP matrix
r_mat=zeros(nG,length(nR),length(dP),length(nP));
Ne=zeros(nG,length(nR),length(dP),length(nP)); 

for a=1:nG % for each geometry
    
    for k=1:length(nR)
        
        for i=1:length(dP)
            
            ns=nR(k);
            ne=((nP(k,:).^2)*(nR(k)^2))'; % total number of elements
            PTw=pE*ne; % PT transmitting antenna power (W)
            PT=10*log10(PTw); % PT transmitting antenna power (dBW)
            D=squeeze(D_mat(a,k,i,:)); % directivity
            BW=squeeze(BW_mat(a,k,i,:)); % beam width
            Ne(a,k,i,:)=ne; % total number of radiating elements
            GT=D; % tx gain (assumed equal to the directivity)
            EIRP(a,k,i,:)=PT+GT; % EIRP (dBW) is the effective isotopic radiated power of the transmitting antenna
            Gr=0; % Gr is the rx gain in dBi
            To=290; % To is the ambient temperature in K
            Ta=290; % Ta is the antenna temperature in K
            NF=9; % NF represents the noise figure in dB
            GrT=Gr-NF-10*log10(To+(Ta-To)*10^(-0.1*NF)); % Gr/T (dBi/K) is the figure of merit at the receiver
            RE=physconst('EarthRadius'); % radius of Earth RE         
            alpha=90; % satellite elevation angle α       
            ho=500e3; % satellite altitude ho
            d=sqrt((RE^2)*(sin(deg2rad(alpha)))^2+(ho.^2)+(2*ho.*RE))-RE*sin(deg2rad(alpha)); % slant range d            
            Ploss=20*log10(4*pi/c)+20*log10(fc)+20*log10(d); % free space path loss
            Aloss=0.5;  % Aloss (dB) atmospheric looses due to gases, rain fades etc.  
            Pmargins=1.50; % Shadow fading margin
            Pad=1; % additional loss, for example degradation due to feeder links in case of non-regenerative systems
            B=30e6; % communication bandwidth 30 MHz (S-band)
            K=-228.6; % K is the Boltzman constant (dBW/K/Hz) 
            % Average C/I within a satellite beam
            PRX(a,k,i,:)=EIRP(a,k,i,:)-Ploss-Aloss-Pmargins-Pad;           
            SNR(a,k,i,:)=EIRP(a,k,i,:)+GrT-K-Ploss-Aloss-Pmargins-Pad-10*log10(B);
            r=(ho.*(tand(BW/2)))/1e3;
            r_mat(a,k,i,:)=r;  
        end
    end
end

spectralEfficiency1=log2(1+(10.^(SNR/10)));
mbps1_mat=B*spectralEfficiency1/1024/1024;
A_km2=pi*(r_mat.^2);
thpDensity1_mat=mbps1_mat./A_km2;

fprintf('Step2 - elapsed time: %s [HH:MM:SS] - SNR\n',duration(0,0,toc, 'Format', 'hh:mm:ss'));

%% Step3: Single Beam Channel Capacity With Interference

SINR_lim=-5;
SINR_db_3=zeros(length(nR),length(dP),length(nP));
iB_mat=zeros(length(nR),length(dP),length(nP));
angIntStop=25;
for k=1:length(nR)
    for i=1:length(dP)
        for j=1:length(nP)
            iB=0;
            while (SINR_db_3(k,i,j)>SINR_lim)
                iB=iB+1;
                [D_max,angSigId]=max(squeeze(Pat_mat(3,k,i,j,:)));
                prxS=PRX(3,k,i,j);
                prXS_lin=(10^(prxS/10))/(iB+1);            
                
                angIntStart=angs(angSigId)+BW_mat(3,k,i,j);
                
                angIntIdxs=find((angs>=angIntStart) & (angs<=angIntStop));
                prxIntRange=PRX(3,k,i,j)-D_max+squeeze(Pat_mat(3,k,i,j,angIntIdxs));
                prxIntRangeAvg_lin=mean(10.^(prxIntRange/10));
                
                SIR_lin_3=prXS_lin/prxIntRangeAvg_lin;
                SIR_db_3=10*log10(SIR_lin_3);
                
                SNR_lin=(10^(SNR(3,k,i,j)/10))/(iB+1);
                SNR_db=10*log10(SNR_lin);
    
                SINR_db_3(k,i,j)=-10*log10(10.^(-0.1*SNR_db)+10^(-0.1*SIR_db_3));
                fprintf("Alternative 3 [%d, %d, %d, %d]: SIR_db %.1f SINR_db %.1f beams=%d\n",k,i,j,Ne(3,k,i,j),SIR_db_3,SINR_db_3(k,i,j),iB);
                iB_mat(k,i,j)=iB;
            end
        end
    end
end
SINR=zeros(nG,length(nR),length(dP),length(nP));
SINR(3,:,:,:)=SINR_db_3;
fprintf('Step3 - elapsed time: %s [HH:MM:SS] - SINR\n',duration(0,0,toc, 'Format', 'hh:mm:ss'));

spectralEfficiency2=log2(1+(10.^(SINR/10)));
mbps2_mat=B*spectralEfficiency2/1024/1024;
thpDensity2_mat=mbps2_mat./A_km2;

%% Save results

fprintf('Script end time: %s [HH:MM:SS]\n',datetime("now"));
filename = sprintf('data\\results_%s', datetime('now','Format','yyyyMMdd_HHmm'));
if saveData == 1
    save(filename);
end
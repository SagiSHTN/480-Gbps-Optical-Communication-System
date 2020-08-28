%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Communication System based on DP-16QAM modulation and      %
% optic fiber with Chromatic Dispersion, PMD, Random axes    %
% axes rotations and Cross talk. The system work at 480 Gbps %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%% Setup

M=16;
numofbits=2e6;%number of bits per polarization
SR=60*1e9;%symbol rate (baud [symbol/sec])
gammab_dB=1:16;
gammab=10.^(gammab_dB./10);
OSNR=(4*gammab*SR)/12.5e9;
OSNR_dB=10*log10(OSNR);
N0=zeros(1,10);
L=numofbits/log2(M);%Number of symbol per pol
BER=zeros(1,length(gammab));
SER=zeros(1,length(gammab));
sps=8;
span=51;
c=2.997e8;
lambda=1550*1e-9;%wavelength [m]
D=17*1e-6;%Dispersion [psec/nm-km]
b2=-(lambda).^2/(2*pi*c)*D;% [psec^2/km] cd coefficient
Lfiber=1000e3; %length of optic fiber [m]
Ts=1/SR;
dt=Ts/sps;%samples interval
t=0:dt:(L*sps-1)*dt;
df=1/max(t);
Fmax=1/dt;
f=-Fmax/2:df:Fmax/2;
fc=0e3;%Dynamic rotations frequency[Hz]
conv_time=0;

setup=[1,1,0]; %[1,0,0]=RDE mode, [1,1,0]=CMA+RDE mode, [0,0,1]=LMS mode, [1,1,1]=CMA+RDE and LMS mode (sepereated systems)
%% Tx

%Tx pol-X
Npilot=32; %number of pilot symbols in each period
pilot_period=1024; %data segment length
[y_x,h,Symbols_x,Bit1_x,Bit2_x,Bit3_x,Bit4_x] = QAM16_Tx(L,span,sps,Npilot,pilot_period);%Transmitter,generate random 16QAM symbols with pilots
%output is an analog signal
%Tx pol-Y
[y_y,h,Symbols_y,Bit1_y,Bit2_y,Bit3_y,Bit4_y] = QAM16_Tx(L,span,sps,Npilot,pilot_period);
Symbols=[Symbols_x;Symbols_y];
y = [y_x ; y_y]; %Vector of transmitted data composed of 2 polarizations (x and y)

%Dynamic rotations based on jones matrix
y_xout=y_x.*cos(2*pi*fc*t) + y_y.*sin(2*pi*fc*t);
y_yout=y_y.*cos(2*pi*fc*t) - y_x.*sin(2*pi*fc*t);
y_new=[y_xout;y_yout];
%% Channel-Optic fiber

DGD=0e-12; %Average DGD along the fiber
Nsection=1;% Slice the fiber into Nsection sections
%Model the optic fiber wih Nsection sections and apply 2 random rotations
%and DGD
for i=1:Nsection
    theta1=2*pi*rand(1); %Random angle [0,2pi]
    theta2=2*pi*rand(1);
    [y_t,~] = QAM16_optic_channel_pmd(y_new,f,Lfiber/Nsection,D,DGD/sqrt(Nsection),lambda,theta1,theta2);%apply optic fiber phenomena
    y_new=y_t;
end

for k=1:length(gammab)
    tic
    %AWGN
    N0=2.5/gammab(k); %Noise psd Eb/gammab
    n=sqrt(N0/2).*(randn(2,length(y_t))+1j*(randn(2,length(y_t)))); %generate gaussian white noise
    
    %% Rx
    
    r_t = y_t + n; %Received signal with AWGN
    [ycdc_t] = QAM16_Rx(r_t,h,span,sps,Ts,L,b2,Lfiber-15e3);%Apply SRRC filter and CDC.Lfiber-Lcdc
    ycdc_t_temp=ycdc_t;
    
    
    if setup(1)==1||setup(2)==1
        %CMA and RDE MIMO EQ
        N_CMA=41;mew_CMA=0.25e-5;mew_CMA_sing=0.1e-5;%CMA&RDE EQ parameters
        mew_RDE=7e-5;
        %Apply CMA and RDE MIMO EQ
        [taps_CMA,taps_RDE,epsilon_CMA,epsilon_RDE,CMA_output,RDE_output,conv_time] = CMA_RDE_eq(ycdc_t,mew_CMA,mew_CMA_sing,N_CMA,mew_RDE,Symbols,setup);
        
        %Organize data
        x_RDE=RDE_output(1,:);y_RDE=RDE_output(2,:);
        epsilon_x_RDE=epsilon_RDE(1,:);epsilon_y_RDE=epsilon_RDE(2,:);
        x_CMA=CMA_output(1,:);y_CMA=CMA_output(2,:);
        epsilon_x_CMA=epsilon_CMA(1,:);epsilon_y_CMA=epsilon_CMA(2,:);
        
        %Pol and syncronization
        r_x = max(abs(xcorr(x_RDE(1:2:200e3),Symbols_x((conv_time/2)+1:(conv_time/2)+1+100e3))));
        r_y = max(abs(xcorr(x_RDE(1:2:200e3),Symbols_y((conv_time/2)+1:(conv_time/2)+1+100e3))));
        
        if r_x > r_y %Polarization match
            offset_xx=finddelay(x_RDE(1:2:end),Symbols_x((conv_time/2)+1:end));%find offset
            x_RDE=circshift(x_RDE,offset_xx*2);%Sync
            offset_yy=finddelay(y_RDE(1:2:end),Symbols_y((conv_time/2)+1:end));%find offset
            y_RDE=circshift(y_RDE,offset_yy*2);%Sync
            %Rotations Compensation
            common_phase_x=angle(mean(x_RDE(end-2048*2:2:end).*conj(Symbols_x(end-1024*2:end))));
            common_phase_y=angle(mean(y_RDE(end-2048*2:2:end).*conj(Symbols_y(end-1024*2:end))));
            x_RDE=x_RDE*exp(-1j*common_phase_x);%apply phase compensation
            y_RDE=y_RDE*exp(-1j*common_phase_y);%apply phase compensation
            RDE_output= [x_RDE(1:2:end);y_RDE(1:2:end)];
            
        else%Polarization mismatch
            offset_xy=finddelay(x_RDE(1:2:end),Symbols_y((conv_time/2)+1:end));%find offset
            x_RDE=circshift(x_RDE,offset_xy*2);%Sync
            offset_yx=finddelay(y_RDE(1:2:end),Symbols_x((conv_time/2)+1:end));%find offset
            y_RDE=circshift(y_RDE,offset_yx*2);%Sync
            %Rotations Compensation
            common_phase_x=angle(mean(x_RDE(end-2048*2:2:end).*conj(Symbols_y(end-1024*2:end))));
            common_phase_y=angle(mean(y_RDE(end-2048*2:2:end).*conj(Symbols_x(end-1024*2:end))));
            x_RDE=x_RDE*exp(-1j*common_phase_x);%apply phase compensation
            y_RDE=y_RDE*exp(-1j*common_phase_y);%apply phase compensation
            RDE_output= [y_RDE(1:2:end);x_RDE(1:2:end)];%Sample-sps=1
        end
    end
    
    
    if setup(3) == 1
        %LMS MIMO EQ
        ycdc_t=ycdc_t_temp;
        N_LMS=41;%LMS filters size [taps]
        delay=fix(N_LMS/2);
        r_x = max(abs(xcorr(ycdc_t(1,1:2:10e4),Symbols_x(1:5e4))));
        r_y = max(abs(xcorr(ycdc_t(2,1:2:10e4),Symbols_y(1:5e4))));
        
        if r_x > r_y %Polarization match
            offset_xx=finddelay(ycdc_t(1,1:2:end),Symbols_x);%find offset
            ycdc_t(1,:)=circshift(ycdc_t(1,:),offset_xx*2-delay);%Sync and include the LMS delay
            offset_yy=finddelay(ycdc_t(2,1:2:end),Symbols_y);
            ycdc_t(2,:)=circshift(ycdc_t(2,:),offset_yy*2-delay);
            
        else    %Polarization mismatch
            offset_xy=finddelay(ycdc_t(1,1:2:end),Symbols_y);%find offset
            ycdc_t(1,:)=circshift(ycdc_t(1,:),offset_xy*2-delay);%Sync and include the LMS delay
            offset_yx=finddelay(ycdc_t(2,1:2:end),Symbols_x);
            ycdc_t(2,:)=circshift(ycdc_t(2,:),offset_yx*2-delay);
            ycdc_t=[ycdc_t(2,:);ycdc_t(1,:)];
        end
        
        %LMS
        mew_LMS=4e-5;%LMS setp size
        pilot_symb=[Symbols_x(1:Npilot);Symbols_y(1:Npilot)];%LMS pilot symbols
        %Apply LMS MIMO EQ
        [taps_LMS,x_LMS,y_LMS,epsilon_x_LMS,epsilon_y_LMS] = LMS_eq(ycdc_t,mew_LMS,N_LMS,pilot_symb,Npilot,pilot_period);
        
        
        %Sync data after LMS
        offset_xx=finddelay(x_LMS(1:2:end),Symbols_x);%find offset
        x_LMS=circshift(x_LMS,offset_xx*2);%Sync
        offset_yy=finddelay(y_LMS(1:2:end),Symbols_y);
        y_LMS=circshift(y_LMS,offset_yy*2);
        LMS_output= [x_LMS(1:2:end);y_LMS(1:2:end)];%Sample-sps=1
    end
    
    %% Detection
    
    %RDE Detect
    if setup(1) == 1
        [Symbols_x_det_RDE,Bit1_x_det_RDE,Bit2_x_det_RDE,Bit3_x_det_RDE,Bit4_x_det_RDE] = symbols_detect(RDE_output(1,:),length(RDE_output)); %Detect symbols and bits of polarization X
        [Symbols_y_det_RDE,Bit1_y_det_RDE,Bit2_y_det_RDE,Bit3_y_det_RDE,Bit4_y_det_RDE] = symbols_detect(RDE_output(2,:),length(RDE_output)); %Detect symbols and bits of polarization Y
    end
    if setup(3) == 1
        %LMS Detect
        [Symbols_x_det_LMS,Bit1_x_det_LMS,Bit2_x_det_LMS,Bit3_x_det_LMS,Bit4_x_det_LMS] = symbols_detect(LMS_output(1,:),length(LMS_output)); %Detect symbols and bits of polarization X
        [Symbols_y_det_LMS,Bit1_y_det_LMS,Bit2_y_det_LMS,Bit3_y_det_LMS,Bit4_y_det_LMS] = symbols_detect(LMS_output(2,:),length(LMS_output)); %Detect symbols and bits of polarization Y
        
    end
    
    
    %% Calculations
    
    alpha_sym_RDE=(conv_time/2)+0e3;%offset symbols to convergance state
    alpha_sym_LMS=100e3;%offset symbols to convergance state
    beta=2e3;%end of the signal slicing
    e=0e3;
    if setup(1) == 1
        %RDE calcs
        %SER
        error_symb_x_RDE = sum(Symbols_x(alpha_sym_RDE+1:end-beta) ~= Symbols_x_det_RDE(1+e:end-beta));%number of symbol errors in pol-X
        error_symb_y_RDE = sum(Symbols_y(alpha_sym_RDE+1:end-beta) ~= Symbols_y_det_RDE(1+e:end-beta));%number of symbol errors in pol-Y
        SER_x_RDE(k) = error_symb_x_RDE/(length(Symbols_x)-alpha_sym_RDE-beta);%numeric SER calc pol-X
        SER_y_RDE(k) = error_symb_y_RDE/(length(Symbols_y)-alpha_sym_RDE-beta);%numeric SER calc pol-Y
        SER_RDE(k) = 0.5*(SER_x_RDE(k)+SER_y_RDE(k));%Total SER
        %BER
        error_bit_x_RDE = sum(Bit1_x(alpha_sym_RDE+1:end-beta) ~= Bit1_x_det_RDE(1+e:end-beta)) + sum(Bit2_x(alpha_sym_RDE+1:end-beta) ~= Bit2_x_det_RDE(1+e:end-beta)) + sum(Bit3_x(alpha_sym_RDE+1:end-beta) ~= Bit3_x_det_RDE(1+e:end-beta)) + sum(Bit4_x(alpha_sym_RDE+1:end-beta) ~= Bit4_x_det_RDE(1+e:end-beta));%number of bits errors in pol-X
        error_bit_y_RDE = sum(Bit1_y(alpha_sym_RDE+1:end-beta) ~= Bit1_y_det_RDE(1+e:end-beta)) + sum(Bit2_y(alpha_sym_RDE+1:end-beta)~= Bit2_y_det_RDE(1+e:end-beta)) + sum(Bit3_y(alpha_sym_RDE+1:end-beta)~= Bit3_y_det_RDE(1+e:end-beta)) + sum(Bit4_y(alpha_sym_RDE+1:end-beta) ~= Bit4_y_det_RDE(1+e:end-beta));%number of bits errors in pol-Y
        BER_x_RDE(k) = error_bit_x_RDE/(numofbits-alpha_sym_RDE-beta);%numeric BER calc pol-X
        BER_y_RDE(k) = error_bit_y_RDE/(numofbits-alpha_sym_RDE-beta);%numeric BER calc pol-Y
        BER_RDE(k) = 0.5*(BER_x_RDE(k) + BER_y_RDE(k));%Total BER
    end
    
    if setup(3) == 1
        %LMS calc
        %SER
        error_symb_x_LMS = sum(Symbols_x(alpha_sym_LMS:end) ~= Symbols_x_det_LMS(alpha_sym_LMS:end));%number of symbol errors in pol-X
        error_symb_y_LMS = sum(Symbols_y(alpha_sym_LMS:end) ~= Symbols_y_det_LMS(alpha_sym_LMS:end));%number of symbol errors in pol-Y
        SER_x_LMS(k) = error_symb_x_LMS/(length(Symbols_x)-alpha_sym_LMS);%numeric SER calc pol-X
        SER_y_LMS(k) = error_symb_y_LMS/(length(Symbols_y)-alpha_sym_LMS);%numeric SER calc pol-Y
        SER_LMS(k) = 0.5*(SER_x_LMS(k)+SER_y_LMS(k));%Total SER
        %BER
        error_bit_x_LMS = sum(Bit1_x(alpha_sym_LMS:end) ~= Bit1_x_det_LMS(alpha_sym_LMS:end)) + sum(Bit2_x(alpha_sym_LMS:end) ~= Bit2_x_det_LMS(alpha_sym_LMS:end)) + sum(Bit3_x(alpha_sym_LMS:end) ~= Bit3_x_det_LMS(alpha_sym_LMS:end)) + sum(Bit4_x(alpha_sym_LMS:end) ~= Bit4_x_det_LMS(alpha_sym_LMS:end)); %number of bits errors in pol-X
        error_bit_y_LMS = sum(Bit1_y(alpha_sym_LMS:end) ~= Bit1_y_det_LMS(alpha_sym_LMS:end)) + sum(Bit2_y(alpha_sym_LMS:end)~= Bit2_y_det_LMS(alpha_sym_LMS:end)) + sum(Bit3_y(alpha_sym_LMS:end)~= Bit3_y_det_LMS(alpha_sym_LMS:end)) + sum(Bit4_y(alpha_sym_LMS:end) ~= Bit4_y_det_LMS(alpha_sym_LMS:end));%number of bits errors in pol-Y
        BER_x_LMS(k) = error_bit_x_LMS/(numofbits-alpha_sym_LMS);%numeric BER calc pol-X
        BER_y_LMS(k) = error_bit_y_LMS/(numofbits-alpha_sym_LMS);%numeric BER calc pol-Y
        BER_LMS(k) = 0.5*(BER_x_LMS(k) + BER_y_LMS(k));%Total BER
        
    end
    
    %Theory calc
    SER_theory(k) = 1-(1-(((2*(sqrt(M)-1))/sqrt(M))*(qfunc(sqrt((3*log2(M)*gammab(k))/(M-1))))))^2;%Theortical SER calc
    BER_theory(k) = (2/log2(M))*(1-1/sqrt(M))*erfc(sqrt((3*log2(M)*gammab(k))/(2*(M-1))));% Theortical BER calc
    
    toc; disp(['gammab=',num2str(gammab_dB(k))]);
end


%% Plots

if setup(1)==1
    %RDE plots
    %BER plot
    figure
    semilogy(OSNR_dB,BER_RDE);
    title('CMA&RDE BER - DGD=0psec')
    grid on
    xlabel('OSNR [dB]')
    ylabel('BER')
    hold on
    semilogy(OSNR_dB,BER_theory);
    hold off
    legend('Numeric','Theory')
    %SER plot
    figure
    semilogy(OSNR_dB,SER_RDE);
    title('CMA&RDE SER - DGD=0psec')
    grid on
    xlabel('OSNR [dB]')
    ylabel('SER')
    hold on
    semilogy(OSNR_dB,SER_theory);
    hold off
    legend('Numeric','Theory')
    %Error function and constellation plot
    figure;
    subplot(2,2,1)
    plot(abs(epsilon_x_RDE.^2));
    title('epsilon x rde')
    subplot(2,2,3)
    plot(abs(epsilon_y_RDE.^2));
    title('epsilon y rde')
    subplot(2,2,2)
    display_constellation_16QAM(RDE_output(1,1:end))
    title('Pol X')
    subplot(2,2,4)
    display_constellation_16QAM(RDE_output(2,1:end))
    title('Pol Y')
    
end
if setup(3)==1
    %LMS Plots
    %BER plot
    figure
    semilogy(OSNR_dB,BER_LMS);
    grid on
    title('LMS BER - DGD=0psec')
    xlabel('OSNR [dB]')
    ylabel('BER')
    hold on
    semilogy(OSNR_dB,BER_theory);
    hold off
    legend('Numeric','Theory')
    %SER plot
    figure
    semilogy(OSNR_dB,SER_LMS);
    grid on
    title('LMS SER - DGD=0psec')
    xlabel('OSNR [dB]')
    ylabel('SER')
    hold on
    semilogy(OSNR_dB,SER_theory);
    hold off
    legend('Numeric','Theory')
    %Error function and constellation plot
    figure;
    subplot(2,2,1)
    plot(abs(epsilon_x_LMS.^2));
    title('epsilon x lms')
    subplot(2,2,3)
    plot(abs(epsilon_y_LMS.^2));
    title('epsilon y lms')
    subplot(2,2,2)
    display_constellation_16QAM(LMS_output(1,alpha_sym_LMS+1:end))
    title('Pol X')
    subplot(2,2,4)
    display_constellation_16QAM(LMS_output(2,alpha_sym_LMS+1:end))
    title('Pol Y')
end

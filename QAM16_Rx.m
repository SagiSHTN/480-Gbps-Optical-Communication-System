function [ycdc_t] = QAM16_Rx(r_t,h,span,sps,Ts,L,b2,Lfiber)
%Simulates a receiver of a communication system based on
%DP16-QAM Modulation, SRRC matched filter and chromatic dispersion compensation.
%Input:             r_t - received signal
%                   h - SRRC matched filter 
%                   span - matched filter width in symbols
%                   sps - number of samples per symbol 
%                   Ts - symbol time
%                   L - number of symbols per polarization
%                   b2 - chromatic dispersion coefficient
%                   Lfiber - optic fiber length [m]
%Output:            ycdc_t - received signal after SRRC filtering CD compensation, signal is in time domain 


%Matched filter - sps=8
ymf_x = conv(r_t(1,:),h); % matched filter (srrc)
ymf_y = conv(r_t(2,:),h); 
ymf_x=ymf_x(sps/2*span+1:end-sps/2*span); %Remove the dealy parts after filtering
ymf_y=ymf_y(sps/2*span+1:end-sps/2*span); 

%CDC
ymf1_x = downsample(ymf_x,4,0); %sps=2
ymf1_y = downsample(ymf_y,4,0); 
sps1=2;
dt1=Ts/sps1;
t1=0:dt1:(L*sps1-1)*dt1;
df1=1/max(t1);
Fmax1=1/dt1;
f1=-Fmax1/2:df1:Fmax1/2;
Hcdc_f = exp(-1*(1j*(b2/2)*(2*pi*f1).^2*Lfiber)); %cdc

Ymfx_f = fftshift((fft(ymf1_x)));%Frequency domain
Ymfy_f = fftshift((fft(ymf1_y)));
Ycdc_x_f = Ymfx_f.*Hcdc_f;%Apply CDC
Ycdc_y_f = Ymfy_f.*Hcdc_f;
ycdc_x_t = ifft(ifftshift(Ycdc_x_f));
ycdc_y_t = ifft(ifftshift(Ycdc_y_f));


ycdc_t = [ycdc_x_t ; ycdc_y_t];%Time domain signal
end


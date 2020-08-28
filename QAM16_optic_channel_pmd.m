function [y_t,b2] = QAM16_optic_channel_pmd(y,f,Lfiber,D,DGD,lambda,theta1,theta2)
%Simulates optic fiber phenomena which includes Chromatic Dispersion, PMD
%and random polarization rotations
%Input:    y - Transmitted signal
%          f - Frequency axis
%          Lfiber - Optic fiber length [m]
%          D - Dispersion [psec/nm-km]
%          DGD - Differential group delay [sec]
%          lambda - Laser wavelength [m]
%          theta1,theta2 - polraization rotations angle [radians]
%Output:   b2 - Chromatic dispersion coefficient [psec^2/km]
%          y_t - signal in time domain with the applied fiber phenomena

c=2.997e8;%speed of light [m/sec]
b2=-(lambda).^2/(2*pi*c)*D;% [psec^2/km] cd coefficient

%SOP 1
RotationMat1 = [cos(theta1) -sin(theta1) ; sin(theta1) cos(theta1)];%jones matrix 
yrot = RotationMat1*y; %applied 1st rotation

%CD
Yx_f = fftshift(fft(yrot(1,:)));
Yy_f = fftshift(fft(yrot(2,:)));
Hcd_f = exp(1j*(b2/2)*((2*pi*f).^2)*Lfiber); %chromatic dispersion
Yx_f1 = Yx_f(1:length(f)).*Hcd_f; %signal with Rotation&CD (freq domain)
Yy_f1 = Yy_f(1:length(f)).*Hcd_f;

%DGD
Yx_f3=Yx_f1.*exp((-1j*2*pi.*f*DGD)/2); %Apply DGD (time domain)
Yy_f3=Yy_f1.*exp((+1j*2*pi.*f*DGD)/2);
Y_f=[Yx_f3;Yy_f3];
Yx_f_dgd = Y_f(1,:);
Yy_f_dgd = Y_f(2,:);
yx_t = ifft(ifftshift(Yx_f_dgd)); %signal with Rotation&CD%DGD (time domain)
yy_t = ifft(ifftshift(Yy_f_dgd));
y_t1 = [yx_t;yy_t];

%SOP 2
RotationMat2 = [cos(theta2) -sin(theta2) ; sin(theta2) cos(theta2)];
yrot2 = RotationMat2*y_t1; %applied 2nd rotation (after CD&DGD)
y_t = yrot2; %signal in time after rotation-cd-dgd-rotation

end


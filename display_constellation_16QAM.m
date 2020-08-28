function display_constellation_16QAM(IQ)

%for display of QPSK constellation diagram with color density
%Ix=real part of teh symbol vector
%Qx=imaginary aprt of the symbol vector

Ix=real(IQ);
Qx=imag(IQ);
data=[Ix',Qx'];
count=hist2d(data,-5.5:0.05:5.5,-5.5:0.05:5.5);
imagesc(-5.5:0.05:5.5,-5.5:0.05:5.5,count);
colormap([1,1,1;jet(128)]) 
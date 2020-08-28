function [Symbols_det,Bit1_det,Bit2_det,Bit3_det,Bit4_det] = symbols_detect(r,L)
%Detect the input symbols based on 16QAM modulation and grey coding
%The detection rule is Maximum Likelihood decision rule
%Input:     r - input signal
%           L - number of symbols to detect
%Output:    Symbols_det - symbols after detection
%           Bit1_det,Bit2_det,Bit3_det,Bit4_det  - bits after detection

Symbols_det=zeros(1,L);
Bit1_det=zeros(1,L);
Bit2_det=zeros(1,L);
Bit3_det=zeros(1,L);
Bit4_det=zeros(1,L);

Symbol_ind = real(r)>0;
Symbols_det(Symbol_ind) = (real(r(Symbol_ind))>2)*2 + 1;
Bit1_det(Symbol_ind) = 0;
Bit2_det(Symbol_ind) = real(r(Symbol_ind))>2;
Symbol_ind = real(r)<0;
Symbols_det(Symbol_ind) = (real(r(Symbol_ind))<-2)*(-2) - 1;
Bit1_det(Symbol_ind) = 1;
Bit2_det(Symbol_ind) = real(r(Symbol_ind))<-2;


Symbol_ind = imag(r)>0;
Symbols_det(Symbol_ind) = (imag(r(Symbol_ind))>2)*2j + 1j + Symbols_det(Symbol_ind);%Imaginary number is 1j or 3j
Bit3_det(Symbol_ind) = 0;
Bit4_det(Symbol_ind) = imag(r(Symbol_ind))>2;
Symbol_ind = imag(r)<0;
Symbols_det(Symbol_ind) = (imag(r(Symbol_ind))<-2)*(-2j) - 1j + Symbols_det(Symbol_ind);
Bit3_det(Symbol_ind) = 1;
Bit4_det(Symbol_ind) = imag(r(Symbol_ind))<-2;

end


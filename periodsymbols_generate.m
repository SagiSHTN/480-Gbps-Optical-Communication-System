function [Symbols,Bit1,Bit2,Bit3,Bit4] = periodsymbols_generate(L,Npilot,pilot_period)
%The function generate L random 16QAM modulated symbols 
%The data is made of segments by length of pilot_period and Npilot number 
%of pilot symbols in each data segment.
%Input:     L - number of symbols
%           Npilot - number of pilot symbols per period
%           pilot_period - data segment size in symbols
%
%Output:    Symbols - output 16QAM symbols with L/pilot_period segments of
%data
%           Bit1,Bit2,Bit3,Bit4 - random generated bits


Bit1_pilot=rand(1,Npilot)>0.5;
Bit2_pilot=rand(1,Npilot)>0.5;
Bit3_pilot=rand(1,Npilot)>0.5;
Bit4_pilot=rand(1,Npilot)>0.5;

Bit1 =  rand(1,L) > 0.5; %generate bits (0,1)[uniformly distributed]
Bit2 =  rand(1,L) > 0.5;
Bit3 =  rand(1,L) > 0.5;
Bit4 =  rand(1,L) > 0.5;
count=0;
for j=1:L
    if (rem(j,pilot_period)<=Npilot ) && (rem(j,pilot_period)>0)
        count=count+1;
        if count>Npilot
            count=1;
        end
        Bit1(j)=Bit1_pilot(count);
        Bit2(j)=Bit2_pilot(count);
        Bit3(j)=Bit3_pilot(count);
        Bit4(j)=Bit4_pilot(count);

    end

end


%Symbols Mapping
I=zeros(1,L); %Inphase
Q=zeros(1,L); %Quadrature

%Inphase 1 or 3
ind=(Bit1==0);
I(ind)=1;
I(ind)=((Bit2(ind)==1)*2)+1;

%Inphase -1 or -3
ind=(Bit1==1);
I(ind)=-1;
I(ind)=((Bit2(ind)==1)*(-2))-1;

%Quadrature 1j or 3j
ind=(Bit3==0);
%     I(ind)=1j;
I(ind)= I(ind)+((Bit4(ind)==1)*2j)+1j;

%Quadrature -1j or -3j
ind=(Bit3==1);
%     I(ind)=1j;
I(ind)=I(ind)+((Bit4(ind)==1)*(-2j))-1j;

Symbols=I+1j*Q;

end



function [y,h,Symbols,Bit1,Bit2,Bit3,Bit4] = QAM16_Tx(L,span,sps,Npilot,pilot_period)
%Simulates a transmitter of a communication system based on
%16-QAM Modulation with grey coding and SRRC pulse shaper.
%The function generate random data and pilots and insert Npilot pilot symbols
%every pilot_period symbols.
%Inputs:         L - number of total symbols generated
%                span - pulse shaper width in symbols
%                sps - number of samples per symbol
%                Npilot - number of pilot symbols per period
%                pilot_period - total number of symbols per period
%                each period contains Npilot pilot symbols +
%                pilot_period-Npilot data
%Output:         y - transmitted signal in time domain
%                h - SRRC filter 
%                Symbols - 16QAM symbols with L/pilot_period segments of
%                data
%                Bit1,Bit2,Bit3,Bit4 - random generated bits

[Symbols,Bit1,Bit2,Bit3,Bit4]=periodsymbols_generate(L,Npilot,pilot_period); %Generate data with pilots
Symbols_up=upsample(Symbols,sps); %upsample Symbols by sps factor
h=rcosdesign(0.1,span,sps,'sqrt'); %Pulse shaper - square root raised cosine
y = conv(Symbols_up,h); % filter the symbols with SRRC
y=y(sps/2*span+1:end-sps/2*span); %Remove the dealy parts after filtering

end


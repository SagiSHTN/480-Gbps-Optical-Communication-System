function [taps,x_out,y_out,epsilon_x,epsilon_y] = LMS_eq(y_mimo,mew,N,pilot_symb,Npilot,pilot_period)
%The function apply LMS MIMO equalization on the received signal
%Input:     y_mimo - input signal
%           mew - step size of LMS algorithm
%           N - LMS h filters size [taps]
%           pilot_symb - LMS pilot symbols
%           Npilot - number of pilot symbols
%           pilot_period - number of symbols in each data segment
%Output:    taps - last 4 filters taps
%           x_out,y_out - signal in pol X and pol Y after LMS equalization
%           epsilon_x,epsilon_y - error function of pol X and pol Y

x_in=y_mimo(1,:);%sps=2
y_in=y_mimo(2,:);
mew1=mew;

%initialize filters
if mod(N,2)==0
    N=N+1;%filters size
end
h(1).xx=zeros(1,N).';%column vector
h(1).xx(ceil(N/2))=1;
h(1).yy=zeros(1,N).';
h(1).yy(ceil(N/2))=1;
h(1).xy=zeros(1,N).';
h(1).yx=zeros(1,N).';

count=0;
pilot_symb_count=0;
for k=1:length(x_in)
    if k<=N
        X_IN = flip([zeros(1,N-k),x_in(1:k)]).';%pol X input buffer,N elements sliding block
        Y_IN = flip([zeros(1,N-k),y_in(1:k)]).';%pol X input buffer,N elements sliding block
    else
        X_IN = flip(x_in(k-N+1:k)).';%N elements sliding block
        Y_IN = flip(y_in(k-N+1:k)).';
    end
    
    x_out(k) = h(k).xx.'*X_IN + h(k).xy.'*Y_IN ;%Apply equalization with h filters based on MIMO butterfly setup
    y_out(k) = h(k).yx.'*X_IN + h(k).yy.'*Y_IN ;
    
    
    if  mod(k,2)>0%sps=1
        count=count+1;
        
        if (rem(count,pilot_period)<=Npilot ) && (rem(count,pilot_period)>0)%pilot symbols stage
            
            pilot_symb_count=pilot_symb_count+1;%what pilot symbol is it
            p=10;
            mew=mew1;
            mew=p*mew;%bigger set size in pilot symblos stage
            Detected_x=symbols_detect(pilot_symb(1,pilot_symb_count),1);
            Detected_y=symbols_detect(pilot_symb(2,pilot_symb_count),1);
        else
            pilot_symb_count=0;
            p=1;
            mew=mew1;
            mew=p*mew;
            [Detected_x,~,~,~,~] = symbols_detect(x_out(k),1);%detect which symbol is it (ML rule)
            [Detected_y,~,~,~,~] = symbols_detect(y_out(k),1);
        end
        epsilon_x(k) = Detected_x - x_out(k);%error function update based on the error between the detected symbol and actual symbol
        epsilon_y(k) = Detected_y - y_out(k);
        
        %Calculate next stage taps
        h(k+1).xx = h(k).xx + mew*epsilon_x(k)*conj(X_IN);
        h(k+1).xy = h(k).xy + mew*epsilon_x(k)*conj(Y_IN);
        h(k+1).yx = h(k).yx + mew*epsilon_y(k)*conj(X_IN);
        h(k+1).yy = h(k).yy + mew*epsilon_y(k)*conj(Y_IN);
    else
        epsilon_x(k) = epsilon_x(k-1);
        epsilon_y(k) = epsilon_y(k-1);
        
        %Calculate next stage taps
        h(k+1).xx = h(k).xx;
        h(k+1).xy = h(k).xy;
        h(k+1).yx = h(k).yx;
        h(k+1).yy = h(k).yy;
    end
    
end


h(1).xx = h(k).xx;
h(1).xy = h(k).xy;
h(1).yy = h(k).yy;
h(1).yx = h(k).yx;


taps=h(k);
end


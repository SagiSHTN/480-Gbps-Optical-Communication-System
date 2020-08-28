function [taps,x_out,y_out,epsilon_x,epsilon_y,Singularity_RDE] = RDE_eq_sps2(ycdc_t,mew,N,taps_CMA)
%The function apply CMA MIMO equalization on the received signal
%Input:     ycdc_t - input signal
%           mew - step size of RDE algorithm
%           N - RDE h filters size [taps]
%           taps_CMA - taps used for initialization
%Output:    taps - last 4 filters taps
%           x_out,y_out - signal in pol X and pol Y after CMA equalization
%           epsilon_x,epsilon_y - error function of pol X and pol Y
%           Signgularity_RDE - flag used to warn of Singularity in the
%           algorithm

x_in=ycdc_t(1,:);%sps=2
y_in=ycdc_t(2,:);


%initialize filters
if mod(N,2)==0
    N=N+1;%filters size
end
if nargin==4
    h(1).xx=taps_CMA.xx;%column vector
    h(1).yy=taps_CMA.yy;
    h(1).xy=taps_CMA.xy;
    h(1).yx=taps_CMA.yx;
else
    h(1).xx=zeros(1,N).';%column vector
    h(1).xx(ceil(N/2))=1;
    h(1).yy=zeros(1,N).';
    h(1).yy(ceil(N/2))=1;
    h(1).xy=zeros(1,N).';
    h(1).yx=zeros(1,N).';
end
for iter=1:1
    for k=1:length(x_in)
        if k<=N
            X_IN = flip([zeros(1,N-k),x_in(1:k)]).';%pol X input buffer,N elements sliding block
            Y_IN = flip([zeros(1,N-k),y_in(1:k)]).';%pol Y input buffer,N elements sliding block
        else
            X_IN = flip([zeros(1,N-k),x_in(k-N+1:k)]).';%N elements sliding block
            Y_IN = flip([zeros(1,N-k),y_in(k-N+1:k)]).';
        end
        
        x_out(k) = h(k).xx.'*X_IN + h(k).xy.'*Y_IN ;%Apply equalization with h filters based on MIMO butterfly setup
        y_out(k) = h(k).yx.'*X_IN + h(k).yy.'*Y_IN ;
        
        R1=sqrt(2); R2=sqrt(10) ; R3=sqrt(18); %16-QAM radii
        error_x=[(R1-(abs(x_out(k)))), (R2-(abs(x_out(k)))),(R3-(abs(x_out(k))))];%calculate current symbol distance from each constellation radius
        error_y=[(R1-(abs(y_out(k)))), (R2-(abs(y_out(k)))),(R3-(abs(y_out(k))))];
        [~,indx]=min(abs(error_x));%find which radius is the closest to the checked symbol
        [~,indy]=min(abs(error_y));
        
        if  mod(k,2)>0%sps=1
            
            epsilon_x(k)=error_x(indx(1));%update error function
            epsilon_y(k)=error_y(indy(1));
            
            %Calculate next stage taps
            h(k+1).xx = h(k).xx + mew*epsilon_x(k)*x_out(k)*conj(X_IN);
            h(k+1).xy = h(k).xy + mew*epsilon_x(k)*x_out(k)*conj(Y_IN);
            h(k+1).yx = h(k).yx + mew*epsilon_y(k)*y_out(k)*conj(X_IN);
            h(k+1).yy = h(k).yy + mew*epsilon_y(k)*y_out(k)*conj(Y_IN);
        else
            epsilon_x(k) =epsilon_x(k-1);
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
end  
%Singularity check
alpha=200e3;%100e3 symbols in sps2
r_xy = max(abs(xcorr(x_out(1:end-alpha),y_out(1:end-alpha))));%Cross correlation
r_xx = max(abs(xcorr(x_out(1:end-alpha),x_out(1:end-alpha))));
r_yy = max(abs(xcorr(y_out(1:end-alpha),y_out(1:end-alpha))));
% RDE_sing=r_xy/sqrt(r_xx*r_yy)
if r_xy/sqrt(r_xx*r_yy)>0.8 %There is singularity, update taps
    Singularity_RDE=1;
    disp('There is Singulairty in RDE');
    for k=1:N%update taps incase of singulairty
        h(length(x_in)).xx(k) = 0.5*(h(length(x_in)).xx(k) + conj(h(length(x_in)).yy(N-k+1)));
        h(length(x_in)).xy(k) = 0.5*(h(length(x_in)).xy(k) - conj(h(length(x_in)).yx(N-k+1)));
    end
    for k=1:N
        h(length(x_in)).yy(k) = conj(h(length(x_in)).xx(N-k+1));
        h(length(x_in)).yx(k) = -conj(h(length(x_in)).xy(N-k+1));
    end
else
    Singularity_RDE=0;%no singulairty
end


taps=h(k);
end




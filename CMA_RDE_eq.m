function [taps_CMA,taps_RDE,epsilon_CMA,epsilon_RDE,CMA_output,RDE_output,conv_time] = CMA_RDE_eq(ycdc_t,mew_CMA,mew_CMA_sing,N_CMA,mew_RDE,Symbols,setup)
%The function apply CMA and RDE MIMO equalization based on the stated algorithms
%Input:     ycdc_t - received signal
%           mew_CMA - step size of CMA
%           mew_CMA_sing - step size of CMA incase of Singularity
%           N_CMA - CMA&RDE filter size,number of taps
%           mew_RDE - step size of CMA
%           Symbols - transmitted symbols -to calculate optimal radius for
%           CMA
%           setup - system mode
%Output:    taps_CMA,taps_RDE - last 4 filter taps of each MIMO EQ
%           epsilon_CMA,epsilon_RDE - error function of each MIMO EQ
%           CMA_output,RDE_output - output signal from each MIMO EQ
%           conv_time - number of symbols used in Pre-Convergence state
%           [symbols]

ycdc_t_temp=ycdc_t;
N_RDE=N_CMA;
%% CMA
%Pre-Convergance state

if setup(2) == 1
    conv_time=600e3;%How many symbols assigned to CMA
    ycdc_t_new=[ycdc_t(1,1:conv_time);ycdc_t(2,1:conv_time)];%input signal of CMA EQ
    [taps_CMA,x_CMA,y_CMA,epsilon_x_CMA,epsilon_y_CMA,Singularity_CMA] = CMA_eq(ycdc_t_new,mew_CMA,N_CMA,Symbols);%Apply CMA MIMO EQ
    CMA_output=[x_CMA;y_CMA];
    if Singularity_CMA==1 %reapeat CMA MIMO EQ if singulairty happened
        [taps_CMA,x_CMA,y_CMA,epsilon_x_CMA,epsilon_y_CMA,Singularity_CMA] = CMA_eq(ycdc_t_new,mew_CMA_sing,N_CMA,Symbols,taps_CMA);
        CMA_output=[x_CMA;y_CMA];
    end
end
%% RDE

if setup(1) == 1
    if setup(2) == 0%System made of RDE only, Setup=[1,0,0]
        conv_time=0;
        epsilon_x_CMA=0;epsilon_y_CMA=0;x_CMA=0;y_CMA=0;taps_CMA=0;
        ycdc_t_temp=[ycdc_t_temp(1,conv_time+1:end);ycdc_t_temp(2,conv_time+1:end)];%RDE input signal start after CMA pre-convergence state
        [taps_RDE,x_RDE,y_RDE,epsilon_x_RDE,epsilon_y_RDE,Singularity_RDE] = RDE_eq_sps2(ycdc_t_temp,mew_RDE,N_RDE);%Apply RDE MIMO EQ
        RDE_output= [x_RDE;y_RDE];
    else%System is made of CMA and RDE EQs
        ycdc_t_temp=[ycdc_t(1,conv_time+1:end);ycdc_t(2,conv_time+1:end)];
        [taps_RDE,x_RDE,y_RDE,epsilon_x_RDE,epsilon_y_RDE,Singularity_RDE] = RDE_eq_sps2(ycdc_t_temp,mew_RDE,N_RDE,taps_CMA);
        RDE_output= [x_RDE;y_RDE];
    end
    if Singularity_RDE==1
        [taps_RDE,x_RDE,y_RDE,epsilon_x_RDE,epsilon_y_RDE,Singularity_RDE] = RDE_eq_sps2(ycdc_t_temp,mew_RDE,N_RDE,taps_RDE);
        RDE_output= [x_RDE;y_RDE];
        
    end
end
%% outputs
%organize outputs
epsilon_CMA=[epsilon_x_CMA;epsilon_y_CMA];
epsilon_RDE=[epsilon_x_RDE;epsilon_y_RDE];
RDE_output=[x_RDE;y_RDE];
CMA_output=[x_CMA;y_CMA];


end

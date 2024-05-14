close all; clear all; clc;

Fs=273*10^6;%ADC sample rate
fc=1*10^9;%True carrier frequency
ps=1*10^6;%True symbol rate

B_e=(1.6+(rand(1,1)-0.5)*0.6)*ps;
fc_e=fc*(1+(rand(1,1)-0.5)*0.0001);



aa=readmatrix('data/QAM16_1M_1G_4_0.2_100_num10.csv');%The actual values collected by ADC
data=aa(:,4)';
data=data-2048;%Data unsigned to signed 
data=data*1/max(abs(data));%Amplitude normalization

figure;
[ps_e,phase_e,rate_e]=sample_estimate(data,Fs,fc_e,B_e);%Joint estimation method
Shape_b=rcosdesign(0.6,100,round(Fs/ps_e),'sqrt');
[demod]=QAM_demode(data,ps_e,Fs,fc_e,Shape_b,1,0);%Baseband demodulation
dd_new=scplot(demod,phase_e,Fs,ps_e);%Constellation diagram
R_max=r_max_get(demod);


subplot(2,3,3);
demode_plot(dd_new,R_max,1.2);
title("The constellation diagram with frequency offset")

[F_zz,F_z,m_z]=fre_correct(dd_new,fc_e,ps_e);%Frequency offset correction


dd_final=fc_cor_demode_plot(dd_new,ps_e,m_z);%Constellation diagram
subplot(2,3,6);
demode_plot(dd_final,R_max,1.2);
title("Final constellation diagram") 











function demode_plot(data_in,R_max,bei)
    hold on;
    for i=1:length(data_in)
        plot(real(data_in(i)),imag(data_in(i)),'r.');
    end
    axis([-bei*R_max bei*R_max -bei*R_max bei*R_max]);
    hold off;
end



function R_max=r_max_get(data_in)
    R_max=0;
    for i=1:length(data_in)
        R=sqrt(real(data_in(i))^2+imag(data_in(i))^2);
        if(R>R_max)
            R_max=R;
        end
    end
end


function dd_final=fc_cor_demode_plot(dd_new,ps_e,m_z)
     dt=1/(ps_e);
    for i=1:length(dd_new)
        R=sqrt(real(dd_new(i))^2+imag(dd_new(i))^2);
        theta=angle(dd_new(i));
        I=R*cos(theta+dt*(i-1)*m_z*2*pi);
        Q=R*sin(theta+dt*(i-1)*m_z*2*pi);
        dd_final(i)=I+sqrt(-1)*Q;
    end
end
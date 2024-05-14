function  [ps_m,phase_m,rate_e,X_ff]=sample_estimate(x,Fs,fc,B)  

    x=x(1:floor(length(x)));
    t=0:length(x)-1;%
    t=t/Fs;
    
    f0_di=cos(2*pi*(fc)*t); 
    f0_dq=sin(2*pi*(fc)*t);  
    %
    demod_mult_i=x.*f0_di;
    demod_mult_q=x.*f0_dq;
    %
    demod=demod_mult_i+sqrt(-1)*demod_mult_q;


    xfft=fft(x.^2+hilbert(x).^2);
    x_begin=B/2;
    x_end=B/1;
    x_begin=fix(x_begin*length(x)/Fs);
    x_end=fix(x_end*length(x)/Fs);
%     plot(Fs*(x_begin:x_end)/length(x),abs(xfft(x_begin:x_end)));
   
    xfft_max=0;
    for i=x_begin:x_end
        if(abs(xfft(i))>xfft_max)
            xfft_max=abs(xfft(i));
            ps_m=i*Fs/length(x);
        end
    end
    
    num3=fix(Fs/ps_m);
    ps_j1=100;
    theta_j1=40;
    Shape_b2=rcosdesign(0.6,100,num3,'sqrt');
    demodp=filter(Shape_b2,1,demod);
    ps_new1=ps_m*0.95;
    ps_new2=ps_m*1.05;

    cnt=floor((ps_new2-ps_new1)/ps_j1);
    cnt_phase=fix((num3)/theta_j1);
    sum123=[];
    summ=[];
    ttt=0;
    mmax=0;
    phase_m=0;
    ps_m=0;
    for phase_new=0:cnt_phase:theta_j1*cnt_phase
      ttt=ttt+1;
      tt=0;
        demod_t=demodp(1+phase_new:length(demodp));
    for ps_new=ps_new1:cnt:cnt*ps_j1+ps_new1
        
        rate_c=Fs/ps_new;
        tt=tt+1;
        changdu_new=floor(length(demod_t)/rate_c);
        demod_new=[];
        for i=1:changdu_new
            demod_new=[demod_new demod_t(fix(i*rate_c))];
        end
        sum123(tt)=(std(real(demod_new),1)^2+std(imag(demod_new),1)^2);
        if sum123(tt)>=mmax
            ps_m=ps_new;
            phase_m=phase_new;
            mmax=sum123(tt);
        end
    end
    summ(ttt,:)=sum123;
    end

    x_la=ps_new1:cnt:cnt*ps_j1+ps_new1;
    y_la=0:cnt_phase:theta_j1*cnt_phase;
    subplot(2,3,1);
    mesh(x_la,y_la,summ)
    title("Joint estimation - Large scale")



% ps_m=10^6;


    num3=fix(Fs/ps_m);
    ps_j1=200;
    theta_j1=num3;
    Shape_b2=rcosdesign(0.6,100,num3,'sqrt');
    demodp=filter(Shape_b2,1,demod);
%     ps_new1=floor(Fs/(num3+1));
%     ps_new2=floor(Fs/(num3-1));

    ps_new1=ps_m*0.995;
    ps_new2=ps_m*1.005;

    cnt=floor((ps_new2-ps_new1)/ps_j1);
    cnt_phase=1;
    sum123=[];
    summ=[];
    ttt=0;
    mmax=0;
    phase_m=0;
    ps_m=0;
    for phase_new=0:cnt_phase:fix(theta_j1/cnt_phase)*cnt_phase
      demod_t=demodp(1+phase_new:length(demodp));
      ttt=ttt+1;
      tt=0;
    for ps_new=ps_new1:cnt:cnt*ps_j1+ps_new1
        rate_c=Fs/ps_new;
        tt=tt+1;
        changdu_new=floor(length(demod_t)/rate_c);
        demod_new=[];
        for i=1:changdu_new
            demod_new=[demod_new demod_t(fix(i*rate_c))];
        end
%         sum123(tt)=(std(real(demod_new),1)^2+std(imag(demod_new),1)^2);
%         sum123(tt)=(std(real(demod_new),1));
        sum123(tt)=(mean(real(demod_new).^2+imag(demod_new).^2));
        if sum123(tt)>=mmax
            ps_m=ps_new;
            phase_m=phase_new;
            mmax=sum123(tt);
        end
    end
    summ(ttt,:)=sum123;
    end
%     
%     figure;

    plot_x=1./(ps_new1:cnt:cnt*ps_j1+ps_new1);
%     plot_y=0:cnt_phase:fix(theta_j1/cnt_phase)*cnt_phase;
    plot_y=0:fix(theta_j1/cnt_phase);
    plot_z=summ;
    subplot(2,3,2);
    mesh(plot_x,plot_y,plot_z)
    title("Joint estimation - Small scale")
%     axis([1/ps_new1 1/ps_new2 -inf inf])
%     surf


    rate_c=Fs/ps_m;
    Shape_b2=rcosdesign(0.6,100,fix(rate_c),'sqrt');
    demod=filter(Shape_b2,1,demod);
     demodpt=demod(1+phase_m:length(demod));
      changdu_new=floor(length(demodpt)/rate_c);
        demod_new=[];
        for i=20:changdu_new
            demod_new=[demod_new demodpt(fix(i*rate_c))];
        end

     rate_e=fix(rate_c);
%      if phase_m>=rate_e*0.95
%          phase_m=0;
%      end
%      ps_m
%     phase_m
end
function  [ttt,demodpp]=QAM_demode(qam,ps,Fs,fc,Shape_b,bei,phase_e)  
    rate=fix(Fs/ps);

    qam=qam(1:length(qam));

     t=phase_e:length(qam)-1+phase_e;
    t=t/Fs;
    
    
    f0_di=cos(2*pi*(fc)*t); 
    f0_dq=sin(2*pi*(fc)*t);  
    %
    demod_mult_i=qam.*f0_di;
    demod_mult_q=qam.*f0_dq;
    %
    demodpp=demod_mult_i+sqrt(-1)*demod_mult_q;
      demod=filter(Shape_b,1,demodpp);
     ttt=demod*bei;

end
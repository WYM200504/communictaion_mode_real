function  [demod_new]=scplot(demod,phase,Fs,ps)  

    rate_c=Fs/ps;
%      demodpt=demod(1+phase+fix(20*rate_c):length(demod));
    demodpt=demod(1+phase:length(demod));
      changdu_new=floor(length(demodpt)/rate_c);
        demod_new=[];
        for i=1:changdu_new
            demod_new=[demod_new demodpt(fix(1+(i-1)*rate_c))];
        end
  
end
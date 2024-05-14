function  [F_zz,F_z,m_z]=fre_correct(dd_new,fc_e,ps_e)  
    
    len=10000;
    err=0.0000001;
    for i=1:length(dd_new)
        R(i)=sqrt(real(dd_new(i))^2+imag(dd_new(i))^2);
        theta(i)=angle(dd_new(i));
    end
    max_R=max(abs(R))*1.5;
    dt=1/(ps_e);
    chu=64;
    min_z=1000000000000000;
    F_e=fc_e;
    m_z=F_e;
    for k=-len:len
        m=F_e*err*k;
        z_z=zeros((chu*2+1)*(chu*2+1),1);
        for i=1:length(dd_new)
              I=R(i)*cos(theta(i)+dt*(i-1)*m*2*pi);
              Q=R(i)*sin(theta(i)+dt*(i-1)*m*2*pi);
              I_c=fix(chu*(I+max_R)/max_R);
              Q_c=fix(chu*(Q+max_R)/max_R);
              z_z(I_c+Q_c*(chu*2+1))=1;
        end
        sum_z(k+len+1)=sum(z_z);
        if sum_z(k+len+1)<=min_z
            min_z=sum_z(k+len+1);
            m_z=m;
        end
    end
    F_zz=[];
    count_n=0;
    avg_m=mean(sum_z);
    min_m=avg_m-min_z;
    for i=1:length(sum_z)
        count_n=count_n+1;
        if count_n>=2*len
                break;
        end
        if sum_z(count_n)<=avg_m-min_m*0.8
            m_kz=count_n;
            m_ky=count_n+len/5;
            if m_ky>=2*len
                m_ky=2*len;
            end
       
            m_zzn=find(sum_z(m_kz:m_ky)==min(sum_z(m_kz:m_ky)));
            
        
            m_zz=F_e+F_e*err*(m_zzn+m_kz-len-1);
            F_zz=[F_zz m_zz];
            count_n=count_n+len/5;
           
        end
    end
    if length(F_zz)>=100
        F_zz=F_zz(1:100);
    end
    sum_z=sum_z/chu/chu;
%     figure

subplot(2,3,4);
    plot(-len*F_e*err+F_e:F_e*err:len*F_e*err+F_e,sum_z,LineWidth=1)
    title("Constellation dispersion - Large scale")
%     axis([-len*F_e*err+F_e len*F_e*err+F_e 0.07 0.185])
    F_e=F_e+m_z;

m_z1=m_z;






    len=10000;
    err=0.000000001;
    for i=1:length(dd_new)
        R(i)=sqrt(real(dd_new(i))^2+imag(dd_new(i))^2);
        theta(i)=angle(dd_new(i))+dt*(i-1)*m_z*2*pi;
    end
    min_z=1000000000000000;
    m_z=F_e;
    for k=-len:len
        m=F_e*err*k;
        z_z=zeros((chu*2+1)*(chu*2+1),1);
        for i=1:length(dd_new)
              I=R(i)*cos(theta(i)+dt*(i-1)*m*2*pi);
              Q=R(i)*sin(theta(i)+dt*(i-1)*m*2*pi);
              I_c=fix(chu*(I+max_R)/max_R);
              Q_c=fix(chu*(Q+max_R)/max_R);
              z_z(I_c+Q_c*(chu*2+1))=1;
        end
        sum_z2(k+len+1)=sum(z_z);
        if sum_z2(k+len+1)<=min_z
            min_z=sum_z2(k+len+1);
            m_z=m;
        end
    end
    F_zz=[];
    count_n=0;
    avg_m=mean(sum_z2);
    min_m=avg_m-min_z;
    for i=1:length(sum_z2)
        count_n=count_n+1;
        if count_n>=2*len
                break;
        end
        if sum_z2(count_n)<=avg_m-min_m*0.8
            m_kz=count_n;
            m_ky=count_n+len/5;
            if m_ky>=2*len
                m_ky=2*len;
            end
       
            m_zzn=find(sum_z2(m_kz:m_ky)==min(sum_z2(m_kz:m_ky)));
            
        
            m_zz=F_e+F_e*err*(m_zzn+m_kz-len-1);
            F_zz=[F_zz m_zz];
            count_n=count_n+len/5;
           
        end
    end
    if length(F_zz)>=100
        F_zz=F_zz(1:100);
    end
    sum_z2=sum_z2/chu/chu;
%     figure
subplot(2,3,5);
    plot(-len*F_e*err+F_e:F_e*err:len*F_e*err+F_e,sum_z2,LineWidth=1)
    title("Constellation dispersion - Small scale")
%     axis([-len*F_e*err+F_e len*F_e*err+F_e 0.07 0.185])
    F_z=F_e+m_z;
    m_z=m_z+m_z1;
end
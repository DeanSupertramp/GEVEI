for i = 2:(length(t))

        x(i) = x(i-1) + dt*(v(i)+v(i-1))/2;
        mu_r=mu_r0*(1+v(i)/100);
        Rr0=mu_r*m*g;
        Rr=Rr0*cosd(alpha); %resistenza al rotolamento
        Ftr(i)=Rr+Rg+(0.5*A*Cd*ro*(v(i)^2)+(m_eff*a(i)));


        Ftr_p(i)=max(0,Ftr(i));
        Ftr_n(i)=abs(min(0,Ftr(i)));

        P_p(i)=Ftr_p(i)*v(i);
        P_n(i)=Ftr_n(i)*v(i);

        P_eff_p(i)=P_p(i)/eff_mech;%P_eff=potenza effettiva da dividere tra IC e EL
        P_eff_n(i)=P_n(i)*eff_mech;

        w(i)=v(i)/r*G0; %w motore
        v_rpm=(60/(2*pi))*w(i);



        if w(i)==0
            T_p(i)=0;
            T_n(i)=0;
            Tmax_el_w=0;
        else
            T_p(i)=P_eff_p(i)/w(i);
            T_n(i)=P_eff_n(i)/w(i);
            if(w(i)>wc_el)
                Tmax_el_w=Pmax_el/w(i);
            else
                Tmax_el_w=Tmax_el;
            end
        end    



        if(Ftr(i)>0)
           if((T_p(i)/(G*G0))>Tmax_el_w)
               Tstop=Tmax_el_w;
           else
               Tstop=T_p(i)/(G*G0);
           end
           T_opt=linspace(0,Tstop,1000);
           T_ice=T_p(i)/(G*G0)-T_opt;     
           eff_motore_ic_vec=mappa_eta_ICE(w(i),T_ice,E0);
           Pfuel_vec=(T_ice*w(i))./eff_motore_ic_vec;
           eff_motore_el_vec=mappa_eta_EL_mot(w(i),T_opt,E0);
           P_batt_p_vec=T_opt*w(i)./eff_motore_el_vec;
           I_p_vec=P_batt_p_vec./V(i-1);
           D_SOC_p_vec=((I_p_vec).^k)*dt/Q0;
           SOC_vec=SOC(i-1)-D_SOC_p_vec;

           V_vec=215+m1.*(SOC_vec);
           P_pk_p_vec=I_p_vec.^k .* V_vec;
           SOCp_vec=-(1/Q0)*(V_vec/(2*R0)-sqrt((V_vec/(2*R0)).^2 - ((P_pk_p_vec)/R0)));
           Pech_vec=-V_vec*Q0.*SOCp_vec;
           w_SOC=5*(0.9-SOC(i-1));
           H_vec=Pfuel_vec+(L1(i)+w_SOC).*Pech_vec;
           

           [Hmin,min_index]=min(H_vec);
           Tel_p(i)=T_opt(min_index);
           Tel_n(i)=0;

        else
           if(T_n(i)>Tmax_el_w)
               Tstop=-Tmax_el_w;
           else
               Tstop=-T_n(i);
           end
           T_opt=linspace(Tstop,0,1000);
           eff_motore_el_vec=mappa_eta_EL_gen(w(i),T_opt,E0);
           P_batt_n_vec=T_opt*w(i)./eff_motore_el_vec;
           I_n_vec=-P_batt_n_vec./V(i-1);
           D_SOC_n_vec=((I_n_vec).^(1/k))*dt/Q0;
           SOC_vec=SOC(i-1)+D_SOC_n_vec; 
           SOC_vec(find(SOC_vec>1))=1;
           V_vec=215+m1.*(SOC_vec);
           P_pk_n_vec=I_n_vec.^(1/k) .* V_vec;
           SOCp_vec=-(1/Q0)*(V_vec/(2*R0)-sqrt((V_vec/(2*R0)).^2 - ((-P_pk_n_vec)/R0)));
           Pech_vec=-V_vec*Q0.*SOCp_vec;
           w_SOC=5*(0.9-SOC(i-1));
           H_vec=(L1(i)+w_SOC).*Pech_vec;

           [Hmin,min_index]=min(H_vec);
           Tel_p(i)=0;
           Tel_n(i)=-T_opt(min_index);

        end

        
        if w(i)==0
            T_p(i)=0;
            T_n(i)=0;
            Tic(i)=0;
            Tel_p(i)=0;
            Tel_n(i)=0;
            eff_motore_el(i)=0;
            eff_gen(i)=0;
        else
            Tic(i)=(T_p(i)/(G*G0))-Tel_p(i);
            eff_motore_el(i)=mappa_eta_EL_mot(w(i),Tel_p(i),E0);
            eff_gen(i)=mappa_eta_EL_gen(w(i),Tel_n(i),E0);
            eff_motore_ic(i)=mappa_eta_ICE(w(i),Tic(i),E0);            
            P_ic(i)=Tic(i)*w(i)/eff_motore_ic(i);          
        end

        
        
        P_batt_p(i)=Tel_p(i)*w(i)/eff_motore_el(i);
        if(P_batt_p(i)>PbattMax)
            P_batt_p(i)=PbattMax;
        end
        if(isnan(P_batt_p(i)))
            P_batt_p(i)=0;
        end
        P_batt_n(i)=Tel_n(i)*w(i)*eff_gen(i); %calcolo la potenza fornita alla batteria
        if(isnan(P_batt_n(i)))
            P_batt_n(i)=0;
        end

        I_p(i)=P_batt_p(i)/V(i-1); % corrente erogata dalla batteria
        if(I_p(i)>IbattMax)
            I_p(i)=IbattMax;
        end

        I_n(i)=P_batt_n(i)/V(i-1); % corrente fornita alla batteria

        D_SOC_p(i)=((I_p(i))^k)*dt/Q0;
        D_SOC_n(i)=(I_n(i)^(1/k))*dt/Q0;
        %D_SOC_n(i)=0; %senza frenata rigenerativa

        SOC(i)=SOC(i-1)-D_SOC_p(i)+D_SOC_n(i);     % Con frenata rigenerativa
        if(SOC(i)>1)
            SOC(i)=1;
        end
    
        
        
        Q(i)=Q0*SOC(i);                 %carica all'istante i
        V(i)=215+m1*(SOC(i));  %tensione nominale all'istante i

        P_pk_n(i)=I_n(i).^(1/k) .* V(i);
        P_pk_p(i)=I_p(i).^(k) .* V(i);

        dSOCp(i)=((5e-4)/(35e3))*((P_pk_p(i)-P_pk_n(i))/1e3);
        L1(i+1)=L1(i)-dt*L1(i)*dSOCp(i);

    end
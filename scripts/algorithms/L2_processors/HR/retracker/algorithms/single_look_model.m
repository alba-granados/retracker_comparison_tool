function stack=single_look_model(x,l, fit_p, nf_p, cnf_p,func_f0,func_f1)
        if nargin < 7
            func_f1=0; % initialize altitude to 0;
        end
        if nargin <6
            func_f0=0; % initialize altitude to 0;
        end
        
        epoch=fit_p(1);
        
        if cnf_p.rou_flag
            sigmaz=1e-2/4; rou=fit_p(2);
        else
            rou=nf_p.rou; sigmaz=fit_p(2);
        end
        %------------dilation term---------------------------------------------
        g       =   sqrt(2*nf_p.alphag_a*nf_p.alphag_r./(nf_p.alphag_a + nf_p.alphag_r * 4 * (nf_p.Lx/nf_p.Ly)^4 * l.^2 + 2 * nf_p.alphag_a * nf_p.alphag_r * (sigmaz/nf_p.Lz)^2 ));
        
        %--------------Varibales  ---------------------------------------------
        k=(x - epoch); xl=nf_p.Lx.*l; alpha_sigma=1/(nf_p.h^2*rou); yk=nf_p.Ly.*abs(sqrt(k)); yk(k<0)=0.0; gk=g.*k;
        %
        %**********************************************************************
        %************************ Antenna & Surface ***************************
        %**********************************************************************
        %Constant Term        
        Bkl=2.0*exp(-nf_p.alphax *(xl-nf_p.xp).^2-alpha_sigma * xl.^2-nf_p.alphay * nf_p.yp^2-(nf_p.alphay + alpha_sigma).*(yk).^2).*cosh(2*nf_p.alphay*nf_p.yp*yk);        
        
        % Linear Term
        %     Tkl(k~=0)=(nf_p.Ly./abs(sqrt(k(k~=0)))).*(nf_p.alphay*nf_p.yp).*tanh(2*nf_p.alphay*nf_p.yp*yk(k~=0))-...
        %                 (nf_p.alphay + alpha_sigma)*nf_p.Ly^2;
        %small angle approximation for 0
        Tkl=nf_p.Ly^2.*(nf_p.alphay*nf_p.yp)^2*2-(nf_p.alphay + alpha_sigma)*nf_p.Ly^2;
        
        %**********************************************************************
        %************************ Dilation functions **************************
        %**********************************************************************
        funcf0=0.*gk;
        if cnf_p.lut_flag
            funcf0(gk==0)=1.077900274770464;
            indexes_1=gk>=cnf_p.LUT_ximin & gk<=cnf_p.LUT_ximax;
            indexes_2=floor((gk(indexes_1)-cnf_p.LUT_ximin)./cnf_p.LUT_step)+1;
            indexes_3=gk>cnf_p.LUT_ximax;
            funcf0(indexes_1)=func_f0(indexes_2);
            if any(indexes_3)
                funcf0(indexes_3)=sqrt(pi./(2.0*gk(indexes_3))).*(1+3./(8*(gk(indexes_3)).^2)+105./(16*8*(gk(indexes_3)).^4));
            end            
            switch cnf_p.power_wfm_model
                case 'complete'
                    funcf1=0.*gk;
                    funcf1(gk==0)=0.515224256147498;
                    funcf1(indexes_1)=func_f1(indexes_2);
                    if any(indexes_3)
                        funcf1(indexes_3)=-1.0*sqrt(pi.*gk(indexes_3)/8).*(1./((gk(indexes_3)).^2)+15./(8.*(gk(indexes_3)).^4));
                    end                    
            end
        else
            funcf0(gk~=0)=pi/4.0*sqrt(abs(gk(gk~=0))).*(besseli(-1/4,1/4*(gk(gk~=0)).^2,1)+sign(gk(gk~=0)).*besseli(1/4,1/4*(gk(gk~=0)).^2,1)); funcf0(gk==0)=1.077900274770464;%2^(1/4)*gamma(5/4);
            switch cnf_p.power_wfm_model
                case 'complete'
                    funcf1=zeros(1,length(x));
                    funcf1(gk~=0)=-1.0*pi/8.0*(abs(gk(gk~=0)).^(3/2)).*(besseli(1/4,1/4*(gk(gk~=0)).^2,1)-besseli(-3/4,1/4*(gk(gk~=0)).^2,1)+sign(gk(gk~=0)).*(besseli(-1/4,1/4*(gk(gk~=0)).^2,1)-besseli(3/4,1/4*(gk(gk~=0)).^2,1)));
                    funcf1(gk==0)=0.515224256147498;%gamma(3/4)/(2.0*(2)^(1/4));
            end
            
        end
        
        %**********************************************************************
        %************************ Common Cte param ****************************
        %**********************************************************************
        %K_cte=sqrt(2*pi)*sqrt(1/(2*alphag_r))*sqrt(1/(2*alphag_a))*;
        
        
        
        %**********************************************************************
        %************************ Power Waveform ******************************
        %**********************************************************************
        switch cnf_p.power_wfm_model
            case 'simple'
                stack         =   sqrt(g).*Bkl.*funcf0;
            case 'complete'
                stack         =   sqrt(g).*Bkl.*(funcf0+Tkl.*g.*((sigmaz/nf_p.Lz)^2).*funcf1);
        end
    end
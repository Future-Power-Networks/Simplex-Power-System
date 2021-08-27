% This class defines the model of a grid-forming voltage source inverter.

% Author(s): Yitong Li

%% Notes
%
% The model is in load convention.
% 
% The model is in admittance form.
%
% dw means the derivative of w

%% Class

classdef GridFormingVSI < SimplusGT.Class.ModelAdvance
    
    % For temporary use
    properties(Access = protected)
        v_od_r;
        v_oq_r;
        P0;
        Q0;
    end
    
    methods
        % constructor
        function obj = GridFormingVSI(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
         	State  = {'i_ld','i_lq','i_ld_i','i_lq_i','v_od','v_oq','v_od_i','v_oq_i','i_od','i_oq','v_d_ref','w','theta'};
            Input  = {'v_d','v_q','P0'};
            Output = {'i_d','i_q','w','theta'};
        end
        
        % Calculate the equilibrium
        function [x_e,u_e,xi] = Equilibrium(obj)
         	% Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
            
            % Get parameter
            xwLf    = obj.Para(1);
            Rf      = obj.Para(2);
            xwCf    = obj.Para(3);
            xwLc    = obj.Para(4);
            Rc      = obj.Para(5);
            Xov     = obj.Para(6);
            Rov     = 0;
            xDw     = obj.Para(7);
            Dv      = 0.05;
            xfdroop = obj.Para(8);
            xfvdc   = obj.Para(9);
            xfidq   = obj.Para(10);
            W0      = obj.Para(11);
           
            % Calculate parameters
            Lf = xwLf/W0;
            Cf = xwCf/W0;
            Lc = xwLc/W0;
            Dw = xDw*W0;
            wf = xfdroop*2*pi;
            w_v_odq = xfvdc*2*pi;
            w_i_ldq = xfidq*2*pi;
            
            v_gd = V;
            v_gq = 0;
            i_od = P/V;
            i_oq = -Q/V;
            
            v_gdq = v_gd + 1i*v_gq;
            i_odq = i_od + 1i*i_oq;
            v_odq = v_gdq - i_odq*(Rc + 1i*w*Lc);
            i_cdq = v_odq*(1i*w*Cf);
            i_ldq = i_odq - i_cdq;
            e_dq  = v_odq - i_ldq*(Rf + 1i*w*Lf);
            
            i_ld = real(i_ldq);
            i_lq = imag(i_ldq);
            i_ld_i = real(e_dq);
            i_lq_i = imag(e_dq);
            v_od = real(v_odq);
            v_oq = imag(v_odq);
            v_od_i = i_ld;
            v_oq_i = i_lq;
            i_od = real(i_od);
            i_oq = imag(i_oq);
            theta = xi;
            
            obj.P0 = P*(-1);
            obj.Q0 = Q*(-1);
            
            % ??? Temp
            v_odq_r = v_odq + (Rov + 1i*Xov)*i_odq*(-1);
            v_od_r = real(v_odq_r);
            v_oq_r = imag(v_odq_r);
            obj.v_od_r = v_od_r;
            obj.v_oq_r = v_oq_r;
            
            % Notes:
            % Ideally, v_oq_r = 0 should be valid. So, this equilibrium
            % calculation has to be corrected.
            % The P-F and Q-V droop has not been considerred, which should
            % also be corrected.
            
            % Get equilibrium
            x_e = [i_ld; i_lq; i_ld_i; i_lq_i; v_od; v_oq; v_od_i; v_oq_i; i_od; i_oq; v_od_r; w; theta];
            u_e = [v_gd; v_gq; 0];
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
        	% Get the power PowerFlow values
            V	= obj.PowerFlow(3);
            
            % Get input
            v_gd   = u(1);
            v_gq   = u(2);

            % Get state
            i_ld   = x(1);
            i_lq   = x(2);
            i_ld_i = x(3);
            i_lq_i = x(4);
            v_od   = x(5);
            v_oq   = x(6); 
            v_od_i = x(7);
            v_oq_i = x(8);
            i_od   = x(9);
            i_oq   = x(10);
            v_od_r = x(11);
            w      = x(12);
            theta  = x(13);
            
            % Get parameters
            xwLf    = obj.Para(1);
            Rf      = obj.Para(2);
            xwCf    = obj.Para(3);
            xwLc    = obj.Para(4);
            Rc      = obj.Para(5);
            Xov     = obj.Para(6);
            Rov     = 0;
            xDw     = obj.Para(7);
            Dv      = 0.05;
            xfdroop = obj.Para(8);
            xfvdq   = obj.Para(9);
            xfidq   = obj.Para(10);
            W0      = obj.Para(11);
            
            % Update paramters
            Lf = xwLf/W0;
            Cf = xwCf/W0;
            Lc = xwLc/W0;
            Dw = xDw*W0;
            wf = xfdroop*2*pi;
            w_v_odq = xfvdq*2*pi;
            w_i_ldq = xfidq*2*pi;
            
            kp_i_ldq = w_i_ldq*Lf;
            ki_i_ldq = w_i_ldq^2*Lf/4;
            kp_v_odq = w_v_odq*Cf;
            ki_v_odq = w_v_odq^2*Cf/4*50;
                % This is a different way of setting voltage PI
                % kp_v_odq = 1/(16*w_i_ldq*Lf);
                % ki_v_odq = 1/(4*Lf);
            
            % Saturation setting
            EnableSaturation = 1;
            
            % Frequency limit and saturation
            w_limit_H = W0*1.1;
            w_limit_L = W0*0.9;
            % Capacitor voltage limit
            v_od_limit_H = 1.5;
            v_od_limit_L = -1.5;
            v_oq_limit_H = 1.5;
            v_oq_limit_L = -1.5;    
            % Current reference limit
            i_ld_limit = 1.5;
            i_lq_limit = 1.5;
            % Ac voltage limit
            e_d_limit_H = 1.5;
            e_d_limit_L = -1.5;
            e_q_limit_H = 1.5;
            e_q_limit_L = -1.5;
       
  
            % State space equations
         	% dx/dt = f(x,u)
            % y     = g(x,u)
            %if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)
            
                % Power measurement
                p = (v_od*i_od + v_oq*i_oq)*(-1);   % (-1) appears because the model is in load convention
                q = (-v_od*i_oq + v_oq*i_od)*(-1);

                % Droop controller
                % Standard droop equations:
                % PF droop: w = w0 + Dw*(P0 - P*LPF)
                % QV droop: v_od_r = v_od_0 + Dv*(Q0 - Q*LPF)
                %           v_oq_r = v_oq_0
                P0     = obj.P0;
                Q0     = obj.Q0;
                if 1
                    dw = (W0 + Dw*(P0 - p) - w)*wf;         % P-w droop
                else
                    dw = (W0 - Dw*(P0/V - i_od) - w)*wf; 	% id-w droop
                end
                
                % Limitation for w
                if EnableSaturation
                     if (w >= w_limit_H && dw >=0 ) || (w <= w_limit_L && dw <= 0 )
                        dw = 0;                  
                     end
                end
                
                switch 3
                    case 1
                        dv_od_r = (obj.v_od_r + Dv*(Q0 - q) - v_od_r)*wf;   % Q-V droop
                        % Limitation for v_od_r
                        if EnableSaturation
                             if (v_od_r >= v_od_limit_H && dv_od_r >=0 ) || (v_od_r <= v_od_limit_L && dv_od_r <= 0 )
                                dv_od_r = 0;                  
                             end
                        end                        
                    case 2
                        v_od_r = obj.v_od_r + Dv*(Q0 - q);                  % Q-V droop without LPF
                        dv_od_r = 0;
                        % Limitation for v_od_r
                        if EnableSaturation
                            v_od_r = min(v_od_r,v_od_limit_H);
                            v_od_r = max(v_od_r,v_od_limit_L);
                        end  
                    case 3
                        v_od_r = obj.v_od_r;                                % No Q-V droop
                        dv_od_r = 0;
                        % Limitation for v_od_r
                        if EnableSaturation
                            v_od_r = min(v_od_r,v_od_limit_H);
                            v_od_r = max(v_od_r,v_od_limit_L);
                        end                        
                    otherwise
                        error(['Error'])
                end
                v_oq_r = obj.v_oq_r;
                % Limitation for v_oq_r
                if EnableSaturation
                    v_oq_r = min(v_oq_r,v_oq_limit_H);
                    v_oq_r = max(v_oq_r,v_oq_limit_L);
                end  
                
                % AC voltage control
                error_v_od = v_od_r - v_od - (i_od*Rov-i_oq*Xov)*(-1);
              	error_v_oq = v_oq_r - v_oq - (i_oq*Rov+i_od*Xov)*(-1);
                i_ld_r = -(error_v_od*kp_v_odq + v_od_i);
                i_lq_r = -(error_v_oq*kp_v_odq + v_oq_i);
                dv_od_i = error_v_od*ki_v_odq;
                dv_oq_i = error_v_oq*ki_v_odq;
                % AC voltage controller anti-windup
                if EnableSaturation
                     if (v_od_i >= i_ld_limit && dv_od_i >=0 ) || (v_od_i <= -i_ld_limit && dv_od_i <= 0)
                        dv_od_i = 0;                  
                     end
                     if (v_oq_i >= i_lq_limit && dv_oq_i >= 0 ) || (v_oq_i <= -i_lq_limit && dv_oq_i <= 0)
                        dv_oq_i = 0;                  
                     end                    
                end
                % Current saturation
                if EnableSaturation
                    i_ld_r = min(i_ld_r,i_ld_limit);
                    i_ld_r = max(i_ld_r,-i_ld_limit);
                    i_lq_r = min(i_lq_r,i_lq_limit);
                    i_lq_r = max(i_lq_r,-i_lq_limit);
                end

                % AC current control
                error_i_ld = i_ld_r-i_ld;
                error_i_lq = i_lq_r-i_lq;
                e_d = -error_i_ld*kp_i_ldq + i_ld_i;
                e_q = -error_i_lq*kp_i_ldq + i_lq_i;
                di_ld_i = -error_i_ld*ki_i_ldq;
                di_lq_i = -error_i_lq*ki_i_ldq;
                % Current controller anti-windup
                if EnableSaturation
                     if (i_ld_i >= e_d_limit_H && di_ld_i >=0 ) || (i_ld_i <= e_d_limit_L && di_ld_i <= 0)
                        di_ld_i = 0;
                     end
                     if (i_lq_i >= e_q_limit_H && di_lq_i >= 0 ) || (i_lq_i <= e_q_limit_L && di_lq_i <= 0)
                        di_lq_i = 0;                  
                     end                    
                end
                % Ac voltage (duty cycle) saturation
                if EnableSaturation
                    e_d = min(e_d,e_d_limit_H);
                    e_d = max(e_d,e_d_limit_L);
                    e_q = min(e_q,e_q_limit_H);
                    e_q = max(e_q,e_q_limit_L);
                end

                % Lf equation
                % e_d - v_od = -(di_ld/dt*Lf + Rf*i_ld - w*Lf*i_lq)
                % e_q - v_oq = -(di_lq/dt*Lf + Rf*i_lq + w*Lf*i_ld)
                di_ld = (v_od - e_d - Rf*i_ld + w*Lf*i_lq)/Lf;
                di_lq = (v_oq - e_q - Rf*i_lq - w*Lf*i_ld)/Lf;

                % Cf equation
                % -(i_ld - i_od) = Cf*dv_cd/dt - w*Cf*v_cq
                % -(i_lq - i_oq) = Cf*dv_cq/dt + w*Cf*v_cd
                dv_od = (-(i_ld - i_od) + w*Cf*v_oq)/Cf;
                dv_oq = (-(i_lq - i_oq) - w*Cf*v_od)/Cf;

                % Lc equation
                % v_od - v_d = -(Lc*di_od/dt + Rc*i_od - w*Lc*i_oq)
                % v_oq - v_q = -(Lc*di_oq/dt + Rc*i_oq + w*Lc*i_od)
                di_od = (v_gd - v_od - Rc*i_od + w*Lc*i_oq)/Lc;
                di_oq = (v_gq - v_oq - Rc*i_oq - w*Lc*i_od)/Lc;

                % Phase angle
                dtheta = w;
            if CallFlag == 1     
                % dx
                f_xu = [di_ld; di_lq; di_ld_i; di_lq_i; dv_od; dv_oq; dv_od_i; dv_oq_i; di_od; di_oq; dv_od_r; dw; dtheta];
                Output = f_xu;
                
            elseif CallFlag == 2     
            % ### Call output equations: y = g(x,u)
                g_xu = [i_od; i_oq; w; theta];
                Output = g_xu;
            end
            
        end
        
    end
end
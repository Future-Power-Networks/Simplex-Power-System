clear all
close all
clc

%%
CaseVoltageControl = 'DoubleLoop';
% CaseVoltageControl = 'SingleLoop';
% CaseVoltageControl = 'OpenLoop';

%% Base value
BaseValue();

%% Set parameters
% Grid voltage
syms vgD vgQ wg

% Line impedance
syms Lg Rg

% LC Filter
syms Lf Cf Rf

% Droop control
syms vm wf Dw Dv Pr Qr V0 W0

% Voltage controller
syms kpv kiv

% Current controller
syms kpi kii

% Cross-decoupling gain
Fcdv = 0;
Fcdi = 0;

%% System states
% Droop controller
syms w delta

% Voltage controller
syms vdi vqi

% Current controller
syms idi iqi

% Passive component
syms id iq vd vq igd igq

%%
% Inverse frame transformation
% vgdq = vgDQ * e^{-j*delta}
vgd = vgD*cos(delta) + vgQ*sin(delta);
vgq = -vgD*sin(delta) + vgQ*cos(delta);

% Power calculation
q = vq*igd - vd*igq;
p = vd*id + vq*iq;

% Droop control
% Equations:
dw = ((Pr-p)*Dw + W0 - w)*wf;
% dvm = ((Qr-q)*Dv + V0 - vm)*wf;

% Angle difference between inverter and inf bus
% s*delta = w - wg;
ddelta = w - wg;
% dtheta = w;

switch CaseVoltageControl
    case 'DoubleLoop'
        % Voltage controller
        % dvdi = vm - vd;
        dvdi = vm - vd;
        dvqi = 0 - vq;
        idr = kpv*dvdi + kiv*vdi - Fcdv*Cf*Wbase*vq;
        iqr = kpv*dvqi + kiv*vqi + Fcdv*Cf*Wbase*vd;

        % Current controller
        didi = idr - id;
        diqi = iqr - iq;
        ed = kpi*didi + kii*idi - Fcdi*Lf*Wbase*iq;
        eq = kpi*diqi + kii*iqi + Fcdi*Lf*Wbase*iq;
    case 'SingleLoop'
        
    case 'OpenLoop'
        % Open loop voltage control
        ed = vm;
        eq = 0;
end

% Inverter-side inductor
did = (ed - vd + w*Lf*iq - Rf*id)/Lf;
diq = (eq - vq - w*Lf*id - Rf*iq)/Lf;

% Filter capacitor
dvd = (id-igd + w*Cf*vq)/Cf;
dvq = (iq-igq - w*Cf*vd)/Cf;

% Grid-side inductor
digd = (vd - vgd + w*Lg*igq - Rg*igd)/Lg;
digq = (vq - vgq - w*Lg*igd - Rg*igq)/Lg;

% Is this modeling right? Or should I model Lf, Cf also in wg?

%% Calculate the state matrix
switch CaseVoltageControl
    case 'DoubleLoop'
        state = [vdi; vqi; idi; iqi; id; iq; vd; vq; igd; igq; w; delta];
        f_xu = [dvdi; dvqi; didi; diqi; did; diq; dvd; dvq; digd; digq; dw; ddelta];
    case 'SingleLoop'
    case 'OpenLoop'
        state = [id; iq; vd; vq; igd; igq; w; delta];
        f_xu = [did; diq; dvd; dvq; digd; digq; dw; ddelta];
end

Amat = jacobian(f_xu,state);

%% Set numerical number
Cf = 0.02/Wbase;
Lf = 0.05/Wbase;
Rf = 0.01;

wf = 2*pi*10;

wv = 250*2*pi;
kpv = Cf*wv;
kiv = Cf*wv^2/4*20;

wi = 1000*2*pi;
kpi = Lf*wi;
kii = Lf*(wi^2)/4;

Dw = 0.05*Wbase/Sbase;
Dv = 0;

% Dw = Dw*10;

vd = 1.1509;
vq = 0;
P = 0.5;
Q = 0.5;
igd = P/vd;
igq = -Q/vd;
igD = igd;
igQ = igq;
id = igd;
iq = igq;
vm = vd;

delta = 5.9874/180*pi;

Xg = 0.3;
Lg = Xg/Wbase;
Rg = Xg/5;

vgD = 1;
vgQ = 0;

%% Replace symbolic by numerical number

Amat = subs(Amat,'kpi',kpi);
Amat = subs(Amat,'kii',kii);

Amat = subs(Amat,'Dw',Dw);
Amat = subs(Amat,'Dv',Dv);

Amat = subs(Amat,'vd',vd);
Amat = subs(Amat,'vq',vq);
Amat = subs(Amat,'delta',delta);

Amat = subs(Amat,'igd',igd);
Amat = subs(Amat,'igq',igq);

Amat = subs(Amat,'vgD',vgD);
Amat = subs(Amat,'vgQ',vgQ);

Amat = subs(Amat,'id',id);
Amat = subs(Amat,'iq',iq);

Amat = subs(Amat,'wf',wf);

Amat = subs(Amat,'Cf',Cf);
Amat = subs(Amat,'Lf',Lf);
Amat = subs(Amat,'Rf',Rf);

Amat = subs(Amat,'w',Wbase);
Amat = subs(Amat,'wg',Wbase);

%% Sweep parameters
Amat = subs(Amat,'Rg',Rg);
Amat = subs(Amat,'Lg',Lg);

Amat = subs(Amat,'kpv',kpv);
Amat = subs(Amat,'kiv',kiv);

EigVec = eig(Amat);
EigVecHz = EigVec/(2*pi);
ZoomInAxis = [-20,10,-60,60];
PlotPoleMap(EigVecHz,ZoomInAxis,9999);


% ScaleFactor = logspace(-1,2,10);
% for i = 1:length(ScaleFactor)
% Amat_ = subs(Amat,'kiv',kiv*ScaleFactor(i));
% 
% Amat_ = double(Amat_);
% 
% % Calculate poles
% EigVec = eig(Amat_);
% EigVecHz = EigVec/(2*pi);
% 
% % Plot poles
% ZoomInAxis = [-20,10,-60,60];
% PlotPoleMap(EigVecHz,ZoomInAxis,9999);
% end

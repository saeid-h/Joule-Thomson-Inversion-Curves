
%% ------------------------------------------------------------------------
%
%                       Project Specifications
%
% -------------------------------------------------------------------------
% This code is provided to draw Joule-Thomson inversion curve based on  
% Peng-Robinson equation of state for binary mixture of C2 & nC5
%
% Reference:
%       G.W. Dilay, R.A. Heidemann, Calculation of Joule?Thomson 
%       inversion curves from equations of state, 
%       Ind. Eng. Chem. Fund. 25 (1986) 152?158.
%
%
% PR EOS:
%       P = RT / (v-b) - a / (v^2 + 2*b*v-b^2)
%
% JT Eq:
%       T dP/dT + v dP/dv == 0
% 
% Author:   Saeid Hosseinipoor
% email:    saied@ou.edu
% Coauthor: Shadi Salahshoor
% email:    s.salahshoor@ou.edu
%
% Created:  March 9, 2017
% Modified: March 11, 2017
% -------------------------------------------------------------------------


clear all
close all
clc
syms T

%% ------------------------------------------------------------------------
%
%                       Constants definition
%
% -------------------------------------------------------------------------

k  = [0 0.020;
      0.020 0];         % binary interaction coefficient
omega = [0.0995 0.2515];% acentric coefficient
Pc = [48.72 33.70];     % critical pressure (Bar)
Tc = [305.32 469.7];    % critical temperature (K)
vc = [0.1455 0.3130];   % critical volume (m3/kg.mol)
MW = [30.070 72.150];   % Molecular weight (kg /kg.mol)

R  = 8.3144598;         % gas constant (m3.bar/K.mol)


%% ------------------------------------------------------------------------
%
%                       Strating Points
%
% -------------------------------------------------------------------------


% loop for fraction of C2
for xi = 0:.2:1
    
    Temperature = [];
    Pressure = [];
    
    % find fractions for binary mixture
    x(1) = xi;
    x(2) = 1- xi;

% ------------------------------------------------------------------------
%                       Parameters calculation
% -------------------------------------------------------------------------

    b = x * (0.077796 .* R .* Tc ./ Pc)';
    ac = 0.457235 .* R.^2 .* Tc.^2 ./ Pc;
    m = [0.3796 1.485 -0.1644 0.01667] * ...
        [ 1 1; omega; omega.^2; omega.^3];
    alpha = (1 + m .* (1 - sqrt(T ./ Tc))).^2;
    ai = ac .* alpha;

    a = 0; 
    for i = 1:length(x)
        for j = 1:length(x)
            a = a + x(i) * x(j) * (1-k(i,j)) ...
                * sqrt(ac(i)*ac(j)*alpha(i)*alpha(j));
        end
    end

    % find dirivatives from EOS to plug in JT equation
    dai_dT = - m .* ac .* sqrt( alpha ./ (T .* Tc));

    da_dT = x(1).^2 .* dai_dT(1) ...
        + x(1).*x(2).*(1-k(1,2))...
        .* (sqrt(ai(2)./ai(1)).*dai_dT(1)+sqrt(ai(1)./ai(2)).*dai_dT(2))...
        + x(2).^2 .* dai_dT(2);

    % loop for assumed specific volumes to find T and P
    for v = [(3.5-2*x(1)):.1:5 6:1:26 30:5:50 60:10:200]
        % Joule-Thompson Eqution by substitution of EOS
        JTPR = R * T / (v-b) ...
        - da_dT * T / (v^2+2*v*b-b^2) ...
        - R*T*v / (v-b)^2 + 2*a*v*(v+b) / (v^2+2*v*b-b^2)^2;
    
        % Peng-Robinson EOS
        PR = R .* T ./ (v-b) - a ./ (v.^2 + 2.*v.*b-b.^2);

        % Solve JT for temperature
        temp = vpasolve (JTPR == 0, T);
        for j = 1:length(temp)
            if ~isreal(temp(j))
                temp(j) = [];
            end
        end
        temp = double(temp);

        % solve PR for P
        P = double(subs(PR, T, temp));

        % consider valid pressures
        if P > 0
            Temperature = [temp Temperature];
            Pressure = [P Pressure];
        end

    end
    
    % plot the results
    plot(Pressure, Temperature)
    ylabel ('Temperature (K)')
    xlabel ('Pressure (Bar)')
    title ('Joule-Thomson Curves for Binary Mixture of C2 and nC5')
    hold on
end

legend '0-100' '2-80' '40-60' '80-20' '100-0'







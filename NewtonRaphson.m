% Método de Newton-Raphson para fluxo de potência
% Disciplina: Análise de Sistemas Elétricos
% Equipe discente: Davi, Devid e Jonathan

clc; clear;

% Admitâncias série e shunt das linhas
Z12 = 0.01060+1j*0.07830;
Y12  = 1/Z12;
Sh12 = 1j*0.1134/2;
Z13 = 0.00850+1j*0.0335;
Y13  = 1/Z13;
Sh13 = 1j*0.0610/2;
Z24  = 0.00780+1j*0.05700;
Y24  = 1/Z24;
Sh24 = 1j*0.0991/2; 
Z34 = 0.01695+1j*0.06850;
Y34  = 1/Z34;
Sh34 = 1j*0.1213/2;

% Matriz Ybus
Y = [(Y12+Y13+Sh12+Sh13) (-Y12) (-Y13) 0;
     (-Y12) (Y12+Y24+Sh12+Sh24) 0 (-Y24);
     (-Y13) 0 (Y13+Y34+Sh13+Sh34) (-Y34);
      0 (-Y24) (-Y34) (Y24+Y34+Sh24+Sh34)];

G = real(Y);
B = imag(Y);

Barras = 4; % número de barras

V = [1; 1; 1; 1.03];
theta = zeros(Barras,1);

% valores base
s_base = 100; % MVA
v_base = 230; % kV
i_base = (s_base*1e6)/(sqrt(3)*v_base*1e3); % A

% cargas e geradores do sistema
L_1 = 53+1j*27.3;   % carga barra 1
L_2 = 100+1j*30.4;  % carga barra 2
L_3 = 250+1j*60;    % carga barra 3
L_4 = 73+1j*41.55;  % carga barra 4
G_4 = 185.5;        % gerador barra 4 

% potências em cada barra
S = (1/s_base)*[-L_1;-L_2;-L_3;G_4-L_4];

P = real(S); % potência ativa
Q = imag(S); % potência reativa

tol = 0.001;   % tolerância máxima
maxIter = 50;  % nº máximo de iterações

% Considerando:
% Barra 1 = slack
% Barras 2,3 = PQ
% Barra 4 = PV

pq = [2 3 ]; % barras PQ
pv = [4];    % barra PV
pvpq = [pq pv]; % vetor com as barras PV e PQ

npq = length(pq); % tamanho do vetor de barras PQ
npv = length(pv); % tamanho do vetor de barras PV
nvar = npq + npv; % tamanho do vetor das barras combinadas

for iter = 1:maxIter
    
    % Cálculo das potências
    Pcal = zeros(Barras,1);
    Qcal = zeros(Barras,1);
    
    for i = 1:Barras
        for j = 1:Barras
            Pcal(i) = Pcal(i) + V(i)*V(j)*...
                (G(i,j)*cos(theta(i)-theta(j))...
                + B(i,j)*sin(theta(i)-theta(j)));
            Qcal(i) = Qcal(i) + V(i)*V(j)*...
                (G(i,j)*sin(theta(i)-theta(j))...
                - B(i,j)*cos(theta(i)-theta(j)));
        end
    end
    
    % diferente entre os valores de P e Q
    delta_P = P - Pcal;
    delta_Q = Q - Qcal;
    
    % matriz com as diferenças de P e Q
    delta_F = [delta_P(pvpq); delta_Q(pq)];
    
    if max(abs(delta_F)) < tol
        fprintf('Iterações: %d \n', iter);
        break;
    end
    
    % Jacobiana (inicialmente zeradas)
    H = zeros(nvar,nvar);
    Nmat = zeros(nvar,npq);
    M = zeros(npq,nvar);
    L = zeros(npq,npq);
    
% H e N
for a = 1:nvar
    
    i = pvpq(a);
    
    % submatriz H
    for b = 1:nvar
        
        j = pvpq(b);
        
        if i == j            
            H(a,a) = -Qcal(i) - B(i,i)*V(i)^2;
        else            
            H(a,b) = V(i)*V(j)* ...
                ( G(i,j)*sin(theta(i)-theta(j)) ...
                - B(i,j)*cos(theta(i)-theta(j)) );
        end
     end
    
    % submatriz N
    for b = 1:npq
        
        j = pq(b);
        
        if i == j
           Nmat(a,b) = Pcal(i)/V(i) + G(i,i)*V(i);
        else            
            Nmat(a,b) = V(i)* ...
                ( G(i,j)*cos(theta(i)-theta(j)) ...
                + B(i,j)*sin(theta(i)-theta(j)) );
        end
    end
end

% submatriz M
for a = 1:npq
    
    i = pq(a);
    
    for b = 1:nvar
        
        j = pvpq(b);
        
        % M
        if i == j            
            M(a,a) = Pcal(i) - G(i,i)*V(i)^2;
        else            
            M(a,b) = -V(i)*V(j)* ...
                ( G(i,j)*cos(theta(i)-theta(j)) ...
                + B(i,j)*sin(theta(i)-theta(j)) );
        end
    end
    
    % submatriz L
    for b = 1:npq
        
        j = pq(b);
        
        if i == j            
            L(a,a) = Qcal(i)/V(i) - B(i,i)*V(i);
        else            
            L(a,b) = V(i)* ...
                ( G(i,j)*sin(theta(i)-theta(j)) ...
                - B(i,j)*cos(theta(i)-theta(j)) );
        end
    end
end

% matriz jacobiana
jacobiana = [H Nmat;
             M   L];

delta_x = jacobiana \ delta_F;
delta_theta = delta_x(1:nvar);
delta_V = delta_x(nvar+1:end);
    
theta(pvpq) = theta(pvpq) + delta_theta;
V(pq) = V(pq) + delta_V; % PV não atualiza V
    
end

fprintf('\nTensão nas barras:\n');
for i = 1:Barras
    fprintf('Barra %d: V = %.3f pu < %.3f°\n', i, V(i), rad2deg(theta(i)));
end

% Convertendo as tensões nas barras para forma polar
V1 = V(1)*exp(1j*theta(1));
V2 = V(2)*exp(1j*theta(2));
V3 = V(3)*exp(1j*theta(3));
V4 = V(4)*exp(1j*theta(4));

% função para calcular os parâmetros em cada linha
function [fluxo_ik, fluxo_ki, perdas, Iik_real] = sist_ele_pot(Vi,Vk,Zik,Shik,Sbase,Ibase) 
Iik = (Vi - Vk)/Zik + (Vi*Shik); %calcula a corrente da barra i para barra k
Iki = (Vk - Vi)/Zik + (Vk*Shik); %calcula a corrente da barra k para barrra i
fluxo_ik = Sbase*(Vi)*conj(Iik);   %calcula o fluxo da barra i para barra k
fluxo_ki = Sbase*(Vk)*conj(Iki);   %calcula o fluxo da barra k para barrra i
perdas = fluxo_ik + fluxo_ki;      %calcula as perdas na linha
Iik_real = Iik*Ibase;              %calcula a corrente na linha   
end

% linha 1-2
fprintf('\nlinha 1-2');
[fluxo_12, fluxo_21, perdas_12, Iik12] = sist_ele_pot(V1, V2, Z12, Sh12,s_base,i_base);
fprintf('\nCorrente: %.3f A', abs(Iik12));
fprintf('\nFluxo: %.3f MW | %.3f MVAr', real(fluxo_12), imag(fluxo_12));
fprintf('\nPerdas: %.3f MW | %.3f MVAr\n', real(perdas_12), imag(perdas_12));

% linha 1-3
fprintf('\nlinha 1-3');
[fluxo_13,fluxo_31, perdas_13, Iik13] = sist_ele_pot(V1, V3, Z13, Sh13,s_base,i_base);
fprintf('\nCorrente: %.3f A', abs(Iik13));
fprintf('\nFluxo: %.3f MW | %.3f MVAr', real(fluxo_13), imag(fluxo_13));
fprintf('\nPerdas: %.3f MW | %.3f MVAr\n', real(perdas_13), imag(perdas_13));

% linha 2-4
fprintf('\nlinha 2-4');
[fluxo_24,fluxo_42, perdas_24, Iik24] = sist_ele_pot(V2, V4, Z24, Sh24,s_base,i_base);
fprintf('\nCorrente: %.3f A', abs(Iik24));
fprintf('\nFluxo: %.3f MW | %.3f MVAr', real(fluxo_24), imag(fluxo_24));
fprintf('\nPerdas: %.3f MW | %.3f MVAr\n', real(perdas_24), imag(perdas_24));

% linha 3-4
fprintf('\nlinha 3-4');
[fluxo_34,fluxo_43, perdas_34, Iik34] = sist_ele_pot(V3, V4, Z34, Sh34,s_base,i_base);
fprintf('\nCorrente: %.3f A', abs(Iik34));
fprintf('\nFluxo: %.3f MW | %.3f MVAr', real(fluxo_34), imag(fluxo_34));
fprintf('\nPerdas: %.3f MW | %.3f MVAr\n', real(perdas_34), imag(perdas_34));

% potência na barra slack
P_slack = real(fluxo_12) + real(fluxo_13) + real(L_1);

Q_slack = imag(fluxo_12) + imag(fluxo_13) + imag(L_1);
fprintf('\nBarra Slack');
fprintf('\nPotência Ativa: %.3f MVA \nPotência Reativa: %.3f MVar\n ', P_slack, Q_slack); 

% potência reativa no gerador 4
Qg_4 = abs(imag(fluxo_42) + imag(fluxo_43) + imag(L_4));
fprintf('\nPotência reativa gerador 4: %.3f MVar', Qg_4);

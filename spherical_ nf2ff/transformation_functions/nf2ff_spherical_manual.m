function [data_nf2ff] = nf2ff_spherical_manual(data_nf,f,theta_range,phi_range)
% Diese Funktion führt die sphärische Nahfeld-Fernfeld-Transformation 
% durch. Sie berechnet die Fernfeldkomponenten des elektrischen Feldes auf 
% Basis einer zuvor rotierten Nahfelddatentabelle.

% Die Funktion kombiniert die spektrale Zerlegung der Felder mittels 
% sphärischer Vektorfunktionen mit den berechneten Ausstrahlungskoeffizienten 
% Q_{smn}, um das Fernfeld in jeder gewünschten Raumrichtung zu bestimmen.

%
% Input Arguments:
%
%       data_nf           Tabelle der rotierten Nahfelddaten mit 8 Spalten:
%                         X, Y, Z, theta, phi, E (3-spaltig), Eabs
%                         
%                         siehe Funktion rotateSphericalNFData für 
%                         Details zur Struktur der Tabelle
%
%       f                 Frequenz in Hz
%
%       theta_range       Theta-Auswertebereich im Fernfeld (Vektor, Bogenmaß)
%
%       phi_range         Phi-Auswertebereich im Fernfeld (Vektor, Bogenmaß)
%
%
% Output Arguments:
%
%       data_nf2ff        Tabelle mit transformierten Fernfelddaten:
%                         theta, phi, Etheta (komplex), Ephi (komplex), Eabs
%
% Die Funktion erstellt ein Gitter aus allen Theta-/Phi-Kombinationen im 
% Fernfeld, berechnet zu jeder Raumrichtung die Wellenfunktionen und bestimmt 
% die Fernfeldkomponenten durch Projektion und Summation aller Modenbeiträge.
% Das Ergebnis ist eine tabellarisch strukturierte Ausgabe, die als 
% Eingabe für grafische Plots und Validierung verwendet werden kann.

% Wave Number
lambda = physconst('LightSpeed')/f;
k0 = 2*pi/lambda;
% Antenna minimal sphere radius
r0 = 0.25;
% Radius of Measurement Sphere
A = mean(data_nf.r);
theta = data_nf.theta;
phi = data_nf.phi;
Etheta = data_nf.E(:,2);
Ephi = data_nf.E(:,3);

% Step Sizes, Rundung wurde verfeinert früher auf 2 jetzt 4 nachkommastellen
delta_theta = diff(unique(round(data_nf.theta,4)));
delta_theta = delta_theta(1);
delta_phi = diff(unique(round(data_nf.phi,4)));
delta_phi = delta_phi(1);

% Theta/Phi values of Far-field points to calculate
[theta_ff,phi_ff] = meshgrid(theta_range,phi_range);
theta_ff = reshape(theta_ff,numel(theta_ff),1);
phi_ff = reshape(phi_ff,numel(phi_ff),1);

% Far-Field radius
r_ff = 10;

% Number of spherical wave modes max 150 sonst zu instabil 
N = round(k0*r0 + 10);
maxidx = 2*N;

% Preallocate space for increased speed
Etheta_ff = zeros(size(theta_ff));
Ephi_ff = zeros(size(phi_ff));

idx = 1;
EthetaTest = zeros(size(phi_ff));
EphiTest = zeros(size(phi_ff));
for s=1:2
    for n = 1:N
        M = repmat((-n:n),length(theta),1);
        % Compute Spherical Expansion Coefficients
        [~,ftheta,fphi,Y] = sphericalVectorWaveFunction(s,M,n,A,theta,phi,k0); 
        ftheta_t = conj(ftheta).*sin(theta);
        fphi_t = conj(fphi).*sin(theta);

        q = (1./Y .*(Etheta.*ftheta_t + Ephi.*fphi_t))*delta_theta*delta_phi;  % NEU
        % q = (1./Y .*(Etheta'*ftheta_t + Ephi'*fphi_t))'*delta_theta*delta_phi;
        
        % Compute Far-Field 
        M = repmat((-n:n),length(theta_ff),1);
        [xtheta,xphi] = sphericalVectorWaveFunctionFarField(s,M,n,r_ff,theta_ff,phi_ff,k0);

        % Etheta_ff = Etheta_ff + sum(repmat(q,1,length(theta_ff))'.*xtheta,2);
        % Ephi_ff = Ephi_ff + sum(repmat(q,1,length(theta_ff))'.*xphi,2);

        Etheta_ff = EthetaTest + sum(q .* xtheta,2);  % NEU 
        Ephi_ff = EphiTest + sum (q .* xphi,2);   % NEU 

        disp ([num2str(idx), '/ ', num2str(maxidx)]);
        idx = idx+1;
        
    end
end


% Create Results Table
s = numel(Etheta_ff);
data_nf2ff = table(zeros(s,1),zeros(s,1),zeros(s,1),zeros(s,1),zeros(s,1));
data_nf2ff.Properties.VariableNames = {'theta','phi','Etheta','Ephi','Eabs'};


data_nf2ff.theta = reshape(theta_ff,s,1);
data_nf2ff.phi = reshape(phi_ff,s,1);
data_nf2ff.Etheta = reshape(Etheta_ff,s,1);
data_nf2ff.Ephi = reshape(Ephi_ff,s,1);
Eabs_nf_ff = sqrt(abs(data_nf2ff.Etheta).^2 + abs(data_nf2ff.Ephi).^2); % Gesamtbetrag des E-Felds
data_nf2ff.Eabs = reshape(Eabs_nf_ff,s,1);


%% Test E Absolut der Nf2ff Korrekt Beginn
% Test ob E absolut korrekt berechnet wurde mit der normalen Formel 
E_r = 0;
E_phi = data_nf2ff.Ephi;
E_theta = data_nf2ff.Etheta;



  Ex = E_r    .* sin(theta_ff) .* cos(phi_ff) ...
       + E_theta .* cos(theta_ff) .* cos(phi_ff) ...
       - E_phi   .* sin(phi_ff);

    % E_y
    Ey = E_r    .* sin(theta_ff) .* sin(phi_ff) ...
       + E_theta .* cos(theta_ff) .* sin(phi_ff) ...
       + E_phi   .* cos(phi_ff);

    % E_z
    Ez = E_r    .* cos(theta_ff) ...
       - E_theta .* sin(theta_ff);
Eabs_test = sqrt(abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2);
tol = 1e-3;
  if all(abs(Eabs_test - data_nf2ff.Eabs) < tol)
    disp('✅ Alle Werte sind gleich (mit Toleranz)/ Test E Absolut der Nf2ff / Variable E_abs');
    else
    disp('❌ Es gibt Unterschiede / Test E Absolut der Nf2ff / Variable E_abs');
    end
%% Test E Absolut der Nf2ff Korrekt Ende
end



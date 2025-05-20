function [] = plotNF2FF_phiCut(data_nf2ff,phi_cut,normalized,logarithmic)
% Diese Funktion erzeugt einen Phi-Cut (bei einem gewünschten phi in Grad (phi_cut)), wobei 
% alle zugehörigen Theta-Werte auf der x-Achse und die entsprechenden Eabs-Werte
% auf der y-Achse dargestellt werden. Die Darstellung von Eabs kann wahlweise 
% logarithmisch in dB (logarithmic = true) oder linear in V/m (logarithmic = false) 
% erfolgen. Zusätzlich kann der Wert entweder auf 0 dB normiert (normalize = true) 
% oder ungenormt (normalize = false) dargestellt werden.
% Wird hauptsächlich dafür genutzt die transformierten Daten zu plotten.
%
% Input Arguments:
%
%       data_nf2ff             Tabelle  mit mindestens 3 Spalten nötig:
%                              phi, theta, Eabs
%
%       phi_cut                gewünschter Winkel Phi(°) bei dem ein Schnitt
%                              gemacht werden soll
%
%       normalized            true für normierung auf 0 Db / false ohne normierung
%
%       logarithmic           true für logarithmisch in Db / false in V/m                    
%
%
% Output Arguments:
%
%       Keine                Es wird ein Plot enstehen ohne Rückgabe wert
%
deg = phi_cut; % für den Titel benötigt
phi_cut = deg2rad(phi_cut); % von Grad in Radian
phi_cut2 = phi_cut + pi; % auch von dem gleichen winkel + 180° 

% Zeilen suchen wo die gesuchten phi und phi+180 Grad sind
index = find(min(abs(data_nf2ff.phi - phi_cut))==abs(data_nf2ff.phi - phi_cut));
index2 = find(min(abs(data_nf2ff.phi - phi_cut2))==abs(data_nf2ff.phi - phi_cut2));

% Alle Theta für Phi und Phi+180°  in X-Achse -180° bis +180° bringen
theta_phiCut_nf2ff = [-fliplr(data_nf2ff.theta(index2)') data_nf2ff.theta(index)']; 

% Maximalen Wert ermitteln für normierung
maxValue1 = (max(data_nf2ff.Eabs(index)));
maxValue2 = (max(data_nf2ff.Eabs(index2)));

if maxValue1>maxValue2
    maxValue =maxValue1;
else
    maxValue=maxValue2;
end
% für die zugehöhrigen Winkeln den Eabs zuordnen
eabs_phiCut_nf2ff = [fliplr(data_nf2ff.Eabs(index2)')'; data_nf2ff.Eabs(index)];

% Normierung auf 0 Db falls normalized true
if normalized == true
    % test_maxVal = max(eabs_phiCut_nf2ff);
   eabs_phiCut_nf2ff = eabs_phiCut_nf2ff ./ maxValue;
end

% Logarithmische Darstellung falls logarithmic true
if logarithmic == true
    eabs_phiCut_nf2ff =  20 * log10(eabs_phiCut_nf2ff);
end

% Plot erzeugen 
figure (Name='Transformierte FF Daten, Phi Schnitt')
plot(rad2deg(theta_phiCut_nf2ff), eabs_phiCut_nf2ff, 'LineWidth', 1.5);
xlabel('Theta [°]');
if logarithmic == true
ylabel('|E| [dB]');
else 
ylabel('|E| [V/m]');  
end
title(['phi = ', num2str(deg) ,' Schnitt'])
end
%% NF2FF Conversion using Spherical Scanner
clear ; close all; clc;

disp('************************************************')
disp('          Near-To-Far-Field Conversion')
disp('          Spherical Scanner')
disp('************************************************')

addpath('misc_functions')
addpath('plot_functions')
addpath('transformation_functions')
addpath ('generator_functions')
%======================================================================
% Simulation von Nahfeld einer Antenne mithilfe der Antenna Toolbox
% und Vorbereitung der Daten für nf2ff_spherical.
%======================================================================

%% 1. Antenne definieren und Nahfeldatensimulation
% Falls gewünscht Parameter anpassen oder eine andere Antenne nehmen.
disp('Load Data...')
%% Hornantenne
% ant= horn;
% freq = 10e9;
% a1 = 4.87 * 0.0254;
% b1 = 3.62 * 0.0254;
% a = 0.9 * 0.0254;
% b = 0.4 * 0.0254;
% p = 10.06 * 0.0254;
% 
% ant.FlareLength =p;
% ant.FlareWidth = a1;
% ant.FlareHeight = b1;
% ant.Width = a;
% ant.Height = b;
% ant.Length =0.1;
% % ant.Tilt=270;
% ant.TiltAxis=[0 1 0];
% % Strahlendes Feld von 0,2 - 1,04 m

%% Dipol Antenne
ant = dipole;    
ant.Length = 0.5;
ant.Width = 0.005;
ant.Tilt=90;
ant.TiltAxis=[1 0 0];
freq = 3e8;

r_nf = 1.5; % z.B. 150 cm um die Antenne herum denn fernfeld bereits ab 2 m 
theta_nf = 0.1:1:179; % Theta darf NICHT 0 oder 180° sein, sont fehler, da durch sin(theta) geteilt wird
phi_nf   = 0.1:1:359; % Schrittweite 1° in Phi

%% 2 Simulation der Nahfelddaten
data_nf_raw = NearfieldDatagenerator(ant,freq,r_nf,theta_nf,phi_nf);


% Jetzt haben wir eine Tabelle mit den richtigen Spaltennamen für rearrangeTables
%% 3a "rearrangeTables" aufrufen 
data_nf_formatted = rearrangeTables(data_nf_raw);

%% 3b "rotateSphericalNFData" aufrufen
data_nf_rotated = rotateSphericalNFData(data_nf_formatted);


%% Jetzt die Daten and die Transformationsfunktion übergeben
% Wir definieren einen Abtastbereich im Fernfeld (theta_range, phi_range),
% dann rufen wir nf2ff_spherical_manual auf:
%% ACHTUNG: Theta darf nicht 0 oder Pi sein sonst kommts zu Nan Werten in der Transformation
theta_range = theta_nf.*pi/180; 
phi_range   = phi_nf.*pi/180;

data_nf2ff = nf2ff_spherical_manual(data_nf_rotated, freq, theta_range, phi_range);
disp('Done! Load Data')


%% 4a Plots der transformierten FF daten
disp('Plotting transformed FF Data...')


%% Theta und Phi Schnitt 
%% Variablen eingeben / Einheit Grad (°)
theta_cut = 90; 
phi_cut = 0;  
normalized = true; % true für normiert auf 0 Db, false für unnormiert 
logarithmic = true; % true für logarithmische, false für lineare Darstellung

%% Plots für Theta und Phi Schnitt
plotNF2FF_thetaCut (data_nf2ff,theta_cut,normalized,logarithmic);
plotNF2FF_phiCut (data_nf2ff,phi_cut,normalized,logarithmic);


%% Transformierten Daten in 3D Plotten 
E_dB = data_nf2ff.Eabs;
E_dB_norm =20* log10(E_dB./ max(E_dB));  % Max auf 0 dB normiert
phi = rad2deg (data_nf2ff.phi);
theta = rad2deg (data_nf2ff.theta);
figure (Name='Transformierte FF Daten')
patternCustom(E_dB_norm,theta, phi);
title('Transformiertes Fernfeldstrahlungsdiagramm');
colorbar;                          % stellt sicher, dass eine Farbskala angezeigt wird
ylabel(colorbar, 'Amplitude [dB]');

%% 4b Nahfeldaten in 3D Plotten 
E_db_2 = data_nf_rotated.Eabs;
E_dB_norm_2 = 20 * log10 (E_db_2./max(E_db_2));
phi_2 = rad2deg (data_nf_rotated.phi);
theta_2 = rad2deg(data_nf_rotated.theta);
figure (Name='Matlab simulierte NF Daten')
patternCustom(E_dB_norm_2,theta,phi);
title('Simuliertes Nahfeldstrahlungsdiagramm');
colorbar;                         
ylabel(colorbar, 'Amplitude [dB]'); 

%% 4c Matlab Simulation der Fernfelddaten
Fernfelddaten_Simulation (ant,freq,theta_cut,phi_cut);
disp('Done! Plotting transformed FF Data')
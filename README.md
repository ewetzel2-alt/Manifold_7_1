clc; clear;

% Variables
freq = 1.57542e9;
% Corrected AZ range for the pattern function
az = -180:1:180; 
el = 0:1:90;
c = physconst('LightSpeed'); % ~3e8

patch_height = 0.0016; % 1.6 mm height
d = dielectric('Name', 'FR4', 'EpsilonR', 4.4, 'LossTangent', 0.02);


% --- 8-Element CRPA Configuration (Adjusted for Spacing) ---
N_ring = 7; 
% Increase radius to ~0.11m to achieve closer to 0.5*lambda spacing
arrayRadius = 0.107; 
angles = linspace(0, 2*pi, N_ring + 1);
angles = angles(1:N_ring);

% Calculate positions of the 7 ring elements (in meters)
x_ring = arrayRadius * cos(angles);
y_ring = arrayRadius * sin(angles);
z_ring = zeros(1, N_ring);

% Define position of the 1 center element (at the origin)
x_center = 0;
y_center = 0;
z_center = 0;

% Combine all element positions into a single 3x8 matrix (X; Y; Z)
elementPositions = [
    [x_ring, x_center];
    [y_ring, y_center];
    [z_ring, z_center]
];

% Create the conformal array object (requires Phased Array System Toolbox)
crpa_array = phased.ConformalArray('ElementPosition', elementPositions);


% Optional: Assign your specific patch antenna 'p' to all elements
% crpa_array.Element = p; 

% --- Add visualization here ---
% figure;
% viewArray(crpa_array, 'Title', '3D Visualization of 8-Element CRPA (7+1 Config)');

figure;
% This plots the default pattern
pattern(crpa_array, freq, az, el);
title(['3D Radiation Pattern at ' num2str(freq/1e9) ' GHz (Default Weights)']);



%Steer the antenna beam (with math)
c = physconst('LightSpeed'); % Get speed of light
lambda = c / freq;
k = 2 * pi / lambda; % Wavenumber

% 2. Define the desired steering direction [Azimuth; Elevation] in degrees
steer_az = 0; 
steer_el = 90;

% Convert the desired direction from spherical (az, el) to a Cartesian unit vector (u_steer)
[x_u, y_u, z_u] = sph2cart(deg2rad(steer_az), deg2rad(steer_el), 1);
u_steer = [x_u; y_u; z_u]; % 3x1 unit vector pointing in the steering direction

% 3. Calculate Path Differences (Delta L) and Phase Shifts (Delta Phi)
% Delta L = dot product of element position vector and the steering unit vector
% elementPositions is 3xN, u_steer is 3x1.
% The result 'path_diffs' is a 1xN vector of Delta L_n values
path_diffs = dot(elementPositions, repmat(u_steer, 1, size(elementPositions, 2)));

% Calculate Phase Shifts: Delta Phi = k * Delta L
phase_shifts_rad = k * path_diffs;

% 4. Create the Complex Weights Vector 'w' (w_n = exp(j * Delta Phi_n))
% The negative sign is crucial to *compensate* for the phase delay
w = exp(-1j * phase_shifts_rad)'; 

% --- Add visualization for the steered pattern ---

figure;
% Pass your manually calculated weights 'w' to the pattern function
pattern(crpa_array, freq, az, el, 'Type', 'directivity', 'Weights', w);
title(['3D Radiation Pattern Manual Calc at ' num2str(freq/1e9) ' GHz, Steered to Az ' num2str(steer_az) ' El ' num2str(steer_el)]);



%Impedance
% --- Conversion for Impedance Analysis (Requires Antenna Toolbox) ---
% 1. Create a physical element (e.g., a dipole or your patch 'p')
resonance_factor = 0.511; 
elem = dipole('Length', lambda * resonance_factor, 'Width', lambda/100); 

% 2. Re-run physical array and impedance
phys_array = conformalArray('Element', elem, 'ElementPosition', elementPositions');
z_matrix = impedance(phys_array, freq); 

% 2. Create an Antenna Toolbox 'conformalArray' using your existing positions
% Transpose elementPositions to get N-by-3
phys_array = conformalArray('Element', elem, ...
                            'ElementPosition', elementPositions');

% 3. Now the impedance function will work correctly
z_matrix = impedance(phys_array, freq); 

% --- Quantifying Mutual Coupling's Impact (New Section) ---

% 'w' is your weights vector from the manual steering section above
% 'z_matrix' is from the Antenna Toolbox section above
% 'Z0' is the desired system impedance (usually 50 Ohms)
Z0 = 50; 

% Calculate the Active Impedance for your specific steered beam
% Formula: Z_active = (Z_matrix * V_terminal) ./ I_terminal
% In MATLAB terms, Z_active = (Z_matrix * w) ./ w 
Z_active = (z_matrix * w) ./ w;

% Calculate the Reflection Coefficient (Gamma)
Gamma = (Z_active - Z0) ./ (Z_active + Z0);

% Calculate the VSWR
VSWR = (1 + abs(Gamma)) ./ (1 - abs(Gamma));

fprintf('\n--- Mutual Coupling Analysis (Steered to Az 45, El 60) ---\n');
fprintf('Center Element (8) Active Impedance: %.2f + j%.2f Ohms\n', real(Z_active(8)), imag(Z_active(8)));
fprintf('Ring Element (1) Active Impedance:    %.2f + j%.2f Ohms\n', real(Z_active(1)), imag(Z_active(1)));
fprintf('Center Element (8) VSWR: %.2f\n', VSWR(8));
fprintf('Ring Element (1) VSWR:    %.2f\n', VSWR(1));

% --- 3D Pattern Visualization ---
figure;
% Display the pattern with the complex weights 'w'
pattern(crpa_array, freq, az, el, 'Type', 'directivity', 'Weights', w);
title(['Optimized 3D Pattern (VSWR ~1.8-2.4) Steered to Az ' ...
      num2str(steer_az) ', El ' num2str(steer_el)]);
colormap('jet'); % High-contrast coloring for clarity


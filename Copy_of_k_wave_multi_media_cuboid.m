clear
%*******************************************************May need update thes below parameters
load('grid_data.mat', 'begin_x_grid', 'begin_y_grid', 'begin_z_grid', 'i','j','ii','num_x_grid','num_y_grid','num_z_grid','x_grid_first', 'y_grid_first', 'z_grid_first');       


% create the computational grid
Nx = num_x_grid;            % number of grid points in the x direction
Ny = num_y_grid;            % number of grid points in the y direction
Nz = num_z_grid;            % number of grid points in the z direction
wavelength = 346 * 25/1000000;  
incresed_pressure = 0;
temper = 18;
fprintf(' : %d Temperature\n', temper);
kk = ii;
% create initial sound source to emit wave
target_x = begin_x_grid; % grid index
target_y = begin_y_grid;
target_z = begin_z_grid;
xyz_first_grid = [x_grid_first, y_grid_first, z_grid_first]/1000; % unit is m

grid_size = wavelength / 12;  %  0.0007208 mm = 0.7208 mm
dx = grid_size; % grid point spacing in the x direction [m]
dy = dx;             % grid point spacing in the y direction [m]
dz = dx;             % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
% load the grid index of ellipse solid and array
ellipse_index = readmatrix('inside_cuboid_indices_3d.csv');
%********************************************************************************************


PML_X_SIZE = 10;PML_Y_SIZE = 10;PML_Z_SIZE = 10; %PML grid points number
% define the properties of the propagation medium
heat = temper + 273.15 ;% trasfer degree to K
R_air = 287.06; % refer to http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/fprops/propsoffluids/node12.html
R_co2 = 188.9; % 
pressure_air = 102000;% unit is Pascal, 101400 is from HK Observatory data
pressure_balloon = pressure_air + incresed_pressure; 
gamma = 1.4;
gamma_co2 = 1.3;


density_balloon = 1.491; % for mixed gases
density_air = pressure_air/R_air/heat; % standard density is 1.1839 kg/m^3 @ 25 degrees.
% density_balloon = density_air;

speed_air = sqrt(gamma*pressure_air/density_air);
% speed_balloon = speed_air;
speed_balloon = 303; % for mixed gases



medium.sound_speed = speed_air * ones(Nx, Ny, Nz);	% [m/s] in air
medium.density = density_air * ones(Nx, Ny, Nz);    % [kg/m^3]
for iii = 1:size(ellipse_index, 1)
    x_idx = ellipse_index(iii, 1); 
    y_idx = ellipse_index(iii, 2); 
    z_idx = ellipse_index(iii, 3); 
    medium.sound_speed(x_idx, y_idx, z_idx) = speed_balloon;
    medium.density(x_idx, y_idx, z_idx) = density_balloon;
end
source_radius = 1; % [grid points]
source.p_mask = zeros(Nx, Ny, Nz);
source.p_mask(target_x, target_y, target_z) = 1; % grid index, with the same dimensions as the computational grid


element_index = readmatrix('element_grid_indices.csv'); % calculate the distance between source point to the four cornor points on array plane.
first_grid_of_array = [1*dx,  element_index(1,2)*dy, 1*dz]; 
secon_grid_of_array = [1*dx,  element_index(1,2)*dy, Nx*dz];
third_grid_of_array = [Nx*dx, element_index(1,2)*dy, 1*dz]; 
fourt_grid_of_array = [Nx*dx, element_index(1,2)*dy, Nx*dz];
dis_grid_first = sqrt((target_x*dx-first_grid_of_array(1))^2+(target_y*dy-first_grid_of_array(2))^2+(target_z*dz-first_grid_of_array(3))^2);
dis_grid_secon = sqrt((target_x*dx-secon_grid_of_array(1))^2+(target_y*dy-secon_grid_of_array(2))^2+(target_z*dz-secon_grid_of_array(3))^2);
dis_grid_third = sqrt((target_x*dx-third_grid_of_array(1))^2+(target_y*dy-third_grid_of_array(2))^2+(target_z*dz-third_grid_of_array(3))^2);
dis_grid_fourt = sqrt((target_x*dx-fourt_grid_of_array(1))^2+(target_y*dy-fourt_grid_of_array(2))^2+(target_z*dz-fourt_grid_of_array(3))^2);
max_dis_grid = max([dis_grid_first,dis_grid_secon,dis_grid_third,dis_grid_fourt]); % [m]

source_freq = 40e3;  % [Hz]
source_mag = 10;     % [Pa]
T = 1/source_freq; % [s]
t_end = max_dis_grid/speed_air + 6*T; % simulation time, [s]
kgrid.makeTime(medium.sound_speed,[],t_end);
% define a time varying sinusoidal source
source_cycle=0.5;
pressure_emission=10;
source.p=pressure_emission*toneBurst(1/kgrid.dt,source_freq,source_cycle); % source
% figure (1);
% plot(0:kgrid.dt:0.5/source_freq, source.p);
% title('Signal of sound source');

% define a series of probes (196 probes) for signal achieved
x_element = element_index(:, 1)';  
y_element = element_index(:, 2)';
z_element = element_index(:, 3)';

sensor.mask = zeros(Nx, Ny, Nz);
%       sensor.mask can be defined in three different ways. (1) As a binary matrix (i.e., a matrix of 1's
%       and 0's with the same dimensions as the computational grid)
%       representing the grid points within the computational grid that
%       will collect the data.
for iii=1:1:196 % refer to https://blog.csdn.net/chen_guowei2000/article/details/136309258
    sensor.mask(x_element(1,iii), y_element(1,iii), z_element(1,iii)) = 1; % grid index
end


figure (2);
hold on;
plot3(x_element, y_element, z_element, 'bo', 'MarkerSize', 5); % index
xlabel('X Index');
ylabel('Y Index');
zlabel('Z Index');
title('196 Probe Positions in 3D Grid');
grid on;
axis equal;
xlim([1, Nx]);
ylim([1, Ny]);
zlim([1, Nz]);
hold off;
% input arguments
input_args = {'PlotLayout', true, 'DataCast', 'single', 'PlotSim',false};


figure(3);
hold on;
axis equal;
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
title('3D View of 196 Circles');
grid on;
for a=1:1:196
    x_3d = (x_element(1,a)- 1)* dx*1000 + xyz_first_grid(1)*1000; % [mm]
    y_3d = (y_element(1,a)- 1)* dy*1000 + xyz_first_grid(2)*1000;
    z_3d = (z_element(1,a)- 1)* dz*1000 + xyz_first_grid(3)*1000; 
  
    theta = linspace(0, 2 * pi, 100);  % 
    x_circle = 4.95 * cos(theta) + x_3d;
    y_circle = y_3d * ones(size(theta));
    z_circle = 4.95 * sin(theta) + z_3d;  %

    plot3(x_circle, y_circle, z_circle, 'b', 'LineWidth', 1.5);
end

xlim([xyz_first_grid(1)*1000, xyz_first_grid(1)*1000 + Nx * dx*1000]);
ylim([xyz_first_grid(2)*1000, xyz_first_grid(2)*1000 + Ny * dy*1000]);
zlim([xyz_first_grid(3)*1000, xyz_first_grid(3)*1000 + Nz * dy*1000]); %[mm]
hold off;

% run the simulation, 
% sensor_data is NxM, N is number of sensors, M is time step
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}); 


sensor_data1 = zeros;
col=0; row=0;
for l=1:1:196
    if(mod(l,14)==0)
        col = 14;
        row = l/14-1;
    else
        col = mod(l,14);
        row = (l-col)/14;
    end
    sensor_data1((col-1)*14+14-row,1:length(sensor_data)) = sensor_data(l,:);
end

xyz_3D_start_x = (target_x - 1) * dx *1000 + xyz_first_grid(1)*1000; % [mm]
xyz_3D_start_y = (target_y - 1) * dy *1000 + xyz_first_grid(2)*1000; % [mm]
xyz_3D_start_z = (target_z - 1) * dz *1000 + xyz_first_grid(3)*1000; % [mm]

iw = 1;
xyz_3D_e1_x = (x_element(1,iw)- 1)* dx*1000 + xyz_first_grid(1)*1000; % [mm]
xyz_3D_e1_y = (y_element(1,iw)- 1)* dy*1000 + xyz_first_grid(2)*1000;
xyz_3D_e1_z = (z_element(1,iw)- 1)* dz*1000 + xyz_first_grid(3)*1000;
dis_e1 = sqrt((xyz_3D_e1_x-xyz_3D_start_x)^2 + (xyz_3D_e1_y-xyz_3D_start_y)^2 + (xyz_3D_e1_z-xyz_3D_start_z)^2);
tof_e1 = dis_e1/speed_air;

ia = 28;
xyz_3D_e28_x = (x_element(1,ia)- 1)* dx*1000 + xyz_first_grid(1)*1000; % [mm]
xyz_3D_e28_y = (y_element(1,ia)- 1)* dy*1000 + xyz_first_grid(2)*1000;
xyz_3D_e28_z = (z_element(1,ia)- 1)* dz*1000 + xyz_first_grid(3)*1000;
dis_e28 = sqrt((xyz_3D_e28_x-xyz_3D_start_x)^2 + (xyz_3D_e28_y-xyz_3D_start_y)^2 + (xyz_3D_e28_z-xyz_3D_start_z)^2);
tof_e28 = dis_e28/speed_air;

ia = 34;
xyz_3D_e34_x = (x_element(1,ia)- 1)* dx*1000 + xyz_first_grid(1)*1000; % [mm]
xyz_3D_e34_y = (y_element(1,ia)- 1)* dy*1000 + xyz_first_grid(2)*1000;
xyz_3D_e34_z = (z_element(1,ia)- 1)* dz*1000 + xyz_first_grid(3)*1000;
dis_e34 = sqrt((xyz_3D_e34_x-xyz_3D_start_x)^2 + (xyz_3D_e34_y-xyz_3D_start_y)^2 + (xyz_3D_e34_z-xyz_3D_start_z)^2);
tof_e34 = dis_e34/speed_air;




% i = 99;
% xyz_3D_e99_x = (x_element(1,i)- 1)* dx*1000 + xyz_first_grid(1)*1000; % [mm]
% xyz_3D_e99_y = (y_element(1,i)- 1)* dy*1000 + xyz_first_grid(2)*1000;
% xyz_3D_e99_z = (z_element(1,i)- 1)* dz*1000 + xyz_first_grid(3)*1000;
% dis_e99 = sqrt((xyz_3D_e99_x-xyz_3D_start_x)^2 + (xyz_3D_e99_y-xyz_3D_start_y)^2 + (xyz_3D_e99_z-xyz_3D_start_z)^2);
% tof_e99 = dis_e99/speed_air;

% i = 98;
% xyz_3D_e98_x = (x_element(1,i)- 1)* dx*1000 + xyz_first_grid(1)*1000; % [mm]
% xyz_3D_e98_y = (y_element(1,i)- 1)* dy*1000 + xyz_first_grid(2)*1000;
% xyz_3D_e98_z = (z_element(1,i)- 1)* dz*1000 + xyz_first_grid(3)*1000;
% dis_e98 = sqrt((xyz_3D_e98_x-xyz_3D_start_x)^2 + (xyz_3D_e98_y-xyz_3D_start_y)^2 + (xyz_3D_e98_z-xyz_3D_start_z)^2);
% tof_e98 = dis_e98/speed_air;

% i = 92;
% xyz_3D_e92_x = (x_element(1,i)- 1)* dx*1000 + xyz_first_grid(1)*1000; % [mm]
% xyz_3D_e92_y = (y_element(1,i)- 1)* dy*1000 + xyz_first_grid(2)*1000;
% xyz_3D_e92_z = (z_element(1,i)- 1)* dz*1000 + xyz_first_grid(3)*1000;
% dis_e92 = sqrt((xyz_3D_e92_x-xyz_3D_start_x)^2 + (xyz_3D_e92_y-xyz_3D_start_y)^2 + (xyz_3D_e92_z-xyz_3D_start_z)^2);
% tof_e92 = dis_e92/speed_air;


%visuliazation of the inition pressure and the sensor mask
%voxelPlot(double(source.p_mask | cart2grid(kgrid, sensor.mask)));% cart2grid(kgrid, cart_data),
%cart_data is 1xN, 2xN or 3xN.

% figure (6);
% plot(0:kgrid.dt:kgrid.dt*(length(sensor_data(1,:))-1),sensor_data(1,:));
% for i=28:1:28
%     hold on
%     plot(0:kgrid.dt:kgrid.dt*(length(sensor_data(1,:))-1),sensor_data(i,:));
% end

% hold off
dtt = kgrid.dt;
save('kgrid_data.mat','dtt','t_end','kk','i','j');
writematrix(sensor_data1, 'sensor_data.csv');

% Time reverse based 3D reconstruction 
% refer to http://www.k-wave.org/documentation/example_pr_3D_tr_planar_sensor.php


% source.p = 0;i
% sensor.time_reversal_boundary_data = sensor_data; % assign the time reversal data
% p0_recon = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
% figure (6);
% % beamPlot(abs(p0_recon));

% % colorbar
% hold on;
% subplot(1,2,1);
% slice_xz = squeeze(abs(p0_recon(:, target_y, :)));
% imagesc(slice_xz);
% xlabel('Z axis');
% ylabel('X axis');
% subplot(1,2,2);
% slice_xy = squeeze(abs(p0_recon(:, :, target_z)));
% imagesc(slice_xy);
% xlabel('Y axis');
% ylabel('X axis');
% hold off






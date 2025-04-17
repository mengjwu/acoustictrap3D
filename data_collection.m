

clear
% physical_coordinate，unit is mm, x is horizonation y is vertical
diary('log-mjwu.txt') ;

speed_air = 346;
wavelength = 346 * 25/1000;  
grid_size = wavelength / 12;  

element_xyz = readmatrix('element_xyz.csv'); % from "physical_coordinate.m"

% cuboid CO2 size
length_cuboid = 210-9.5*2-30*2; 
width_cuboid = 193-4.5*2-9.5-30;   
height_cuboid = 125-9.5-30; 
fprintf('CO2 Volume: %f mm3，can fill %f 100 mL \n', length_cuboid*width_cuboid *height_cuboid,length_cuboid*width_cuboid *height_cuboid/100000);


z_14 = element_xyz(14,3) - 10; % Z dirction
y_14 = element_xyz(14,2) - 17; 
x_14 = element_xyz(14,1) + 5; 
cuboid_bott_left = [x_14, y_14, z_14];

range_total_cuboid = readmatrix('range_total_cuboid.csv'); % from "physical_coordinate.m"

x_range = [range_total_cuboid(1,1), range_total_cuboid(1,2)];  % x, mm
y_range = [range_total_cuboid(2,1), range_total_cuboid(2,2)];  % y
z_range = [range_total_cuboid(3,1), range_total_cuboid(3,2)]; % z, mm


v1 = cuboid_bott_left;         
v2 = v1 - [0, height_cuboid, 0]; 
v3 = v1 + [length_cuboid, 0, 0]; 
v4 = v2 + [length_cuboid, 0, 0]; 
v5 = v1 + [0, 0, width_cuboid]; 
v6 = v2 + [0, 0, width_cuboid]; 
v7 = v3 + [0, 0, width_cuboid];  
v8 = v4 + [0, 0, width_cuboid];  


figure;
hold on;
axis equal;
xlabel('X [mm]');
ylabel('Y [mm]');
zlabel('Z [mm]');
title('Cuboid with Outer Frame and Internal Volume');


cube_vertices = [x_range(1), y_range(1), z_range(1);  
                 x_range(1), y_range(1), z_range(2);
                 x_range(1), y_range(2), z_range(1);
                 x_range(1), y_range(2), z_range(2);
                 x_range(2), y_range(1), z_range(1);
                 x_range(2), y_range(1), z_range(2);
                 x_range(2), y_range(2), z_range(1);
                 x_range(2), y_range(2), z_range(2)];


line([cube_vertices(1,1), cube_vertices(5,1)], [cube_vertices(1,2), cube_vertices(5,2)], [cube_vertices(1,3), cube_vertices(5,3)], 'Color', 'k', 'LineWidth', 2);
line([cube_vertices(2,1), cube_vertices(6,1)], [cube_vertices(2,2), cube_vertices(6,2)], [cube_vertices(2,3), cube_vertices(6,3)], 'Color', 'k', 'LineWidth', 2);
line([cube_vertices(3,1), cube_vertices(7,1)], [cube_vertices(3,2), cube_vertices(7,2)], [cube_vertices(3,3), cube_vertices(7,3)], 'Color', 'k', 'LineWidth', 2);
line([cube_vertices(4,1), cube_vertices(8,1)], [cube_vertices(4,2), cube_vertices(8,2)], [cube_vertices(4,3), cube_vertices(8,3)], 'Color', 'k', 'LineWidth', 2);
line([cube_vertices(1,1), cube_vertices(3,1)], [cube_vertices(1,2), cube_vertices(3,2)], [cube_vertices(1,3), cube_vertices(3,3)], 'Color', 'k', 'LineWidth', 2);
line([cube_vertices(2,1), cube_vertices(4,1)], [cube_vertices(2,2), cube_vertices(4,2)], [cube_vertices(2,3), cube_vertices(4,3)], 'Color', 'k', 'LineWidth', 2);
line([cube_vertices(5,1), cube_vertices(7,1)], [cube_vertices(5,2), cube_vertices(7,2)], [cube_vertices(5,3), cube_vertices(7,3)], 'Color', 'k', 'LineWidth', 2);
line([cube_vertices(6,1), cube_vertices(8,1)], [cube_vertices(6,2), cube_vertices(8,2)], [cube_vertices(6,3), cube_vertices(8,3)], 'Color', 'k', 'LineWidth', 2);
line([cube_vertices(1,1), cube_vertices(2,1)], [cube_vertices(1,2), cube_vertices(2,2)], [cube_vertices(1,3), cube_vertices(2,3)], 'Color', 'k', 'LineWidth', 2);
line([cube_vertices(3,1), cube_vertices(4,1)], [cube_vertices(3,2), cube_vertices(4,2)], [cube_vertices(3,3), cube_vertices(4,3)], 'Color', 'k', 'LineWidth', 2);
line([cube_vertices(5,1), cube_vertices(6,1)], [cube_vertices(5,2), cube_vertices(6,2)], [cube_vertices(5,3), cube_vertices(6,3)], 'Color', 'k', 'LineWidth', 2);
line([cube_vertices(7,1), cube_vertices(8,1)], [cube_vertices(7,2), cube_vertices(8,2)], [cube_vertices(7,3), cube_vertices(8,3)], 'Color', 'k', 'LineWidth', 2);

fill3([v1(1), v3(1), v7(1), v5(1)], [v1(2), v3(2), v7(2), v5(2)], [v1(3), v3(3), v7(3), v5(3)], 'r', 'FaceAlpha', 0.3); % 前面
fill3([v2(1), v4(1), v8(1), v6(1)], [v2(2), v4(2), v8(2), v6(2)], [v2(3), v4(3), v8(3), v6(3)], 'r', 'FaceAlpha', 0.3); % 后面
fill3([v1(1), v2(1), v6(1), v5(1)], [v1(2), v2(2), v6(2), v5(2)], [v1(3), v2(3), v6(3), v5(3)], 'r', 'FaceAlpha', 0.3); % 左面
fill3([v3(1), v4(1), v8(1), v7(1)], [v3(2), v4(2), v8(2), v7(2)], [v3(3), v4(3), v8(3), v7(3)], 'r', 'FaceAlpha', 0.3); % 右面
fill3([v1(1), v2(1), v4(1), v3(1)], [v1(2), v2(2), v4(2), v3(2)], [v1(3), v2(3), v4(3), v3(3)], 'r', 'FaceAlpha', 0.3); % 顶面
fill3([v5(1), v6(1), v8(1), v7(1)], [v5(2), v6(2), v8(2), v7(2)], [v5(3), v6(3), v8(3), v7(3)], 'r', 'FaceAlpha', 0.3); % 底面


num_elements = size(element_xyz, 1); 
radius = 5; % [mm]
theta = linspace(0, 2*pi, 100); 

for i = 1:num_elements
  
    center_x = element_xyz(i, 1);
    center_y = element_xyz(i, 2);
    center_z = element_xyz(i, 3);
    
   
    circle_x = center_x + radius * cos(theta);
    circle_z = center_z + radius * sin(theta);
    
   
    fill3(circle_x, center_y * ones(size(circle_x)), circle_z, 'k', 'EdgeColor', 'none'); 
end
hold off;


x_grid = adjust_range(x_range(1), x_range(2), grid_size); % mm
y_grid = adjust_range(y_range(1), y_range(2), grid_size);
z_grid = adjust_range(z_range(1), z_range(2), grid_size);
x_grid_first = x_grid(1);
y_grid_first = y_grid(1);
z_grid_first = z_grid(1);

num_x_grid = length(x_grid);
num_y_grid = length(y_grid);
num_z_grid = length(z_grid);


[X, Y, Z] = meshgrid(x_grid, y_grid, z_grid);

grid_coordinates = [X(:), Y(:), Z(:)];% the coordinate value in the physcial world
xyz_first_grid=[x_grid(1), y_grid(1), z_grid(1)];

element_grid_indices = zeros;
for idx = 1:num_elements
    element_pos = element_xyz(idx, :); % [x, y, z]
    grid_index_x = round((element_pos(1) - xyz_first_grid(1)) / grid_size) + 1;%[mm] to index
    grid_index_y = round((element_pos(2) - xyz_first_grid(2)) / grid_size) + 1;
    grid_index_z = round((element_pos(3) - xyz_first_grid(3)) / grid_size) + 1;
    element_grid_indices(idx, 1:3) = [grid_index_x, grid_index_y, grid_index_z];
end

writematrix(element_grid_indices, 'element_grid_indices.csv');


for i=3:1:3 %4mm
    for j=6:1:6 %for j
    fprintf('*********************************************************************************************************************************\n');
    element_xyz = readmatrix('element_xyz.csv'); % from "physical_coordinate.m"
    length_cuboid = 131;
    width_cuboid = 144.5;
    height_cuboid = 85.5;
    grid_size = 0.7208;

    
    z_14 = element_xyz(14,3) - 10; % Z dirction
    y_14 = element_xyz(14,2) - 16 - i;  
    x_14 = element_xyz(14,1) - 1 + j;   
    cuboid_bott_left = [x_14, y_14, z_14];
    
 
    v1 = cuboid_bott_left;            
    v2 = v1 - [0, height_cuboid, 0];  
    v3 = v1 + [length_cuboid, 0, 0];  
    v4 = v2 + [length_cuboid, 0, 0];  
    v5 = v1 + [0, 0, width_cuboid];  
    v6 = v2 + [0, 0, width_cuboid];  
    v7 = v3 + [0, 0, width_cuboid];   
    v8 = v4 + [0, 0, width_cuboid];  

 
    x_min_cuboid = v1(1,1);  
    x_max_cuboid = v3(1,1); 
    y_min_cuboid = v2(1,2); 
    y_max_cuboid = v1(1,2);  
    z_min_cuboid = v1(1,3);  
    z_max_cuboid = v5(1,3);  
    xyz_first_grid=[-52.199, -59.883583, 156.1];
     
    grid_min_x = round((x_min_cuboid - xyz_first_grid(1)) / grid_size) + 1; % [mm] to index
    grid_max_x = round((x_max_cuboid - xyz_first_grid(1)) / grid_size) + 1;

    grid_min_y = round((y_min_cuboid - xyz_first_grid(2)) / grid_size) + 1; % [mm] to index
    grid_max_y = round((y_max_cuboid - xyz_first_grid(2)) / grid_size) + 1;

    grid_min_z = round((z_min_cuboid - xyz_first_grid(3)) / grid_size) + 1; % [mm] to index
    grid_max_z = round((z_max_cuboid - xyz_first_grid(3)) / grid_size) + 1;

    number_voexl = (grid_max_x-grid_min_x+1)*(grid_max_y-grid_min_y+1)*(grid_max_z-grid_min_z+1);
 
    inside_cuboid_indices_3d = zeros(number_voexl,3);
 
    counter = 1;
    for x_idx = grid_min_x:grid_max_x
        % fprintf('x_idex: %d\n', x_idx);
      for y_idx = grid_min_y:grid_max_y
        % fprintf('y_idex: %d\n', y_idx);
            for z_idx = grid_min_z:grid_max_z
              inside_cuboid_indices_3d(counter, 1:3) = [x_idx, y_idx, z_idx];
              counter = counter + 1;
            end
      end
    end
 
    writematrix(inside_cuboid_indices_3d, 'inside_cuboid_indices_3d.csv');

    for ii=33:1:63  %% vertical direction: 34mm to 94 mm

        xyz_first_grid=[-52.199, -59.883583, 156.1];
        num_x_grid = 217;
        num_y_grid = 175;
        num_z_grid = 231;
        x_grid_first = xyz_first_grid(1);
        y_grid_first = xyz_first_grid(2);
        z_grid_first = xyz_first_grid(3);

 
        element_xyz = readmatrix('element_xyz.csv');
        grid_size = 0.7208;% [mm]
        fprintf('=============================================================================================\n');
        fprintf('KK=%d,\n', ii);
        if (mod(ii,7)==0) % 7
            begin_x = 0.5 * element_xyz(85,1) + 0.5 * element_xyz(99,1) -3*1.5 + (7-1)*1.5; %  
            begin_y = element_xyz(85,2) - 34-(floor((ii)/7)-1)*1.5; %  
        elseif(mod(ii,7)) % 1,2,3,4,5,6
            begin_x = 0.5 * element_xyz(85,1) + 0.5 * element_xyz(99,1) -3*1.5 + (mod(ii,7)-1)*1.5;  
            begin_y = element_xyz(85,2) - 34-floor((ii)/7)*1.5; % 正上方34mm处 start
        end
        begin_z = 0.5 * element_xyz(91,3) + 0.5 *element_xyz(92,3);% center
        
        begin_x_grid =  round((begin_x - xyz_first_grid(1)) / grid_size) + 1;
        begin_y_grid =  round((begin_y - xyz_first_grid(2)) / grid_size) + 1;
        begin_z_grid =  round((begin_z - xyz_first_grid(3)) / grid_size) + 1;
        
        save('grid_data.mat', 'begin_x_grid', 'begin_y_grid', 'begin_z_grid','i','j','ii','num_x_grid','num_y_grid','num_z_grid','x_grid_first', 'y_grid_first', 'z_grid_first');
        
        run('Copy_of_k_wave_multi_media_cuboid.m');
        run('extract_maxpoint.m');
    end
    end % for j
end

diary off;
 
function adjusted_range = adjust_range(range_start, range_end, grid_size)
    total_range = range_end - range_start;
    remainder = mod(total_range, grid_size);
    if remainder >= (grid_size / 2)
        range_end = range_end + (grid_size - remainder); 
    else
        range_end = range_end - remainder;   
    end
    adjusted_range = (range_start + grid_size/2) : grid_size: range_end;
end

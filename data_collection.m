% --我们先把思路写在这，再编程序
% --首先我们根据PCB上的那个mark点，来定义大空气立方体的范围，然后再定义co2的范围，最后改变目标点位置

clear
% physical_coordinate，unit is mm, x is horizonation 往右为正, y is vertical 往下为正
diary('log-mjwu.txt') ;
%*********************************************每次计算需要更新，保持上下计算一致******************************************
speed_air = 346;
wavelength = 346 * 25/1000;  % air 波长 8.65 mm,% 网格大小，波长的1/12, [mm]
grid_size = wavelength / 12;  % 网格大小 0.7208 [mm]
% 读取MAT文件的element数据
element_xyz = readmatrix('element_xyz.csv'); % from "physical_coordinate.m"

% cuboid CO2 的已知实际尺寸
length_cuboid = 210-9.5*2-30*2; % 水平方向长度 (x方向) [mm], 9.5是蓝色海绵的厚度，30是黄色海绵厚度
width_cuboid = 193-4.5*2-9.5-30;  % 深度 (z方向) [mm]
height_cuboid = 125-9.5-30; % 竖直方向高度 (y方向) [mm]
fprintf('CO2的空间体积是: %f mm3，相当于%f个100 mL注射器\n', length_cuboid*width_cuboid *height_cuboid,length_cuboid*width_cuboid *height_cuboid/100000);

% 计算 cuboid 左下角顶点位置 (通过 point1 修正)
% 求解第14个element的位置，然后利用这个信息，求解左上角顶点 [x, y, z]，从而得到整个co2立方体的位置
z_14 = element_xyz(14,3) - 10; % Z dirction
y_14 = element_xyz(14,2) - 17; %当盒子竖直移动时候，初始态y值，初始值17，后面箱子会变
x_14 = element_xyz(14,1) + 5;  %当盒子水平移动时候，x值变化
cuboid_bott_left = [x_14, y_14, z_14];
fprintf('亚克力盒子内部有效腔体的 左下角 相对于element 14偏移(mm)： %f，%f， %f，初始态(mm)是5，17，10. \n', x_14-element_xyz(14,1), element_xyz(14,2)-y_14, element_xyz(14,3)-z_14);

%读取整个大立方体范围
range_total_cuboid = readmatrix('range_total_cuboid.csv'); % from "physical_coordinate.m"
%********************************************************每次计算需要更新，保持上下计算一致******************************************

%首先做大立方体范围定义
x_range = [range_total_cuboid(1,1), range_total_cuboid(1,2)];  % x方向, mm
y_range = [range_total_cuboid(2,1), range_total_cuboid(2,2)];  % y方向
z_range = [range_total_cuboid(3,1), range_total_cuboid(3,2)]; % z方向, mm

%首先是Co2腔体没移动时候，我们把箱子初始态和array，以及整个大立方体范围画出来，这样确认下，它是否正确
v1 = cuboid_bott_left;           % 左下点
v2 = v1 - [0, height_cuboid, 0]; % 左上点
v3 = v1 + [length_cuboid, 0, 0]; % 右下点
v4 = v2 + [length_cuboid, 0, 0]; % 右上点
v5 = v1 + [0, 0, width_cuboid];  % 后左下
v6 = v2 + [0, 0, width_cuboid];  % 后左上
v7 = v3 + [0, 0, width_cuboid];  % 后右下
v8 = v4 + [0, 0, width_cuboid];  % 后右上

% 可视化立方体和elements
figure;
hold on;
axis equal;
xlabel('X [mm]');
ylabel('Y [mm]');
zlabel('Z [mm]');
title('Cuboid with Outer Frame and Internal Volume');

% 绘制外围立方体框架 (12条边，用一种颜色)
cube_vertices = [x_range(1), y_range(1), z_range(1);  
                 x_range(1), y_range(1), z_range(2);
                 x_range(1), y_range(2), z_range(1);
                 x_range(1), y_range(2), z_range(2);
                 x_range(2), y_range(1), z_range(1);
                 x_range(2), y_range(1), z_range(2);
                 x_range(2), y_range(2), z_range(1);
                 x_range(2), y_range(2), z_range(2)];

% 外部框架 (灰色)
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

% 绘制内部小立方体 (用一种颜色填充)
fill3([v1(1), v3(1), v7(1), v5(1)], [v1(2), v3(2), v7(2), v5(2)], [v1(3), v3(3), v7(3), v5(3)], 'r', 'FaceAlpha', 0.3); % 前面
fill3([v2(1), v4(1), v8(1), v6(1)], [v2(2), v4(2), v8(2), v6(2)], [v2(3), v4(3), v8(3), v6(3)], 'r', 'FaceAlpha', 0.3); % 后面
fill3([v1(1), v2(1), v6(1), v5(1)], [v1(2), v2(2), v6(2), v5(2)], [v1(3), v2(3), v6(3), v5(3)], 'r', 'FaceAlpha', 0.3); % 左面
fill3([v3(1), v4(1), v8(1), v7(1)], [v3(2), v4(2), v8(2), v7(2)], [v3(3), v4(3), v8(3), v7(3)], 'r', 'FaceAlpha', 0.3); % 右面
fill3([v1(1), v2(1), v4(1), v3(1)], [v1(2), v2(2), v4(2), v3(2)], [v1(3), v2(3), v4(3), v3(3)], 'r', 'FaceAlpha', 0.3); % 顶面
fill3([v5(1), v6(1), v8(1), v7(1)], [v5(2), v6(2), v8(2), v7(2)], [v5(3), v6(3), v8(3), v7(3)], 'r', 'FaceAlpha', 0.3); % 底面

% 绘制 elements (黑色圆)
num_elements = size(element_xyz, 1); % 获取 element 数量
radius = 5; % 圆的半径 [mm]
theta = linspace(0, 2*pi, 100); % 圆的角度分布

for i = 1:num_elements
    % 获取每个 element 的中心位置
    center_x = element_xyz(i, 1);
    center_y = element_xyz(i, 2);
    center_z = element_xyz(i, 3);
    
    % 计算圆的 XZ 平面坐标
    circle_x = center_x + radius * cos(theta);
    circle_z = center_z + radius * sin(theta);
    
    % 绘制圆 (位于 XZ 平面，Y 坐标固定为 center_y)
    fill3(circle_x, center_y * ones(size(circle_x)), circle_z, 'k', 'EdgeColor', 'none'); % 黑色填充圆
end
hold off;

% ----------------------------------------------生成x, y,z轴的网格点,开始mesh处理----------------------------------------------
x_grid = adjust_range(x_range(1), x_range(2), grid_size); % mm
y_grid = adjust_range(y_range(1), y_range(2), grid_size);
z_grid = adjust_range(z_range(1), z_range(2), grid_size);
x_grid_first = x_grid(1);
y_grid_first = y_grid(1);
z_grid_first = z_grid(1);
% 输出x, y, z方向的grid数量
num_x_grid = length(x_grid);
num_y_grid = length(y_grid);
num_z_grid = length(z_grid);
fprintf('X方向的grid数量: %d\n', num_x_grid);
fprintf('Y方向的grid数量: %d\n', num_y_grid);
fprintf('Z方向的grid数量: %d\n', num_z_grid);
fprintf('整个计算环境的Total grid数量: %d\n', num_z_grid*num_y_grid*num_x_grid);

% 创建三维网格
[X, Y, Z] = meshgrid(x_grid, y_grid, z_grid);
% 将三维网格展开为坐标列表，并生成索引
grid_coordinates = [X(:), Y(:), Z(:)];% the coordinate value in the physcial world
xyz_first_grid=[x_grid(1), y_grid(1), z_grid(1)]; % 这个立方体，左上前这个角是第一个mesh (1 ，1，1), 三个轴上都是最小值
fprintf('整个环境的第一个mesh voxel的3D坐标(mm)是： %f, %f, %f\n\n', xyz_first_grid(1),xyz_first_grid(2),xyz_first_grid(3));

% 遍历每个 element，计算其落在哪个网格
element_grid_indices = zeros;
for idx = 1:num_elements
    % 当前 element 的坐标
    element_pos = element_xyz(idx, :); % [x, y, z]
    % 计算网格索引
    grid_index_x = round((element_pos(1) - xyz_first_grid(1)) / grid_size) + 1;%[mm] to index
    grid_index_y = round((element_pos(2) - xyz_first_grid(2)) / grid_size) + 1;
    grid_index_z = round((element_pos(3) - xyz_first_grid(3)) / grid_size) + 1;
    % 存储索引
    element_grid_indices(idx, 1:3) = [grid_index_x, grid_index_y, grid_index_z];
end
% 将结果保存为 CSV 文件
writematrix(element_grid_indices, 'element_grid_indices.csv');

% 初始化存储 Co2 cuboid 内部的 grid 索引，因为这个腔体是移动的，所以我们需要把它放入到一个for循环内部
for i=3:1:3 %竖直方向移动范围是4mm
    for j=6:1:6 %for j,  水平方向移动范围是10mm，从-5到0到+5，共11 (1:11).当前先不移动，只计算x未偏移的,处于中心位置
    fprintf('*********************************************************************************************************************************\n');
    fprintf('当前i=%d，j=%d,说明CO2腔体在x方向上是偏移中心轴%d mm，在y方向上偏移了初始态%d mm\n', i,j,j -6,i-1);

    % 求解第14个element的位置，然后利用这个信息，求解左上角顶点 [x, y, z]，从而得到整个co2立方体的位置
    element_xyz = readmatrix('element_xyz.csv'); % from "physical_coordinate.m"
    length_cuboid = 131;
    width_cuboid = 144.5;
    height_cuboid = 85.5;
    grid_size = 0.7208;

    
    z_14 = element_xyz(14,3) - 10; % Z dirction
    y_14 = element_xyz(14,2) - 16 - i; %当盒子竖直移动时候，初始态y值，初始值17，后面箱子会变，每次改变1mm
    x_14 = element_xyz(14,1) - 1 + j;  %当盒子水平移动时候，x值变化，当前先不考虑x变化，只在y变化，初始值5
    cuboid_bott_left = [x_14, y_14, z_14];
    
    %得到Co2立方体的8个定点
    v1 = cuboid_bott_left;           % 左下点,[mm]
    v2 = v1 - [0, height_cuboid, 0]; % 左上点
    v3 = v1 + [length_cuboid, 0, 0]; % 右下点
    v4 = v2 + [length_cuboid, 0, 0]; % 右上点
    v5 = v1 + [0, 0, width_cuboid];  % 后左下
    v6 = v2 + [0, 0, width_cuboid];  % 后左上
    v7 = v3 + [0, 0, width_cuboid];  % 后右下
    v8 = v4 + [0, 0, width_cuboid];  % 后右上

    % 创建 CO2 cuboid 的顶点矩阵，确定 cuboid 的边界范围
    x_min_cuboid = v1(1,1); % x 方向的最小值
    x_max_cuboid = v3(1,1); % x 方向的最大值
    y_min_cuboid = v2(1,2); % y 方向的最小值
    y_max_cuboid = v1(1,2); % y 方向的最大值
    z_min_cuboid = v1(1,3); % z 方向的最小值
    z_max_cuboid = v5(1,3); % z 方向的最大值
    xyz_first_grid=[-52.199, -59.883583, 156.1];
    fprintf('  Cuboid CO2 第一个点grid的3D坐标是: %f, %f, %f\n', x_min_cuboid,y_min_cuboid,z_min_cuboid);
    fprintf('  Cuboid CO2 第一个点grid的对应索引: %f, %f, %f\n', round((x_min_cuboid - xyz_first_grid(1)) / grid_size) + 1,round((y_min_cuboid - xyz_first_grid(2)) / grid_size) + 1,round((z_min_cuboid - xyz_first_grid(3)) / grid_size) + 1);

    %  计算xyz三个方向，Co2腔体的最小和最大的那个点的grid index
    grid_min_x = round((x_min_cuboid - xyz_first_grid(1)) / grid_size) + 1; % [mm] to index
    grid_max_x = round((x_max_cuboid - xyz_first_grid(1)) / grid_size) + 1;

    grid_min_y = round((y_min_cuboid - xyz_first_grid(2)) / grid_size) + 1; % [mm] to index
    grid_max_y = round((y_max_cuboid - xyz_first_grid(2)) / grid_size) + 1;

    grid_min_z = round((z_min_cuboid - xyz_first_grid(3)) / grid_size) + 1; % [mm] to index
    grid_max_z = round((z_max_cuboid - xyz_first_grid(3)) / grid_size) + 1;

    number_voexl = (grid_max_x-grid_min_x+1)*(grid_max_y-grid_min_y+1)*(grid_max_z-grid_min_z+1);
    fprintf('  Cuboid CO2 三个方向的grid点数: %f, %f, %f\n', grid_max_x-grid_min_x+1,grid_max_y-grid_min_y+1,grid_max_z-grid_min_z+1);
    % 初始化存储 cuboid Co2内部的 grid 索引
    inside_cuboid_indices_3d = zeros(number_voexl,3);
    % 遍历 cuboid 范围，生成索引
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

    % 输出 cuboid 内部的 grid 数量
    fprintf('  Cuboid CO2 内部的 grid 数量: %d\n', size(inside_cuboid_indices_3d, 1));
    % 将结果保存为 CSV 文件
    writematrix(inside_cuboid_indices_3d, 'inside_cuboid_indices_3d.csv');

    for ii=33:1:63  %% vertical direction: 34mm to 94 mm

        xyz_first_grid=[-52.199, -59.883583, 156.1];
        num_x_grid = 217;
        num_y_grid = 175;
        num_z_grid = 231;
        x_grid_first = xyz_first_grid(1);
        y_grid_first = xyz_first_grid(2);
        z_grid_first = xyz_first_grid(3);

        %计算array的中心位置正上方44mm处位置为起始出发点
        element_xyz = readmatrix('element_xyz.csv');
        grid_size = 0.7208;% [mm]
        fprintf('=============================================================================================\n');
        fprintf('当第KK=%d轮计算,\n', ii);
        if (mod(ii,7)==0) % 7
            begin_x = 0.5 * element_xyz(85,1) + 0.5 * element_xyz(99,1) -3*1.5 + (7-1)*1.5; % 85是第7列第一个，99是第8列第1个
            begin_y = element_xyz(85,2) - 34-(floor((ii)/7)-1)*1.5; % 正上方34mm处 start
        elseif(mod(ii,7)) % 1,2,3,4,5,6
            begin_x = 0.5 * element_xyz(85,1) + 0.5 * element_xyz(99,1) -3*1.5 + (mod(ii,7)-1)*1.5; % 85是第7列第一个，99是第8列第1个
            begin_y = element_xyz(85,2) - 34-floor((ii)/7)*1.5; % 正上方34mm处 start
        end
        begin_z = 0.5 * element_xyz(91,3) + 0.5 *element_xyz(92,3);% center
        fprintf('  array中心点空间位置 (%f mm,%f mm,%f mm)\n', 0.5 * element_xyz(85,1) + 0.5 * element_xyz(99,1),element_xyz(85,2),0.5 * element_xyz(91,3) + 0.5 *element_xyz(92,3))
        fprintf('  出发目标点空间位置 (%f mm,%f mm,%f mm)\n', begin_x,begin_y,begin_z)
        fprintf('  出发目标点空间位置y值是远离array(y=%f mm)是%f mm, x远离center(x=%f mm)是 %f mm\n', element_xyz(85,2),element_xyz(85,2)-begin_y,0.5 * element_xyz(85,1) + 0.5 * element_xyz(99,1),begin_x-0.5 * element_xyz(85,1) - 0.5 * element_xyz(99,1))
        %把这些值转化为grid
        begin_x_grid =  round((begin_x - xyz_first_grid(1)) / grid_size) + 1;
        begin_y_grid =  round((begin_y - xyz_first_grid(2)) / grid_size) + 1;
        begin_z_grid =  round((begin_z - xyz_first_grid(3)) / grid_size) + 1;
        fprintf('  出发点对应的grid索引是: %d,%d,%d\n', begin_x_grid,begin_y_grid,begin_z_grid);
        save('grid_data.mat', 'begin_x_grid', 'begin_y_grid', 'begin_z_grid','i','j','ii','num_x_grid','num_y_grid','num_z_grid','x_grid_first', 'y_grid_first', 'z_grid_first');
        %调用k—wave函数计算
        run('Copy_of_k_wave_multi_media_cuboid.m');
        run('extract_maxpoint.m');
    end
    end % for j
end

diary off;

% 调整网格范围的函数，处理边缘体素
function adjusted_range = adjust_range(range_start, range_end, grid_size)
    total_range = range_end - range_start;
    remainder = mod(total_range, grid_size);
    if remainder >= (grid_size / 2)
        range_end = range_end + (grid_size - remainder);  % 拓展以适应完整的体素
    else
        range_end = range_end - remainder;  % 剔除边缘较小的部分
    end
    adjusted_range = (range_start + grid_size/2) : grid_size: range_end;
end
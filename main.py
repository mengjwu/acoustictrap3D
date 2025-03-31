# -*- coding: utf-8 -*-
import numpy as np
import serial
import time
import numpy as np
import pandas as pd
from datetime import datetime
import pickle
from sklearn.preprocessing import StandardScaler


# Load the trained model from the file
with open('my_model.pkl', 'rb') as f:
    loaded_model = pickle.load(f)
with open('scalerx.pkl', 'rb') as f:
    scalerx = pickle.load(f)
with open('scalery.pkl', 'rb') as f:
    scalery = pickle.load(f)


# 初始化串口
ser = serial.Serial(
    port='/dev/ttyUSB1',
    baudrate=460800,
    parity=serial.PARITY_NONE,
    stopbits=serial.STOPBITS_ONE,
    bytesize=serial.EIGHTBITS,
    timeout=1
)

df = pd.read_csv('mark_start.csv', header=None)
mark2 = df.iloc[:, 0:3].values
mark3 = df.iloc[:, 3:6].values
# 对每一列进行梳理，保留小数点后2位
media_start2_x = [round(x, 2) for x in mark2[:, 0]]
media_start2_y = [round(y, 2) for y in mark2[:, 1]]
media_start2_z = [round(z, 2) for z in mark2[:, 2]]
media_start3_x = [round(x, 2) for x in mark3[:, 0]]
media_start3_y = [round(y, 2) for y in mark3[:, 1]]
media_start3_z = [round(z, 2) for z in mark3[:, 2]]


df1 = pd.read_csv('expected.csv', header=None)
expected_target=np.zeros((1, 3))
expected_target[0,0] = df1.iloc[0, 0]
expected_target[0,1] = df1.iloc[0, 1]
expected_target[0,2] = df1.iloc[0, 2]

# 提取每列的最频繁值作为代表值
def most_frequent(lst):
    return max(set(lst), key=lst.count)
# 获取每列最频繁的值
media_start2_x_val = most_frequent(media_start2_x)
media_start2_y_val = most_frequent(media_start2_y)
media_start2_z_val = most_frequent(media_start2_z)
media_start3_x_val = most_frequent(media_start3_x)
media_start3_y_val = most_frequent(media_start3_y)
media_start3_z_val = most_frequent(media_start3_z)

# 比较上面求得的z值是否一致，误差小于1,大于1则终止整个程序
z_difference = abs(media_start2_z_val - media_start3_z_val)
if z_difference > 1:
    print(f"Error: 两个标记点Z值差别大，CO2位置倾斜需要调整，终止程序.")
    exit()  # 终止程序
else:
    print(f"CO2 box is vertical ， difference of two marks is [mm]: {z_difference}.\n")


# 首先发送 phase，让 ball 到达 target 点
# 初始化延迟矩阵
delay = np.zeros((16, 203))  # 15行203列的零矩阵
for i in range(1, 17):
    filename = f"/home/mjwu/PixelXYZ_20250221/trap_py_learning/target4/{i}.csv"
    num = np.loadtxt(filename, delimiter=',', dtype=int).T  # 读取并转置
    delay[i - 1, :] = num[1, :]  # 第二列的数据赋值到延迟矩阵

# 循环发送每一行
for ii in range(16):
    # 创建发送数据：前后添加标志位
    # 确保串口已打开
    if not ser.is_open:
        ser.open()
    print(f"Sending data for line {ii + 1}...")

    # 创建发送数据：前后添加标志位
    data_to_send = bytearray([170] + [int(delay[ii, j]) for j in range(203)])
    # 按照数据包发送数据，比但个字节速度快很多
    ser.write(data_to_send)
    ser.flush()  # 确保204 bytes已发送完毕

    # hex_values = [f"{byte:02X}" for byte in data_to_send]
    # print("Hex data:", " ".join(hex_values))

    # # 按字节发送数据
    # for j in range(0, 204):
    #     if j == 0:  # 如果是第一个字节，发送标志位170
    #         byte_int = 170
    #     elif j > 0 and j <= 203:  # 否则发送 delay 的值
    #         byte_int = int(delay[ii, j - 1]) # 获取 delay(ii, j-1) 的值并转换为整数
    #     # 发送数据
    #     ser.write(bytes([byte_int]))  # 直接发送单字节
    #     time.sleep(4e-6) # 230400 每个byte发送大约48us，然后延迟 4 us
    #     # hex_values.append(f"{byte_int:02X}")  # 将字节值转换为两位的十六进制字符串
    #     # 打印 204 个字节的十六进制值


    print(f"Done Line {ii + 1}...")
    if ii == 0:
        time.sleep(3) #预留3s时间用于放置小球

    time.sleep(1.5) # 每1轮phase update后，都要延迟1s

# 关闭串口
ser.close()
xy = np.zeros((196, 3))  # 假设需要14x14的网格
k = 39.37
pitch = 10.5  # mm

# 手动设置几个特定的值
xy[0, 0] = (-3198.8 - 3002.0) / (2 * k)
xy[0, 1] = 3513.8 / k
xy[1, 0] = (-3198.8 - 3002.0) / (2 * k)
xy[1, 1] = 3100.4 / k
xy[13, 0] = (-3198.8 - 3002.0) / (2 * k)
xy[13, 1] = -1860.2 / k
xy[14, 0] = (-2785.4 - 2588.6) / (2 * k)
xy[14, 1] = 3513.8 / k

# 构建其余的坐标矩阵
index = 0  # 索引变量，从0开始
for i in range(14):  # 14行
    for j in range(14):  # 14列
        xy[index, 0] = -6.5*pitch + i*pitch
        xy[index, 1] = 6.5*pitch - j*pitch
        index += 1

# 设置第三列为0（保持为二维平面）
xy[:, 2] = 0

# 保存为 CSV 文件
np.savetxt("xyz_coordinates.csv", xy, delimiter=',', comments='')

# 读取 calibration.csv 和 sort_index_edited.csv 文件
f3 = pd.read_csv('calibration.csv')  # 读取 calibration.csv
reson_calibration = np.zeros((f3.shape[0], f3.shape[1]))  # 初始化数组，假设与 f3 形状一致
f2 = pd.read_csv('sort_index_edited.csv', header=None)  # 读取 sort_index_edited.csv
F3 = f3.iloc[:, 0:4].to_numpy()  # 提取前四列并转换为 NumPy 数组，单位为微秒（us）
# 初始化变量
reson_calibration = np.zeros(196, dtype=int)
# 计算 reson_calibration，单位转换为 ns，步长为 100 ns
for i in range(196):
    reson_calibration[i] = round(F3[i, 3] * 1000 / 100)  # 第四列是延迟时间，单位 us 转换为 ns

# 定义常量
T = 25000  # ns
quarter_cycle_ns = 6250  # 四分之一个周期，单位 ns

# 初始化上一次的 err_x 和 err_y
previous_err_x = None
previous_err_y = None
mark_x_drif = None
mark_y_drif = None
mark_z_drif = None

time.sleep(26)

df = pd.read_csv('final.csv', header=None)
target = np.zeros((1, 3))  # 1x3 数组
target_new = np.zeros((1, 3))  # 1x3 数组
target[0, 0] = df.iloc[0, 0]  # X
target[0, 1] = df.iloc[0, 1]  # Y
target[0, 2] = 244  # Z, 这个Z值需要看下dataset，修改下，他们之间有1mm以下误差
print(f"\nEnter Compensation Mode")  # 打印当前时间

file_index = 1

while True:
    # 获取当前时间并格式化为 HH:mm:ss.SSS
    # current_time = datetime.now().strftime('%H:%M:%S.%f')[:-3]
    # print(f"Current time in loop: {current_time}")  # 打印当前时间

    # file_name = f'drift_data{file_index}.csv'
    file_name = f'drift_data.csv'
    # print(f"\nStart New Round of Update.")

    try:
            data = pd.read_csv(file_name, header=None)  # 不使用表头
            # file_index = 2 if file_index == 1 else 1
            # print(data)

            if data.empty:
               raise ValueError("File is empty")

            if data.shape[1] != 3:
                raise ValueError("File does not have exactly 3 columns")

            data1 = pd.read_csv("mark.csv", header=None)  # 没表头

            if (data1.isin([float('inf')])).any().any():  # 如果有inf值
                raise ValueError("Data Contains of 'Inf'，Wait 2ms.")

            # 检查列数是否小于3
            if data1.shape[1] < 3:
                raise ValueError("Data col is less than 3，Wait2ms.")

            startTime0 = time.time()

            if data.iloc[0, 2] != 1:
                print("Data is Invalid, Flag bit is 0\n")
                break

            err_x = data.iloc[0, 0]  # 第一列
            err_y = data.iloc[0, 1]  # 第二列
            flag  = data.iloc[0, 2]  # 第三列

            if previous_err_x == err_x and previous_err_y == err_y:
                # print("Qt Has Not Updated Error Values, Skip This Round.")
                time.sleep(0.001)  # 等待 1 毫秒
                continue  # 跳过本次循环


            # 更新上一轮的 err_x 和 err_y
            previous_err_x = err_x
            previous_err_y = err_y

            # 打印 x 和 y 的值
            print(f"err_x:  {err_x}, err_y:  {err_y}")

            # 判断 err_y 是否小于 1 的绝对值
            if abs(err_y) < 1.0 and abs(err_x) < 0.7:
                print("Error < Threshold, Skip This Round.")
                time.sleep(0.001)  # 等待 5 毫秒
                continue  # 跳过本次循环
            # else:
            #     print('误差大于最低阈值，往下执行.')

            if flag == 0:  # 假设我们只关心第一行的 flag
                print('Flag bit is 0，Skip This Round.')
                time.sleep(0.01)  # 短暂等待，避免高频无意义读取
                break
            # else:
            #     print('误差值合法，往下执行.')

            mark2_x = data1.iloc[0, 0]  # 读取动态的mark坐标， 第一列X
            mark2_y = data1.iloc[0, 1]  # 第二列Y
            mark2_z = data1.iloc[0, 2]  # 第三列Z
            mark3_x = data1.iloc[0, 3]  # 第一列
            mark3_y = data1.iloc[0, 4]  # 第二列
            mark3_z = data1.iloc[0, 5]  # 第三列

            mark_x_drif = (mark2_x + mark3_x)/2 - (media_start2_x_val + media_start3_x_val)/2
            mark_y_drif = (mark2_y + mark3_y)/2 - (media_start2_y_val + media_start3_y_val)/2
            mark_z_drif = (mark2_z + mark3_z)/2 - (media_start2_z_val + media_start3_z_val)/2
            # print(f"CO2 initial Position is : x = {(media_start2_x_val + media_start3_x_val)/2} mm, y = {(media_start2_y_val + media_start3_y_val)/2} mm, z = {(media_start2_z_val + media_start3_z_val)/2} mm.")
            # print(f"CO2 current Position is : x = {(mark2_x + mark3_x)/2} mm, y = {(mark2_y + mark3_y)/2} mm, z = {(mark2_z + mark3_z)/2} mm.")

            media_x = -37.6  # CO'2箱子会动，所以x也会动起来，但是x移动不影响模型，而且我们采样当前数据集不包含x的移动，所以暂且是fixed
            media_y = round(-46.8 + mark_y_drif,1)
            media_z = 165.7  # 箱子只在x和y方向移动，并不会存在z的移动
            print(f"CO2 空间位移 is : x = {mark_x_drif} mm, y = {mark_y_drif} mm, z = {mark_z_drif} mm.")

            num_steps = max(round(abs(err_y)/1), round(abs(err_x)/0.7))
            k = 1

            if num_steps > 10:
                print("Error too big, Skip this round.\n")
                continue  # 跳过本次循环

            else:
                for k in range(1, 2):  #(1，n) 表示执行 n-1 次
                    # print("误差属于合理范围，往下执行，开始计算Phase.")
                    x_need = 0
                    y_need = 0

                    if err_x >= 0 and 0.7 * k < err_x:  # 当误差是大于0.7，则调整1次后后不再调整
                        target_new[0, 0] = round(expected_target[0, 0] - 0.7 * k, 1)  # X 坐标
                    elif err_x >= 0 and 0.7 * k >= err_x:  # 当误差是小于0.7，不调整
                        x_need = 1  # 不需要在x方向上调整
                        target_new[0, 0] = expected_target[0, 0]
                    elif err_x < 0 and 0.7 * k < abs(err_x):
                        target_new[0, 0] = round(expected_target[0, 0] + 0.7 * k, 1)
                    elif err_x < 0 and 0.7 * k >= abs(err_x):  # 当误差是小于0.7，不调整
                        x_need = 1  # 不需要在x方向上调整
                        target_new[0, 0] = expected_target[0, 0]

                    if err_y >= 0 and err_y <= 1.5 and 1 * k <= err_y:  # 当误差是大于等于1，小于等于1.5，则每次调整1 mm
                        target_new[0, 1] = round(expected_target[0, 1] - 1 * k, 1)
                    elif err_y >= 0 and err_y > 1.5 and err_y < 2:  # 当误差是大于1.5小于2，则每次调整err_y，不再续调整第2次
                        target_new[0, 1] = round(expected_target[0, 1] - err_y, 1)
                    elif err_y >= 2:  # 当误差是大于等于2，则调整2 mm，不再续调整第2次
                        target_new[0, 1] = round(expected_target[0, 1] - 2, 1)
                    elif err_y >= 0 and 1 * k > err_y:  # 当误差是小于1.0，不调整
                        y_need = 1  # 不需要在y方向上调整
                        target_new[0, 1] = round(expected_target[0, 1], 1)

                    elif err_y < 0 and 1 * k <= abs(err_y) and abs(err_y) < 1.5:  # 当误差是大于等于1，小于等于1.5，则每次调整1 mm
                        target_new[0, 1] = round(expected_target[0, 1] + 1 * k, 1)
                    elif err_y < 0 and abs(err_y) >= 1.5 and abs(err_y) < 2:  # 当误差是大于1.5小于2，则每次调整err_y, 不再续调整第2次
                        target_new[0, 1] = round(expected_target[0, 1] + abs(err_y), 1)
                    elif err_y < 0 and abs(err_y) >= 2:  # 当误差是大于等于2，则每次调整2 mm, 不再续调整第2次
                        target_new[0, 1] = round(expected_target[0, 1] + 2, 1)
                    elif err_y < 0 and 1 * k > abs(err_y):  # 当误差是小于1.0，不做调整
                        y_need = 1  # 不需要在y方向上调整
                        target_new[0, 1] = round(expected_target[0, 1], 1)



                    if x_need == 1 and y_need == 1:
                        break  # 终止当前的 for 循环

                    target_new[0, 2] = expected_target[0, 2]  # z值不变化，和训练样本的一致

                    # 计算每个点到目标点的距离和 TOF
                    tof = np.zeros(196)
                    dist = np.zeros(196)

                    print(f"target理论期待的目标位置 is : x = {expected_target[0, 0]} mm, y = {expected_target[0, 1]} mm, z = {expected_target[0, 2]} mm.")
                    print(f"target修正后输入ML预测位置 is : x = {target_new[0, 0]} mm, y = {target_new[0, 1]} mm, z = {target_new[0, 2]} mm.")

                    startTime = time.time()
                    data_testing_x = np.array([[media_x, media_y, media_z, target_new[0, 0], target_new[0, 1], target_new[0, 2]]])  # 这里预测的结果是us，不是ns
                    data_testing = scalerx.transform(data_testing_x)
                    data_testing_y = loaded_model.predict(data_testing)
                    tof = scalery.inverse_transform(data_testing_y) # prediction is us
                    tof = tof.flatten() #将ToF的二维数组，转化为1维
                    tof = 1000 * tof  # 单位转换为ns这样，这样可以确保下面的代码不要每一处都修改
                    stopTime = time.time()  # 结束计时
                    print("ToF prediction time is [units: ms]", 1000*(stopTime - startTime))  # 打印用时
                    # print("Prediction Results are [units: nanosecond]",tof)

                    # Time reversed处理，250个step，每个步长100ns，也就是0.1us
                    tof_r = np.zeros(196)
                    tof_r_code = np.zeros(196)

                    # tmax = np.max(tof) # for air computation
                    tmax = np.max(tof)
                    index_ele = np.argmax(tof)
                    tmax = tmax + 100 # 加上100ns
                    for i in range(196):
                        tof_r[i] = tmax - tof[i]  # air
                        tof_r_code[i] = np.round(np.mod(tmax - tof[i], 250 * 100) / 100)  # 40KHz means 周期是25 us = 250000 ns

                    # 另外twintrap还需要在左右2侧加一个signature，左侧加pi / 2, 右侧减去pi / 2
                    # 四分之一个周期，就是25000 / 4 = 6250 ns
                    # 初始化变量
                    trap_code = np.zeros(196, dtype=int)
                    sign_code = np.zeros(196, dtype=int)
                    sign = np.zeros(196, dtype=int)

                    # Twin trap pattern
                    for i in range(98):
                        sign_code[i] = round(quarter_cycle_ns / 100)  # 左侧 +π/2
                        sign_code[i + 98] = -round(quarter_cycle_ns / 100)  # 右侧 -π/2
                        sign[i] = quarter_cycle_ns
                        sign[i + 98] = -quarter_cycle_ns

                    # 加入 holographic signature
                    for i in range(196):
                        trap_code[i] = round((tof_r[i] + sign[i]) % T / 100)  # 范围 1~250， tof_r单位ns，sign也是ns，T 是25000 ns
                    trap_code_calibration = np.zeros(196, dtype=int)

                    # 校正 trap_code
                    for i in range(196):
                        trap_code_calibration[i] = trap_code[i] - reson_calibration[i]
                        if trap_code_calibration[i] < 0:  # 如果校正后的值小于 0，则加上周期 250
                            trap_code_calibration[i] += 250

                    index = f2.iloc[:, 0:2].to_numpy()  # 提取前两列作为索引
                    if not ser.is_open:
                        ser.open()
                    # 发送字符 170 (0xAA)+ 203 bytes
                    data_to_send = bytearray([170] + [int(trap_code_calibration[index[i][0] - 1]) for i in range(203)])  # 按照数据包发送
                    ser.write(data_to_send)
                    ser.flush()
                    # ser.write(bytes([170]))
                    # for i in range(len(index)):
                    #     value = trap_code_calibration[index[i][0] - 1]  # 获取对应值（MATLAB 索引从 1 开始）
                    #     ser.write(bytes([value]))  # 写入单字节数据
                    #     time.sleep(5e-6) #每个byte后，都停5us

                    # 关闭串口
                    ser.close()

                    # 打印当前步数和更新后的目标坐标
                    print(f"Step {k}: The expected target (mm): {target[0, 0]}, {target[0, 1]}, {target[0, 2]}. The updated target: {target_new[0, 0]}, {target_new[0, 1]}, {target_new[0, 2]}. UART Done.")
                    # time.sleep(5e-3)  # 发送完一次phase update, 则休息以下，确保小球稳定再移动
                    stopTime0 = time.time()  # 结束计时
                    print("Python frame cost time is [units: ms]", 1000 * (stopTime0 - startTime0))  # 打印用时

    except (pd.errors.EmptyDataError, ValueError) as e:
            print(f"\nWarning: {e}. Retrying in 2ms...\n")
            time.sleep(0.002)  # 暂停2ms
            continue  # 跳过此次循环，继续下一次

    # 切换到另一个文件
    file_index = 2 if file_index == 1 else 1


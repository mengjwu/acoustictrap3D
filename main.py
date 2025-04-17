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


def most_frequent(lst):
    return max(set(lst), key=lst.count)

media_start2_x_val = most_frequent(media_start2_x)
media_start2_y_val = most_frequent(media_start2_y)
media_start2_z_val = most_frequent(media_start2_z)
media_start3_x_val = most_frequent(media_start3_x)
media_start3_y_val = most_frequent(media_start3_y)
media_start3_z_val = most_frequent(media_start3_z)


z_difference = abs(media_start2_z_val - media_start3_z_val)
if z_difference > 1:
     
    exit()  
else:
    print(f"CO2 box is vertical ， difference of two marks is [mm]: {z_difference}.\n")


 
delay = np.zeros((16, 203))   
for i in range(1, 17):
    filename = f"/home/mjwu/PixelXYZ_20250221/trap_py_learning/target4/{i}.csv"
    num = np.loadtxt(filename, delimiter=',', dtype=int).T 
    delay[i - 1, :] = num[1, :]  

 
for ii in range(16):
 
    if not ser.is_open:
        ser.open()
    print(f"Sending data for line {ii + 1}...")
    data_to_send = bytearray([170] + [int(delay[ii, j]) for j in range(203)])
    ser.write(data_to_send)
    ser.flush() 

    # hex_values = [f"{byte:02X}" for byte in data_to_send]
    # print("Hex data:", " ".join(hex_values))


    # for j in range(0, 204):
    #     if j == 0: 
    #         byte_int = 170
    #     elif j > 0 and j <= 203: 
    #         byte_int = int(delay[ii, j - 1]) 
    #     ser.write(bytes([byte_int]))  
    #     time.sleep(4e-6)
    #     # hex_values.append(f"{byte_int:02X}")  


    print(f"Done Line {ii + 1}...")
    if ii == 0:
        time.sleep(3) 

    time.sleep(1.5) 

 
ser.close()
xy = np.zeros((196, 3))  
k = 39.37
pitch = 10.5   

 
xy[0, 0] = (-3198.8 - 3002.0) / (2 * k)
xy[0, 1] = 3513.8 / k
xy[1, 0] = (-3198.8 - 3002.0) / (2 * k)
xy[1, 1] = 3100.4 / k
xy[13, 0] = (-3198.8 - 3002.0) / (2 * k)
xy[13, 1] = -1860.2 / k
xy[14, 0] = (-2785.4 - 2588.6) / (2 * k)
xy[14, 1] = 3513.8 / k

 
index = 0   
for i in range(14):  
    for j in range(14):   
        xy[index, 0] = -6.5*pitch + i*pitch
        xy[index, 1] = 6.5*pitch - j*pitch
        index += 1

xy[:, 2] = 0

 
np.savetxt("xyz_coordinates.csv", xy, delimiter=',', comments='')

f3 = pd.read_csv('calibration.csv')  
reson_calibration = np.zeros((f3.shape[0], f3.shape[1]))   
f2 = pd.read_csv('sort_index_edited.csv', header=None)  
F3 = f3.iloc[:, 0:4].to_numpy()   

reson_calibration = np.zeros(196, dtype=int)
for i in range(196):
    reson_calibration[i] = round(F3[i, 3] * 1000 / 100)  

 
T = 25000  # ns
quarter_cycle_ns = 6250  #  ns

 
previous_err_x = None
previous_err_y = None
mark_x_drif = None
mark_y_drif = None
mark_z_drif = None

time.sleep(26)

df = pd.read_csv('final.csv', header=None)
target = np.zeros((1, 3))  # 1x3  
target_new = np.zeros((1, 3))  # 1x3  
target[0, 0] = df.iloc[0, 0]  # X
target[0, 1] = df.iloc[0, 1]  # Y
target[0, 2] = 244  # Z 
print(f"\nEnter Compensation Mode")   

file_index = 1

while True:
  
    # current_time = datetime.now().strftime('%H:%M:%S.%f')[:-3]
    # print(f"Current time in loop: {current_time}")  

    # file_name = f'drift_data{file_index}.csv'
    file_name = f'drift_data.csv'
    # print(f"\nStart New Round of Update.")

    try:
            data = pd.read_csv(file_name, header=None)  
            # file_index = 2 if file_index == 1 else 1
            # print(data)

            if data.empty:
               raise ValueError("File is empty")

            if data.shape[1] != 3:
                raise ValueError("File does not have exactly 3 columns")

            data1 = pd.read_csv("mark.csv", header=None)   

            if (data1.isin([float('inf')])).any().any():  
                raise ValueError("Data Contains of 'Inf'，Wait 2ms.")
 
            if data1.shape[1] < 3:
                raise ValueError("Data col is less than 3，Wait2ms.")

            startTime0 = time.time()

            if data.iloc[0, 2] != 1:
                print("Data is Invalid, Flag bit is 0\n")
                break

            err_x = data.iloc[0, 0]   
            err_y = data.iloc[0, 1]   
            flag  = data.iloc[0, 2]  

            if previous_err_x == err_x and previous_err_y == err_y:
                # print("Qt Has Not Updated Error Values, Skip This Round.")
                time.sleep(0.001)   
                continue   

            previous_err_x = err_x
            previous_err_y = err_y

 
            print(f"err_x:  {err_x}, err_y:  {err_y}")

 
            if abs(err_y) < 1.0 and abs(err_x) < 0.7:
                print("Error < Threshold, Skip This Round.")
                time.sleep(0.001)  # 5 ms
                continue   
 
 

            if flag == 0:   
                print('Flag bit is 0，Skip This Round.')
                time.sleep(0.01)   
                break
            # else:
            #     print('wait.')

            mark2_x = data1.iloc[0, 0]  # 
            mark2_y = data1.iloc[0, 1]  #  Y
            mark2_z = data1.iloc[0, 2]  #  Z
            mark3_x = data1.iloc[0, 3]  #  
            mark3_y = data1.iloc[0, 4]  #  
            mark3_z = data1.iloc[0, 5]  #  

            mark_x_drif = (mark2_x + mark3_x)/2 - (media_start2_x_val + media_start3_x_val)/2
            mark_y_drif = (mark2_y + mark3_y)/2 - (media_start2_y_val + media_start3_y_val)/2
            mark_z_drif = (mark2_z + mark3_z)/2 - (media_start2_z_val + media_start3_z_val)/2
            # print(f"CO2 initial Position is : x = {(media_start2_x_val + media_start3_x_val)/2} mm, y = {(media_start2_y_val + media_start3_y_val)/2} mm, z = {(media_start2_z_val + media_start3_z_val)/2} mm.")
            # print(f"CO2 current Position is : x = {(mark2_x + mark3_x)/2} mm, y = {(mark2_y + mark3_y)/2} mm, z = {(mark2_z + mark3_z)/2} mm.")

            media_x = -37.6  
            media_y = round(-46.8 + mark_y_drif,1)
            media_z = 165.7  
            

            num_steps = max(round(abs(err_y)/1), round(abs(err_x)/0.7))
            k = 1

            if num_steps > 10:
                print("Error too big, Skip this round.\n")
                continue  

            else:
                for k in range(1, 2):  
                   
                    x_need = 0
                    y_need = 0

                    if err_x >= 0 and 0.7 * k < err_x:  
                        target_new[0, 0] = round(expected_target[0, 0] - 0.7 * k, 1)  # X  
                    elif err_x >= 0 and 0.7 * k >= err_x:   
                        x_need = 1  
                        target_new[0, 0] = expected_target[0, 0]
                    elif err_x < 0 and 0.7 * k < abs(err_x):
                        target_new[0, 0] = round(expected_target[0, 0] + 0.7 * k, 1)
                    elif err_x < 0 and 0.7 * k >= abs(err_x):  
                        x_need = 1  
                        target_new[0, 0] = expected_target[0, 0]

                    if err_y >= 0 and err_y <= 1.5 and 1 * k <= err_y:   
                        target_new[0, 1] = round(expected_target[0, 1] - 1 * k, 1)
                    elif err_y >= 0 and err_y > 1.5 and err_y < 2:  
                        target_new[0, 1] = round(expected_target[0, 1] - err_y, 1)
                    elif err_y >= 2:  
                        target_new[0, 1] = round(expected_target[0, 1] - 2, 1)
                    elif err_y >= 0 and 1 * k > err_y:  
                        y_need = 1  
                        target_new[0, 1] = round(expected_target[0, 1], 1)

                    elif err_y < 0 and 1 * k <= abs(err_y) and abs(err_y) < 1.5:  
                        target_new[0, 1] = round(expected_target[0, 1] + 1 * k, 1)
                    elif err_y < 0 and abs(err_y) >= 1.5 and abs(err_y) < 2:  
                        target_new[0, 1] = round(expected_target[0, 1] + abs(err_y), 1)
                    elif err_y < 0 and abs(err_y) >= 2:  
                        target_new[0, 1] = round(expected_target[0, 1] + 2, 1)
                    elif err_y < 0 and 1 * k > abs(err_y):   
                        y_need = 1  
                        target_new[0, 1] = round(expected_target[0, 1], 1)



                    if x_need == 1 and y_need == 1:
                        break   

                    target_new[0, 2] = expected_target[0, 2]  

                    tof = np.zeros(196)
                    dist = np.zeros(196)

                    startTime = time.time()
                    data_testing_x = np.array([[media_x, media_y, media_z, target_new[0, 0], target_new[0, 1], target_new[0, 2]]])  
                    data_testing = scalerx.transform(data_testing_x)
                    data_testing_y = loaded_model.predict(data_testing)
                    tof = scalery.inverse_transform(data_testing_y) # prediction is us
                    tof = tof.flatten() 
                    tof = 1000 * tof  
                    stopTime = time.time()  
                    print("ToF prediction time is [units: ms]", 1000*(stopTime - startTime))   
                    # print("Prediction Results are [units: nanosecond]",tof)

                  
                    tof_r = np.zeros(196)
                    tof_r_code = np.zeros(196)

                    # tmax = np.max(tof) # for air computation
                    tmax = np.max(tof)
                    index_ele = np.argmax(tof)
                    tmax = tmax + 100 #  100ns
                    for i in range(196):
                        tof_r[i] = tmax - tof[i]  # air
                        tof_r_code[i] = np.round(np.mod(tmax - tof[i], 250 * 100) / 100)  # 40KHz means  T= 25 us = 250000 ns

              
                    trap_code = np.zeros(196, dtype=int)
                    sign_code = np.zeros(196, dtype=int)
                    sign = np.zeros(196, dtype=int)

                    # Twin trap pattern
                    for i in range(98):
                        sign_code[i] = round(quarter_cycle_ns / 100)  #   +π/2
                        sign_code[i + 98] = -round(quarter_cycle_ns / 100)  #   -π/2
                        sign[i] = quarter_cycle_ns
                        sign[i + 98] = -quarter_cycle_ns

                    #   holographic signature
                    for i in range(196):
                        trap_code[i] = round((tof_r[i] + sign[i]) % T / 100)  #   1~250， tof_r ns，sign ns，T  25000 ns
                    trap_code_calibration = np.zeros(196, dtype=int)

                    #   trap_code
                    for i in range(196):
                        trap_code_calibration[i] = trap_code[i] - reson_calibration[i]
                        if trap_code_calibration[i] < 0:   
                            trap_code_calibration[i] += 250

                    index = f2.iloc[:, 0:2].to_numpy()  
                    if not ser.is_open:
                        ser.open()
                        
                    
                    data_to_send = bytearray([170] + [int(trap_code_calibration[index[i][0] - 1]) for i in range(203)])   
                    ser.write(data_to_send)
                    ser.flush()
 
                    ser.close()

                   
                    print(f"Step {k}: The expected target (mm): {target[0, 0]}, {target[0, 1]}, {target[0, 2]}. The updated target: {target_new[0, 0]}, {target_new[0, 1]}, {target_new[0, 2]}. UART Done.")
                    # time.sleep(5e-3)  
                    stopTime0 = time.time()  
                    print("Python frame cost time is [units: ms]", 1000 * (stopTime0 - startTime0))   

    except (pd.errors.EmptyDataError, ValueError) as e:
            print(f"\nWarning: {e}. Retrying in 2ms...\n")
            time.sleep(0.002)   
            continue  
 
    file_index = 2 if file_index == 1 else 1


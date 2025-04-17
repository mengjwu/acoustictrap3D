**This repository includes PCB layout of phased-array, FPGA code (Verilog HDL) for actuating elements, k-Wave projects for dataset collection, Python scripts for model training, and Arduino code for gas monitoring and XY-table actuation**

**1. FPGA Code**

We use two FPGA (Cyclone III) to control 196 elements. "MRI" is for the mian FPGA and "mjwu_ultrasound_array_secondary" is for the secondary FPGA.

**2. PCB desgin**

The "PCB artwork" is the Geber layout of PCB, and it can directly be used to manufacture PCB board.

**3. XY-axis motor platform**

"xyztable.ino" file is for controlling XY-axis motor platform

**4. Python for model training**

"mian.py" is for leaning model training using the PyCharm IDE.

**5. K-Wave model for data collection**

"data_collection.m" and "Copy_of_k_wave_multi_media_cuboid.m" are for dataset collection. The scripts are performed in Matlab.
   

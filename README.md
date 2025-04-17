**This repository includes all essential resources for replicating the proposed acoustic trapping, including PCB layout files for the phased-array hardware, Verilog HDL code for FPGA-based transducer actuation, k-Wave simulation projects for dataset generation, Python scripts for model training, and Arduino code for gas monitoring and XY-table control.**

**1. FPGA Code (Verilog HDL)**

We use two Cyclone III FPGAs to together control a 196-element phased-array transducer array.

The "MRI" is Verilog HDL code for the main FPGA, which handles timing control and ultrasound actuation logic.

The ultrasound_array_secondary folder contains the code for the secondary FPGA, which assists in synchronization and signal routing across channels.

**2. PCB desgin**

The "PCB artwork" is the Geber layout of PCB, and it can directly be used to manufacture PCB board.

**3. XY-axis motor platform**

"xyztable.ino" file is for controlling XY-axis motor platform

**4. Python for model training**

"mian.py" is for leaning model training using the PyCharm IDE.

**5. K-Wave model for data collection**

"data_collection.m" and "Copy_of_k_wave_multi_media_cuboid.m" are for dataset collection. The scripts are performed in Matlab.
   

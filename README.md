**This repository includes all essential resources for replicating the proposed acoustic trapping, including PCB layout files for the phased-array hardware, Verilog HDL code for FPGA-based transducer actuation, k-Wave simulation projects for dataset generation, Python scripts for model training, and Arduino code for gas monitoring and XY-table control.**

**1. FPGA Code (Verilog HDL)**

We use two Cyclone III FPGAs to together control a 196-element phased-array transducer array.

The "MRI.qpf" is Verilog HDL code for the main FPGA and the "ultrasound_array_secondary.qpf" is the code file for the secondary FPGA. Both handles timing control and ultrasound actuation logic.

**2. PCB desgin**

The "PCB artwork" includes the Gerber layout files required for manufacturing the custom phased-array driving board. These files are compatible with most PCB fabrication services and include silkscreen, soldermask, and drill layers.

**3. XY-axis motor platform**

"xyztable.ino" file contains the Arduino sketch for controlling the XY-axis motorized platform.

**4. Python for model training**

"mian.py" is the entry point for training the learning-based acoustic field prediction model. It is written in Python and tested using the PyCharm IDE. The script includes dataset loading, model definition, training pipeline, and evaluation steps.

**5. K-Wave model for data collection**

"data_collection.m" and "Copy_of_k_wave_multi_media_cuboid.m" are used to simulate acoustic wave propagation in heterogeneous media and collect time-of-flight (ToF) datasets. These scripts are executed in the MATLAB environment with the k-Wave toolbox installed.
   

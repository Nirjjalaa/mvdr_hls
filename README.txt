PROJECT: Hardware-Accelerated MVDR Adaptive Beamforming
         for 64-Element Antenna Arrays

AUTHOR:  Renu
GUIDE:   Sathyabama Sir
DURATION: 3 weeks

TOOLS:
- MATLAB/Octave (simulation)
- Vitis HLS (C++ to RTL)
- Versal FPGA (deployment)

PIPELINE:
MATLAB → C++ → Vitis HLS → Auto RTL → Versal FPGA

PARAMETERS:
- N = 64 (8x8 URA)
- d = 0.5 (half wavelength spacing)
- Target:      az=40°, el=0°
- Interferer1: az=30°, el=10°
- Interferer2: az=-60°, el=-25°

FILES:
- src/mvdr_top.cpp   → Vitis HLS top function
- tb/mvdr_tb.cpp     → Testbench for C simulation
- matlab/            → MATLAB/Octave simulation code

RESULTS:
- Target gain:  1.0000
- I1 rejection: -42.25 dB
- I2 rejection: -56.04 dB

MVDR Beamformer - Golden C++ Reference

Description
-----------
This repository contains the verified C++ implementation of the MVDR beamformer.
The implementation has been validated against the MATLAB reference model.

Verification Results
--------------------
Stage 1 : Steering Vector          PASS
Stage 2 : Covariance Matrix        PASS
Stage 3 : QR Decomposition         PASS
Stage 4 : Back Substitution        PASS
Stage 5 : Beamformer Gain          PASS

Accuracy
--------
Steering vector error        : 2.618e-4
QR reconstruction error      : 9.621e-16
Q orthogonality error        : 1.559e-16
Weight vector error          : 1.455e-12

Platform
--------
Compiler : g++
Language : C++17

# FRSICP: Fast and Robust Symmetric Iterative Closest Point Algorithm

## Overview
FRSICP is a MATLAB implementation of a **Fast and Robust Symmetric Iterative Closest Point (FRSICP)** algorithm designed for point cloud registration. This algorithm addresses limitations in the traditional Iterative Closest Point (ICP) algorithm, specifically improving convergence speed and robustness against noise, outliers, and partial overlaps. FRSICP is suitable for applications in computer vision, robotics, and 3D reconstruction.

## Features
- **Adaptive Robust Loss Function**: Increases resilience to noise and outliers.
- **Symmetric Objective Function**: Ensures robust rotation estimation based on normal vectors.
- **Anderson Acceleration**: Improves convergence speed, reducing computation time.
- **Efficient Linearization with Rodrigues Formula**: Optimizes small rotation increments for faster computation.

## Requirements
- MATLAB R2023b or newer
- Statistics and Machine Learning Toolbox (if applicable)
# FRSICP

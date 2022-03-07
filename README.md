# Lax-Friedrich-Scheme-Stability-Analysis
The Lax–Friedrichs method, named after Peter Lax and Kurt O. Friedrichs, is a numerical method for the solution of hyperbolic partial differential equations based on finite differences. The method can be described as the FTCS (forward in time, centered in space) scheme with a numerical dissipation term of 1/2.

# Workflow
1. Run lax-f.cpp - ``` g++ lax_f.cpp -o advection ```
2. Run the advection executable - ``` ./advection.exe ```
3. After this, you'll obtain dat files for the schemes you specified.
4. To obtain animations, you need to run animate.py or import your data to Google Drive and run the animation block in the Colab Notebook - (https://colab.research.google.com/drive/1F1glJJNm6vnMAyjQKtjErVDeTigqfUw1?usp=sharing).

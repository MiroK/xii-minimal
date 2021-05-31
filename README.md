# Minimal 3d1d with FEniCS_ii

- See Dockerfile for how to use docker/setup the minimal environment for using the code
- If you're using docker, once the container run navigate to `shared` and run `python3 daq_3d1d.py`.
  This demo used formulation from D'Angelo&Quarteroni and is more thoroughly documented
- `python3 laurino_lm.py` uses the 3d-1d formulation with a Lagrange multiplier
- FEniCS_ii source tree is included in the container and relevant functions for
  the reduction operators can be found in `/home/fenics/fenics_ii/xii/assembler/average_matrix.py`
# HP1paralogs-simulations
Scripts to run the simulations of HP1 paralogs

### 3D structures
The `PDB_files` folder contains all the PDB files used for the simulations in the paper (e.g., HP1 homo- and heterodimer, chimeras, ...).

### Single-chain simulations
1. All-atom simulations
The `AA_sim` folder contains the Amber99SBws-STQ force field and sample script to run an all-atom simulation of $\rm HP1\alpha$ homodimer using OpenMM (7.5.1).

```python
python run_equilibration.py -t 50 -f hp1a_dimer_amber99sbws_stq_tip4p2005s  # 50ns equilibration.
python run_omm.py -t 50 -f hp1a_dimer_amber99sbws_stq_tip4p2005s  -idx1 1 -idx2 1  # 350ns production from checkpoint.
```
2. Coarse-grained simulations
Sample scripts in the `CG_sim` folder are used to run the single-chain coarse-grained simulation of $\rm HP1\alpha$ homodimer in LAMMPS.

```bash
mpirun lmp_mpi -in hp1a.lmp > log,out
```

### Coarse-grained phase coexistence simulations

Sample scripts in the `CG_sim` folder are used to run the coarse-grained slab simulation of $\rm HP1\alpha$ homodimer in HOOMD.

```python
python run_nvt.py 170x170x170_box_min.gsd --mode=gpu
```

### Reference:

Phan et al., Interplay between charge distribution and DNA in shaping HP1 paralog phase separation and localization (eLife 2024), [https://doi.org/10.7554/eLife.90820]

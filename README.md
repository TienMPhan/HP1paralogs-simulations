# HP1paralogs-simulations
Scripts to run the simulations of HP1 paralogs

### 3D structures
The `PDB_files` folder contains all the PDB files used for the simulations in the paper (e.g., HP1 homo- and heterodimer, chimeras, ...).

### Single-chain simulations
1. All-atom simulations
The `AA_sim` folder contains the Amber99SBws-STQ force field and sample script to run an all-atom simulation of HP1$\alpha$ using OpenMM (7.5.1).

```python
python run_equilibration -t 50 -f hp1a_dimer_amber99sbws_stq_tip4p2005s  # 50ns equilibration.
python run_omm -t 50 -f hp1a_dimer_amber99sbws_stq_tip4p2005s  -idx1 1 -idx2 1  # 350ns production from checkpoint.
```

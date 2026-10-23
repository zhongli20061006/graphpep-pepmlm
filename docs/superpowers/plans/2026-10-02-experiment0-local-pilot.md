# Experiment 0: local RTX 5060 pilot

User approved trying AFM locally after the two-candidate cloud inputs were delivered.
Use one existing candidate, one AFM v3 model, one seed and the same NCAM1 20–718 receptor (709 aa complex).

- [x] Inspect driver, GPU memory, WSL memory, disk space and existing Python environment.
- [x] Create isolated `/home/g/colabfold_local` environment; preserve existing packages and driver.
- [ ] Install ColabFold 1.6.3 and CUDA 12.8-compatible JAX 0.6.2 dependencies.
- [ ] Verify GPU detection, compiled matrix computation and AFM import.
- [ ] Obtain official AFM v3 parameters and an online MSA.
- [ ] Run one real complex with constrained MSA depth, three recycles and no relaxation.
- [ ] Record runtime, resource usage, confidence, failure evidence and actual parameter count.
- [ ] Import the real result and run existing GraphPep if a valid structure is produced.
- [ ] Deliver an evidence-based local feasibility result and reproducible instructions.

Hardware observed: RTX 5060 Laptop, 8151 MiB total VRAM, driver 573.22, approximately 7.3 GiB free VRAM.
WSL initially provides about 7.6 GiB RAM. An isolated pip CUDA runtime avoids replacing the system toolkit.
Cloud inputs and results remain separate from this local trial.

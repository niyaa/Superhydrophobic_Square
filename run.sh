 mpirun -np 1 ~/../sgepner/Projects/nektar/build_hyp/dist/bin/ADRSolver geom.xml base.xml
 mv geom.fld geomB2D.fld
 mpirun -np 1 ~/../sgepner/Projects/nektar/build_hyp/dist/bin/IncNavierStokesSolver geom.xml base3d.xml
 mv geom.fld geom.bse
 mpirun -np 1 ~/../sgepner/Projects/nektar/build_hyp/dist/bin/IncNavierStokesSolver geom.xml stab.xml

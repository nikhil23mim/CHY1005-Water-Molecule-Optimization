# CHY1005 Project: Water Molecule Geometry and Energy Optimization Using Python

## Aim
To demonstrate molecular geometry and energy minimization for H2O using a Python application.

## Syllabus connection
This project supports the syllabus hands-on topics of molecular visualization, energy minimization, geometry optimization of polyatomic molecules, potential-energy concepts, and scientific simulation/modeling.

## Method
The application represents one oxygen and two hydrogen atoms in Cartesian coordinates. It calculates the two O-H bond lengths and H-O-H angle. A simplified harmonic potential penalizes deviations from reference values. Numerical finite differences estimate the energy gradient, and gradient descent with backtracking seeks a lower-energy geometry.

## Model and assumptions
Reference targets: O-H = 0.9572 Å and H-O-H = 104.52°. Illustrative force constants are used for the bond and angle terms. The energy is an educational model energy, not an absolute molecular energy. The model is not a quantum-chemical calculation and is not a validated production force field.

## Results
Run the Python application and copy the actual values shown in its results table. Insert a screenshot of the 3D molecule and energy graph. Do not label the values as experimental or ab-initio results.

| Parameter | Initial | Optimized | Unit |
|---|---:|---:|---|
| O-H bond 1 | From application | From application | Å |
| O-H bond 2 | From application | From application | Å |
| H-O-H angle | From application | From application | degrees |
| Model potential energy | From application | From application | illustrative kcal/mol |

## Discussion
Discuss whether the geometry moved toward the selected reference targets and whether the illustrative model energy decreased. The result depends on the chosen potential and optimization settings. The application is designed for teaching the optimization workflow rather than predicting all physical properties of water.

## Conclusion
This project demonstrates how computational chemistry software can represent a molecule, calculate geometric descriptors, evaluate a model potential, and optimize coordinates. It connects molecular geometry with potential energy and makes the optimization process visible through a 3D structure, results table, and energy graph.

## References
- Atkins, P. and de Paula, J., *Physical Chemistry*, Oxford University Press, 2011.
- Leach, A. R., *Molecular Modelling: Principles and Applications*.
- CHY1005 Introduction to Computational Chemistry syllabus supplied for this project.

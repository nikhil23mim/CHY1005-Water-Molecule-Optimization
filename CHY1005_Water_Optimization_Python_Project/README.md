# CHY1005 Computational Chemistry Project
## Water Molecule Geometry and Energy Optimization

### What this application does
- Displays a 3D ball-and-stick style water molecule.
- Starts from a deliberately distorted H-O-H geometry.
- Performs numerical-gradient optimization using a simple harmonic bond/angle model.
- Displays initial and optimized bond lengths, bond angle, and model energy.
- Plots the model energy over accepted optimization steps.
- Exports the results table as CSV.

### Important scientific limitation
This is an **educational illustrative model**, not a quantum-chemistry package and not an experimental measurement. It uses a simplified harmonic potential with chosen reference targets (O-H = 0.9572 Å; H-O-H = 104.52°) and illustrative force constants. Energy values are model energies relative to the reference geometry, not absolute molecular energies. In your report and presentation, identify them as simulated/illustrative.

### Requirements
- Python 3.9 or newer recommended
- NumPy
- Matplotlib
- Tkinter (usually included with standard Python on Windows; on some Linux systems install the OS package `python3-tk`)

### Install and run
Open a terminal in this folder and run:

```bash
python -m pip install -r requirements.txt
python water_optimization_app.py
```

On some systems, use `python3` instead of `python`.

### How to use
1. Start the application.
2. Review the distorted initial water geometry.
3. Click **Run geometry optimization**.
4. Compare the table values and energy plot.
5. Click **Save results as CSV** to export your run.

### Project workflow
1. Build a starting H2O geometry.
2. Calculate bond lengths and bond angle.
3. Evaluate the illustrative potential energy.
4. Estimate energy gradients numerically.
5. Use gradient descent with a backtracking step-size search.
6. Compare initial and optimized geometry and model energy.

### Suggested report statement
"The program uses a simplified harmonic potential to demonstrate the principle of geometry optimization. The values are illustrative simulation outputs and should not be interpreted as results from a validated water force field, density functional theory, or experiment."

"""
CHY1005 Computational Chemistry Project
Water Molecule Geometry and Energy Optimization

IMPORTANT SCIENTIFIC NOTE
This educational application uses a simple illustrative harmonic molecular potential.
It is NOT a quantum-chemical calculation and does NOT report experimentally measured
or ab-initio energies. Energy values are model energies in kcal/mol relative to the
chosen equilibrium geometry. The model is intended to demonstrate optimization.
"""

import math
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Reference geometry (angstrom and degrees)
R0 = 0.9572
THETA0_DEG = 104.52
THETA0 = math.radians(THETA0_DEG)

# Illustrative force constants (not a validated water force field)
K_BOND = 450.0       # kcal mol^-1 angstrom^-2, illustrative
K_ANGLE = 55.0       # kcal mol^-1 radian^-2, illustrative


def initial_coordinates():
    """Return a deliberately distorted H-O-H geometry in angstrom."""
    theta = math.radians(125.0)
    r1, r2 = 1.18, 0.82
    oxygen = np.array([0.0, 0.0, 0.0])
    h1 = np.array([r1, 0.0, 0.0])
    h2 = np.array([r2 * math.cos(theta), r2 * math.sin(theta), 0.0])
    return np.array([oxygen, h1, h2], dtype=float)


def geometry(coords):
    """Return O-H bond lengths and H-O-H angle."""
    o, h1, h2 = coords
    v1 = h1 - o
    v2 = h2 - o
    r1 = float(np.linalg.norm(v1))
    r2 = float(np.linalg.norm(v2))
    cosine = np.dot(v1, v2) / max(r1 * r2, 1e-12)
    angle = math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))
    return r1, r2, angle


def potential_energy(coords):
    """
    Illustrative harmonic energy:
    E = 1/2*k_b*(r1-r0)^2 + 1/2*k_b*(r2-r0)^2
        + 1/2*k_theta*(theta-theta0)^2
    Angle difference is in radians. Energy is relative to model equilibrium.
    """
    r1, r2, angle_deg = geometry(coords)
    angle_rad = math.radians(angle_deg)
    e_bond = 0.5 * K_BOND * ((r1 - R0) ** 2 + (r2 - R0) ** 2)
    e_angle = 0.5 * K_ANGLE * (angle_rad - THETA0) ** 2
    return float(e_bond + e_angle)


def numerical_gradient(coords, delta=1e-5):
    grad = np.zeros_like(coords)
    for i in range(coords.shape[0]):
        for j in range(3):
            plus = coords.copy()
            minus = coords.copy()
            plus[i, j] += delta
            minus[i, j] -= delta
            grad[i, j] = (potential_energy(plus) - potential_energy(minus)) / (2 * delta)
    return grad


def optimize_geometry(start, max_steps=250, tolerance=1e-5):
    """
    Gradient descent with backtracking line search.
    Returns optimized coordinates and energy history.
    """
    coords = start.copy()
    history = [potential_energy(coords)]
    step_scale = 0.03

    for _ in range(max_steps):
        grad = numerical_gradient(coords)
        norm = float(np.linalg.norm(grad))
        if norm < tolerance:
            break

        direction = -grad / max(norm, 1e-12)
        current_e = potential_energy(coords)
        step = step_scale
        accepted = False

        for _trial in range(25):
            candidate = coords + step * direction
            candidate_e = potential_energy(candidate)
            if candidate_e < current_e:
                coords = candidate
                history.append(candidate_e)
                step_scale = min(step * 1.25, 0.08)
                accepted = True
                break
            step *= 0.5

        if not accepted:
            # Try a smaller step next iteration; stop if no useful movement.
            step_scale *= 0.5
            if step_scale < 1e-8:
                break

        if len(history) > 2 and abs(history[-2] - history[-1]) < 1e-12:
            break

    return coords, history


class WaterOptimizationApp:
    def __init__(self, root):
        self.root = root
        self.root.title("CHY1005 | Water Molecule Geometry & Energy Optimization")
        self.root.geometry("1120x760")
        self.root.minsize(900, 650)

        self.start = initial_coordinates()
        self.final = self.start.copy()
        self.history = [potential_energy(self.start)]

        self._build_ui()
        self._refresh_all()

    def _build_ui(self):
        top = ttk.Frame(self.root, padding=(12, 10))
        top.pack(fill="x")
        ttk.Label(top, text="Water Molecule Geometry & Energy Optimization",
                  font=("TkDefaultFont", 17, "bold")).pack(anchor="w")
        ttk.Label(top, text="CHY1005 — Introduction to Computational Chemistry",
                  font=("TkDefaultFont", 10)).pack(anchor="w", pady=(3, 0))

        note = ttk.Label(
            self.root,
            text=("EDUCATIONAL MODEL: illustrative harmonic potential only; not a quantum-chemistry "
                  "or experimental result. Energy is a model value relative to the reference geometry."),
            foreground="#8a3b12", wraplength=1060
        )
        note.pack(fill="x", padx=14, pady=(0, 8))

        controls = ttk.Frame(self.root, padding=(12, 4))
        controls.pack(fill="x")
        ttk.Button(controls, text="Reset distorted starting geometry",
                   command=self.reset).pack(side="left", padx=(0, 8))
        ttk.Button(controls, text="Run geometry optimization",
                   command=self.run_optimization).pack(side="left", padx=(0, 8))
        ttk.Button(controls, text="Save results as CSV",
                   command=self.save_csv).pack(side="left")

        self.status = ttk.Label(controls, text="")
        self.status.pack(side="right")

        main = ttk.Frame(self.root, padding=(10, 4))
        main.pack(fill="both", expand=True)

        self.fig = Figure(figsize=(9, 5.5), dpi=100, constrained_layout=True)
        self.ax3d = self.fig.add_subplot(121, projection="3d")
        self.axenergy = self.fig.add_subplot(122)
        self.canvas = FigureCanvasTkAgg(self.fig, master=main)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        lower = ttk.LabelFrame(self.root, text="Results table", padding=8)
        lower.pack(fill="x", padx=12, pady=(2, 8))

        columns = ("parameter", "initial", "optimized", "unit")
        self.tree = ttk.Treeview(lower, columns=columns, show="headings", height=5)
        for col, title, width in [
            ("parameter", "Parameter", 280),
            ("initial", "Initial geometry", 210),
            ("optimized", "Optimized geometry", 220),
            ("unit", "Unit / meaning", 250),
        ]:
            self.tree.heading(col, text=title)
            self.tree.column(col, width=width, anchor="center")
        self.tree.pack(fill="x")

        conclusion = (
            "Interpretation: optimization should lower the illustrative model energy and move the "
            "geometry toward the model's reference O–H bond length (0.9572 Å) and H–O–H angle "
            "(104.52°). These reference geometry values are used as model targets; the calculated "
            "energy is not an absolute molecular energy."
        )
        ttk.Label(self.root, text=conclusion, wraplength=1080, justify="left",
                  padding=(12, 0, 12, 10)).pack(fill="x")

    def reset(self):
        self.start = initial_coordinates()
        self.final = self.start.copy()
        self.history = [potential_energy(self.start)]
        self.status.config(text="Starting geometry reset.")
        self._refresh_all()

    def run_optimization(self):
        self.final, self.history = optimize_geometry(self.start)
        self.status.config(text=f"Finished: {len(self.history)-1} accepted steps")
        self._refresh_all()

    def _refresh_all(self):
        self.ax3d.clear()
        self.axenergy.clear()

        # Plot optimized geometry (or starting geometry before running)
        coords = self.final
        labels = ["O", "H", "H"]
        colors = ["red", "white", "white"]
        sizes = [260, 170, 170]
        for idx, (point, label, color, size) in enumerate(zip(coords, labels, colors, sizes)):
            edge = "black"
            self.ax3d.scatter(point[0], point[1], point[2], s=size, c=color,
                              edgecolors=edge, depthshade=True)
            self.ax3d.text(point[0], point[1], point[2], " " + label, fontsize=11)
        for h_idx in [1, 2]:
            self.ax3d.plot([coords[0, 0], coords[h_idx, 0]],
                           [coords[0, 1], coords[h_idx, 1]],
                           [coords[0, 2], coords[h_idx, 2]], linewidth=3)
        self.ax3d.set_title("3D molecular structure (current)")
        self.ax3d.set_xlabel("X (Å)")
        self.ax3d.set_ylabel("Y (Å)")
        self.ax3d.set_zlabel("Z (Å)")
        self.ax3d.set_box_aspect((1, 1, 0.65))

        self.axenergy.plot(range(len(self.history)), self.history, marker="o", markersize=3)
        self.axenergy.set_title("Energy during optimization")
        self.axenergy.set_xlabel("Accepted optimization step")
        self.axenergy.set_ylabel("Illustrative model energy (kcal/mol)")
        self.axenergy.grid(True, alpha=0.3)

        r1i, r2i, ai = geometry(self.start)
        r1f, r2f, af = geometry(self.final)
        ei = potential_energy(self.start)
        ef = potential_energy(self.final)
        rows = [
            ("O–H bond 1", f"{r1i:.4f}", f"{r1f:.4f}", "angstrom (Å)"),
            ("O–H bond 2", f"{r2i:.4f}", f"{r2f:.4f}", "angstrom (Å)"),
            ("H–O–H bond angle", f"{ai:.3f}", f"{af:.3f}", "degrees"),
            ("Potential energy", f"{ei:.6f}", f"{ef:.6f}", "illustrative kcal/mol"),
            ("Energy change (final − initial)", "—", f"{ef-ei:.6f}", "illustrative kcal/mol"),
        ]
        for item in self.tree.get_children():
            self.tree.delete(item)
        for row in rows:
            self.tree.insert("", "end", values=row)
        self.canvas.draw_idle()

    def save_csv(self):
        from tkinter.filedialog import asksaveasfilename
        path = asksaveasfilename(
            title="Save water optimization results",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")]
        )
        if not path:
            return
        import csv
        r1i, r2i, ai = geometry(self.start)
        r1f, r2f, af = geometry(self.final)
        ei = potential_energy(self.start)
        ef = potential_energy(self.final)
        rows = [
            ("Parameter", "Initial", "Optimized", "Unit / note"),
            ("O-H bond 1", r1i, r1f, "angstrom"),
            ("O-H bond 2", r2i, r2f, "angstrom"),
            ("H-O-H angle", ai, af, "degrees"),
            ("Potential energy", ei, ef, "illustrative kcal/mol; not ab-initio"),
            ("Energy change", "", ef-ei, "illustrative kcal/mol; final minus initial"),
            ("Model", "", "", "harmonic bond-angle educational model"),
            ("Reference bond length", "", R0, "angstrom"),
            ("Reference angle", "", THETA0_DEG, "degrees"),
            ("Bond force constant", "", K_BOND, "illustrative model parameter"),
            ("Angle force constant", "", K_ANGLE, "illustrative model parameter"),
        ]
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerows(rows)
        messagebox.showinfo("Saved", f"Results saved to:\n{path}")


def main():
    root = tk.Tk()
    try:
        ttk.Style().theme_use("clam")
    except tk.TclError:
        pass
    WaterOptimizationApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()

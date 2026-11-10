#!/usr/bin/env python3
"""Plot STRATA-NOVA staged density and luminosity without hiding low-density CSM.

This code reads the stitched ``rho.xg`` and ``radius.xg`` pair by snapshot,
so density and radius always come from the same time record.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

DAY = 86400.0


def configure_style() -> None:
    mpl.rcParams["font.family"] = "STIXGeneral"
    mpl.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams["xtick.direction"] = "in"
    mpl.rcParams["ytick.direction"] = "in"
    mpl.rcParams["xtick.top"] = True
    mpl.rcParams["ytick.right"] = True
    mpl.rcParams["xtick.minor.visible"] = True
    mpl.rcParams["ytick.minor.visible"] = True
    mpl.rcParams["axes.linewidth"] = 1.0


def read_xg(path: Path) -> tuple[np.ndarray, list[np.ndarray], list[np.ndarray]]:
    times: list[float] = []
    x_profiles: list[np.ndarray] = []
    y_profiles: list[np.ndarray] = []
    current_x: list[float] = []
    current_y: list[float] = []
    current_time: float | None = None

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            line = raw.lstrip()
            if line.startswith('"Time ='):
                if current_time is not None:
                    times.append(current_time)
                    x_profiles.append(np.asarray(current_x, dtype=float))
                    y_profiles.append(np.asarray(current_y, dtype=float))
                current_time = float(line.split("=", 1)[1])
                current_x, current_y = [], []
                continue
            fields = line.split()
            if len(fields) < 2:
                continue
            try:
                current_x.append(float(fields[0]))
                current_y.append(float(fields[1]))
            except ValueError:
                continue

    if current_time is not None:
        times.append(current_time)
        x_profiles.append(np.asarray(current_x, dtype=float))
        y_profiles.append(np.asarray(current_y, dtype=float))
    return np.asarray(times), x_profiles, y_profiles



def read_key_value_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if "=" not in raw:
            continue
        key, value = raw.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def as_float(values: dict[str, str], key: str) -> float:
    try:
        return float(values[key].replace("D", "E").replace("d", "e"))
    except (KeyError, ValueError):
        return float("nan")


def plot_density(run_dir: Path, output: Path, days: list[float]) -> None:
    # Read each stitched XG file only once.  Long staged runs can contain many
    # snapshots, so repeatedly parsing the files for every requested day is
    # needlessly expensive.
    t_rho, _, rho_profiles = read_xg(run_dir / "rho.xg")
    t_rad, _, radius_profiles = read_xg(run_dir / "radius.xg")
    if t_rho.size == 0 or t_rad.size == 0:
        raise ValueError("No density/radius snapshots found")

    snapshots: list[tuple[float, np.ndarray, np.ndarray]] = []
    for target in days:
        distance = np.abs(t_rho / DAY - target)
        candidates = np.flatnonzero(np.isclose(distance, np.nanmin(distance), rtol=0.0, atol=1e-12))
        i_rho = int(candidates[-1])
        i_rad = int(np.argmin(np.abs(t_rad - t_rho[i_rho])))
        rho = np.asarray(rho_profiles[i_rho], dtype=float)
        radius = np.asarray(radius_profiles[i_rad], dtype=float)
        n = min(rho.size, radius.size)
        rho, radius = rho[:n], radius[:n]
        valid = np.isfinite(rho) & np.isfinite(radius) & (rho > 0.0) & (radius > 0.0)
        snapshots.append((float(t_rho[i_rho] / DAY), radius[valid], rho[valid]))

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    for actual, radius, rho in snapshots:
        order = np.argsort(radius)
        ax.plot(np.log10(radius[order]), np.log10(rho[order]), lw=1.15, label=f"{actual:.3f} d")

    ax.set_xlabel(r"$\log_{10}(R/\mathrm{cm})$")
    ax.set_ylabel(r"$\log_{10}(\rho/\mathrm{g\,cm^{-3}})$")
    ax.set_title("STRATA-NOVA staged density evolution")
    ax.grid(True, which="major", alpha=0.22, linestyle="--")
    ax.grid(True, which="minor", alpha=0.10, linestyle=":")
    ax.legend(loc="best", fontsize=8)
    positive = np.concatenate([rho for _, _, rho in snapshots]) if snapshots else np.array([])
    if positive.size:
        lo = max(-31.0, np.floor(np.log10(np.min(positive))) - 0.5)
        hi = np.ceil(np.log10(np.max(positive))) + 0.5
        ax.set_ylim(lo, hi)
    fig.tight_layout()
    fig.savefig(output, dpi=500)
    plt.close(fig)


def plot_attachment_profiles(run_dir: Path, output: Path) -> None:
    """Plot the generated density profile at each finite-shell attachment.

    This is the correct place to check whether ``constant_density`` is actually
    flat.  Once a shell has been evolved, shock compression is expected to
    destroy the rectangular initial shape.
    """
    parent = run_dir.parent
    profiles = sorted(parent.glob("stage*_shell/generated_csm/wdprofile_generated.dat"))
    if not profiles:
        return

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    all_rho: list[np.ndarray] = []
    for index, profile in enumerate(profiles, start=1):
        data = np.loadtxt(profile, skiprows=1, ndmin=2)
        radius = np.asarray(data[:, 2], dtype=float)
        rho = np.asarray(data[:, 4], dtype=float)
        valid = np.isfinite(radius) & np.isfinite(rho) & (radius > 0.0) & (rho > 0.0)
        radius, rho = radius[valid], rho[valid]
        order = np.argsort(radius)
        radius, rho = radius[order], rho[order]
        all_rho.append(rho)

        diagnostics = read_key_value_file(profile.parent / "boundary_diagnostics.txt")
        r_in = as_float(diagnostics, "attachment_radius_cm")
        r_out = as_float(diagnostics, "outer_radius_cm")
        expected = as_float(diagnostics, "outer_density_normalization")
        # The diagnostic radii are formatted to six significant figures, so
        # selecting by those rounded bounds can accidentally include the old
        # surface zone.  ``extension_zones`` identifies the appended shell
        # exactly and is therefore the robust flatness test.
        try:
            n_shell = int(float(diagnostics.get("extension_zones", "0")))
        except ValueError:
            n_shell = 0
        shell = np.zeros(radius.size, dtype=bool)
        if 0 < n_shell <= radius.size:
            shell[-n_shell:] = True
        elif np.isfinite(r_in) and np.isfinite(r_out):
            shell = (radius > r_in) & (radius <= r_out)
        if np.any(shell):
            shell_rho = rho[shell]
            measured = float(np.median(shell_rho))
            rel_scatter = float(np.max(np.abs(shell_rho - measured)) / max(abs(measured), 1.0e-300))
            mass = as_float(diagnostics, "csm_mass_msun")
            print(
                f"shell {index}: R=[{r_in:.6e}, {r_out:.6e}] cm, "
                f"M={mass:.6e} Msun, rho_expected={expected:.6e}, "
                f"rho_median={measured:.6e}, max_fractional_scatter={rel_scatter:.3e}"
            )
            # Shade only the newly appended finite shell.  The profile itself
            # remains unsmoothed and is plotted exactly as written by STRATA-NOVA.
            ax.axvspan(np.log10(r_in), np.log10(r_out), alpha=0.07)

        ax.plot(
            np.log10(radius), np.log10(rho), lw=1.25, drawstyle="steps-post",
            label=f"Shell {index} attachment",
        )

    ax.set_xlabel(r"$\log_{10}(R/\mathrm{cm})$")
    ax.set_ylabel(r"$\log_{10}(\rho/\mathrm{g\,cm^{-3}})$")
    ax.set_title("STRATA-NOVA finite-shell attachment profiles (new shell shaded)")
    ax.grid(True, which="major", alpha=0.22, linestyle="--")
    ax.grid(True, which="minor", alpha=0.10, linestyle=":")
    ax.legend(loc="best", fontsize=8)
    if all_rho:
        merged = np.concatenate(all_rho)
        lo = max(-31.0, np.floor(np.log10(np.min(merged))) - 0.5)
        hi = np.ceil(np.log10(np.max(merged))) + 0.5
        ax.set_ylim(lo, hi)
    fig.tight_layout()
    fig.savefig(output, dpi=500)
    plt.close(fig)


def load_luminosity(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, usecols=(0, 1), ndmin=2)
    valid = np.isfinite(data[:, 0]) & np.isfinite(data[:, 1]) & (data[:, 1] > 0.0)
    return data[valid, 0] / DAY, data[valid, 1]


def plot_luminosity(run_dir: Path, output: Path, hesma: Path | None) -> None:
    time, lum = load_luminosity(run_dir / "lum_observed.dat")
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.plot(time, np.log10(lum), lw=1.2, label="STRATA-NOVA staged")
    if hesma is not None and hesma.is_file():
        h = np.loadtxt(hesma, usecols=(0, 1), ndmin=2)
        valid = np.isfinite(h[:, 0]) & np.isfinite(h[:, 1]) & (h[:, 1] > 0.0)
        ax.plot(h[valid, 0], np.log10(h[valid, 1]), ls="None", marker="o", markersize=3,
                label="HESMA n1")
    ax.set_xlim(0.0, max(16.0, float(np.nanmax(time))))
    ax.set_xlabel("Time since explosion (days)")
    ax.set_ylabel(r"$\log_{10}(L/\mathrm{erg\ s^{-1}})$")
    ax.set_title("STRATA-NOVA grey bolometric luminosity with HESMA reference")
    ax.grid(True, which="major", alpha=0.22, linestyle="--")
    ax.grid(True, which="minor", alpha=0.10, linestyle=":")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=500)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path, help="staged final directory (contains rho.xg and lum_observed.dat)")
    parser.add_argument("--out", type=Path, default=Path("Analysis/results"))
    parser.add_argument("--hesma", type=Path, default=None, help="optional HESMA *_lightcurve.dat")
    parser.add_argument("--density-days", type=float, nargs="+", default=[2.0, 7.0, 7.5, 16.0])
    args = parser.parse_args()
    configure_style()
    args.out.mkdir(parents=True, exist_ok=True)
    plot_density(args.run_dir, args.out / "staged_density_smoke.png", args.density_days)
    plot_attachment_profiles(args.run_dir, args.out / "staged_attachment_density.png")
    plot_luminosity(args.run_dir, args.out / "staged_luminosity_hesma.png", args.hesma)


if __name__ == "__main__":
    main()

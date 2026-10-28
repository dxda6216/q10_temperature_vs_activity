#!/usr/bin/env python3
"""
================================================================================
Q10 calculator for temperature vs. activity - desktop (Python/Tkinter) version
================================================================================
Fits activity(T) = A_b * Q10 ** ((T - T_b) / 10) to temperature/activity data
and shows the fit as Activity-vs-Temperature, ln(Activity)-vs-Temperature and
Arrhenius plots, with an optional Arrhenius regression line.

Requirements: Python 3.9+, numpy, scipy, matplotlib (tkinter ships with Python)
Run:          python q10_temperature_vs_activity.py

A browser version of this calculator (single HTML file, no installation) is
provided as index.html.

==============================================================================
1. Background
==============================================================================
[Q10 (temperature coefficient)]
  The factor by which a rate R changes for a 10 °C rise in temperature.
  From rates R1, R2 measured at temperatures T1, T2:
      Q10 = (R2 / R1) ^ (10 / (T2 - T1))
  Typical enzyme-catalysed reactions have Q10 of about 2 to 3.

  Q10 > 1 : activity rises with temperature
  Q10 = 1 : activity does not depend on temperature
  Q10 < 1 : activity falls as temperature rises

[Model that is fitted]
  With A_b the activity at the base temperature T_b:
      A(T) = A_b * Q10 ^ ((T - T_b) / 10)                       ... (1)
  Taking logarithms:
      ln A(T) = ln A_b + ((T - T_b) / 10) * ln Q10              ... (2)
  So ln A is a straight line in T with slope ln(Q10) / 10.

[Meaning of the base temperature T_b]
  T_b is only a reparameterisation: it decides at which temperature the
  fitted activity A_b is reported. Q10, its standard error, R^2 and the
  fitted curve itself are identical for every choice of T_b; only A_b and
  its error change.

[Arrhenius plot]
  The Arrhenius equation k = A * exp(-Ea / (R * T_K)) gives a straight line
  when ln k is plotted against 1/T_K (T_K = absolute temperature), with
  slope -Ea/R (Ea: activation energy, R = 8.314 J mol^-1 K^-1).
  This program plots ln(activity) against 1000/T_K (K^-1). Model (2) is
  linear in T, not in 1/T, so the fitted Q10 curve is slightly bent on the
  Arrhenius plot. The two descriptions are approximately related by
      Ea ~= R * T1 * T2 * ln(Q10) / 10        (T1, T2 in kelvin)

[Arrhenius regression line]
  Optionally (checkbox, on by default) the Arrhenius plot also shows the
  ordinary least-squares straight line of ln(activity) against 1000/T_K,
  fitted to the same data used for Q10:
      ln A = intercept + slope * (1000 / T_K),   Ea = -slope * 1000 * R

==============================================================================
2. Calculation method
==============================================================================
  (a) Non-linear least squares (scipy.optimize.curve_fit) minimising the sum
      of squared residuals of the activity itself (not its logarithm):
          SSE = sum_i (A_i - A_hat(T_i))^2
  (b) Starting values: linear regression of ln(A) on T, converted to A_b and
      Q10 through relation (2).
  (c) Standard errors: curve_fit's covariance matrix (absolute_sigma=False),
      i.e. Cov = s^2 * (J^T J)^-1 with s^2 = SSE / (n - 2). The "±" values
      are the square roots of its diagonal.
  (d) Coefficient of determination: R^2 = 1 - SSE / SST,
      SST = sum_i (A_i - mean(A))^2
  (e) Fitted curve: the full temperature range of the data, extended on each
      side by 33.3 % of its width (end points truncated to integers), drawn
      with 200 points. When a manual temperature axis is wider, the curve is
      drawn across that range too.
  (f) Arrhenius regression: x_i = 1000 / (T_i + 273.15), y_i = ln(A_i),
      simple linear regression y = intercept + slope * x on the fitted data,
      SE(slope) = sqrt(s^2 / Sxx), SE(intercept) = sqrt(s^2 (1/n + mean(x)^2/Sxx)),
      Ea = -slope * 1000 * R, SE(Ea) = SE(slope) * 1000 * R.

==============================================================================
3. Axis scales
==============================================================================
  The axis scales follow the all-data fit automatically. To set them by hand,
  tick "Temperature" and/or "Activity" under "Axis scales" and enter Min and
  Max; press Enter, click the arrows or move to another field to apply.
    Temperature : 0 to 100 °C, step 0.5 °C. The Arrhenius x axis
                  (1000/T) is converted from this range.
    Activity    : 0 to 100 x the largest activity in the data. The step of
                  the arrows is set from the data: (data maximum / 50)
                  rounded to 1, 2 or 5 x 10^k. The ln(Activity) and Arrhenius
                  y axes are converted (ln); with Min = 0 their lower limit
                  stays automatic, because ln(0) is undefined.

==============================================================================
4. How to use
==============================================================================
  1) On start-up the window shows a demo calculation with example data.
  2) Click "Clear" to empty the form.
  3) Enter one "temperature(°C) activity" pair per line (commas, spaces, tabs
     or semicolons), or use "Load CSV..." (temperature in column 1, activity
     in column 2). Optionally enter the activity unit.
  4) Optionally tick the temperature-range box and set the range used for the
     fit, or drag across any plot after a calculation. Double-click a plot or
     click "Use all data" to return to all data.
  5) Choose the base temperature and click "Calculate" (or Ctrl+Enter).
  6) "Save figure..." saves the three plots stacked in one image (PNG, PDF or
     SVG); "Save results..." saves the results text.

==============================================================================
5. Privacy
==============================================================================
  All calculations run locally. The program makes no network connections and
  stores nothing; files are written only when you use the Save buttons.
================================================================================
"""

import csv
import math
import re
import statistics
import tkinter as tk
import tkinter.font as tkfont
from tkinter import filedialog, messagebox, ttk

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector
from scipy.optimize import curve_fit

# ======================================================================
# Constants
# ======================================================================
DEG = "\u00B0"
R_GAS = 8.314462618          # gas constant, J mol^-1 K^-1
DATA_ROWS = 14               # visible lines in the data input box
TITLE_SIZE = 18              # plot-title font size (single-plot tabs)
TITLE_SIZE_COMBINED = 20     # plot-title font size (saved 3-plot figure)

# Plot colours
C_USED, C_EXCL, C_FIT, C_ARR = "#FF0000", "#5DC85D", "#2B8CD6", "dimgray"
C_BAND = "#6F8BAE"           # fit-range band and drag selection
# Drawing order: fitted lines (Z_LINE) below the data points (Z_POINT), so lines never cover points.
Z_LINE, Z_POINT = 2, 3

# Demo data, loaded and calculated at start-up. Temperature in °C, activity in U/mg.
EXAMPLE_PAIRS = [(10, 14.8), (12.5, 18.9), (15, 23.4), (17.5, 27.1), (20, 32.9), (22.5, 42.3),
                 (25, 49.1), (27.5, 59.8), (30, 75.6), (32.5, 88.7), (35, 111.9), (37.5, 131.2),
                 (40, 166.4)]
EXAMPLE_DATA = "\n".join(f"{t}, {a}" for t, a in EXAMPLE_PAIRS)
EXAMPLE_UNIT = "U/mg"

CUSTOM_BASE = "Custom value (enter below)"
BASE_CHOICES = [f"0{DEG}C", f"4{DEG}C", f"20{DEG}C", f"25{DEG}C", f"30{DEG}C", f"37{DEG}C",
                "Minimum", "Maximum", "Average", "Median", CUSTOM_BASE]
FIXED_BASE = {BASE_CHOICES[0]: 0.0, BASE_CHOICES[1]: 4.0, BASE_CHOICES[2]: 20.0,
              BASE_CHOICES[3]: 25.0, BASE_CHOICES[4]: 30.0, BASE_CHOICES[5]: 37.0}

# Manual axis scales
T_AXIS_RANGE = (0.0, 100.0)  # °C
T_AXIS_STEP = 0.5            # °C
A_AXIS_FACTOR = 100          # activity axis may go up to 100 x the data maximum

NOTICE_DEMO = ("You're looking at a demo calculation using example data. To analyse your own "
               "data, click Clear and enter your temperature and activity pairs in the form.")


# ======================================================================
# Calculation core (no GUI dependencies)
# Flow: parse_pairs() -> analyze() -> fit_q10() and arrhenius_fit()
# ======================================================================
def parse_pairs(text):
    """Parse one 'temperature activity' pair per line.

    - separators: comma / space / tab / semicolon (may be mixed)
    - blank lines and lines starting with '#' are ignored
    - one leading non-numeric line is skipped as a header
    Returns (temperatures, activities) as numpy arrays.
    """
    temps, acts = [], []
    header_skipped = False
    for n, line in enumerate(text.splitlines(), start=1):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        tokens = [t for t in re.split(r"[,;\s]+", s) if t]
        try:
            vals = [float(t) for t in tokens]
        except ValueError:
            if not temps and not header_skipped:
                header_skipped = True
                continue
            raise ValueError(f"Line {n}: cannot read numbers: {s!r}") from None
        if len(vals) != 2:
            raise ValueError(f"Line {n}: expected 2 values (temperature, activity), "
                             f"found {len(vals)}: {s!r}")
        temps.append(vals[0])
        acts.append(vals[1])
    return np.array(temps, dtype=float), np.array(acts, dtype=float)


def resolve_base_temperature(selection, custom_value, xa):
    """Base temperature T_b as a number. Minimum/Maximum/Average/Median come from
    xa, the data USED in the fit."""
    if selection in FIXED_BASE:
        return FIXED_BASE[selection]
    if selection == "Minimum":
        return float(np.min(xa))
    if selection == "Maximum":
        return float(np.max(xa))
    if selection == "Average":
        return float(statistics.mean(xa))
    if selection == "Median":
        return float(statistics.median(xa))
    return float(custom_value)


def make_q10_model(base_t):
    """Return model (1), A(T) = A_b * Q10^((T - T_b)/10), for a given T_b."""
    def q10_model(temperature, a_b, q10):
        return a_b * q10 ** ((temperature - base_t) * 0.1)
    return q10_model


def initial_guess(xa, ya, base_t):
    """Starting values from ln(A) = c + s*T:  ln Q10 = 10 s,  ln A_b = c + s*T_b."""
    try:
        slope, intercept = np.polyfit(xa, np.log(ya), 1)
        p0 = (float(np.exp(intercept + slope * base_t)), float(np.exp(10.0 * slope)))
        if all(np.isfinite(p0)) and p0[0] > 0 and p0[1] > 0:
            return p0
    except Exception:
        pass
    return (float(np.mean(ya)), 1.0)


def fit_q10(xa, ya, base_t):
    """Non-linear least-squares fit of model (1). Returns (popt, pcov, R^2, model)."""
    model = make_q10_model(base_t)
    p0 = initial_guess(xa, ya, base_t)
    popt, pcov = curve_fit(model, xa, ya, p0=p0, maxfev=20000)
    residuals = ya - model(xa, *popt)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((ya - np.mean(ya)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return popt, pcov, r2, model


def arrhenius_fit(temps, acts):
    """Least-squares line of ln(activity) on 1000/T_K; Ea = -slope * 1000 * R."""
    t, a = np.asarray(temps, float), np.asarray(acts, float)
    x, y = 1000.0 / (t + 273.15), np.log(a)
    n, mx, my = len(x), x.mean(), y.mean()
    sxx = np.sum((x - mx) ** 2)
    syy = np.sum((y - my) ** 2)
    slope = np.sum((x - mx) * (y - my)) / sxx
    intercept = my - slope * mx
    sse = np.sum((y - intercept - slope * x) ** 2)
    s2 = sse / (n - 2)
    se_slope = math.sqrt(s2 / sxx)
    se_int = math.sqrt(s2 * (1.0 / n + mx * mx / sxx))
    return dict(
        slope=slope, intercept=intercept, se_slope=se_slope, se_int=se_int,
        r2=1 - sse / syy if syy > 0 else float("nan"),
        Ea=-slope * 1000 * R_GAS, se_Ea=se_slope * 1000 * R_GAS,   # J/mol
    )


def arrhenius_lna(arr, temps_c):
    """ln(activity) predicted by the Arrhenius regression for temperatures in °C."""
    return arr["intercept"] + arr["slope"] * 1000.0 / (np.asarray(temps_c, float) + 273.15)


def analyze(temps, acts, use_range, low, high, base_sel, base_custom):
    """Run the whole calculation; xa/ya are the data used for Q10, x/y all data."""
    x, y = np.asarray(temps, float), np.asarray(acts, float)
    if len(x) != len(y):
        raise ValueError(f"Number of temperatures ({len(x)}) and activities ({len(y)}) differ.")
    if len(x) < 3:
        raise ValueError("Enter at least 3 data points.")
    if np.any(y <= 0):
        raise ValueError("All activities must be positive (logarithms are used).")
    if np.any(x <= -273.15):
        raise ValueError("Temperatures must be above absolute zero (-273.15 °C).")
    if use_range:
        if not low < high:
            raise ValueError("Range low limit must be smaller than the high limit.")
        mask = (x >= low) & (x <= high)
        xa, ya = x[mask], y[mask]
    else:
        xa, ya = x, y
    if len(xa) < 3:
        raise ValueError(f"Only {len(xa)} data point(s) in the selected range; need at least 3.")

    base_t = resolve_base_temperature(base_sel, base_custom, xa)
    popt, pcov, r2, model = fit_q10(xa, ya, base_t)
    perr = np.sqrt(np.abs(np.diag(pcov)))

    margin = (x.max() - x.min()) * 0.333
    span = (int(x.min() - margin), int(x.max() + margin) + 1)
    fcx = np.linspace(span[0], span[1], 200)
    fcy = model(fcx, *popt)
    arr = arrhenius_fit(xa, ya)

    # Reference curves fitted to ALL data: the automatic axis limits come from these,
    # so restricting the fit range never changes the axis scales.
    ref = dict(fcy=fcy, arr=arr)
    if use_range:
        try:
            popt_all, _, _, model_all = fit_q10(x, y, base_t)
            ref = dict(fcy=model_all(fcx, *popt_all), arr=arrhenius_fit(x, y))
        except (ValueError, RuntimeError):
            pass

    return dict(x=x, y=y, xa=xa, ya=ya, use_range=use_range, low=low, high=high, base_t=base_t,
                popt=popt, perr=perr, r2=r2, model=model, span=span, fcx=fcx, fcy=fcy,
                arr=arr, ref=ref, unit="")


# ---------------------------------------------------------------------
# Activity-axis scale: 0 .. 100 x data maximum, step derived from the data
# ---------------------------------------------------------------------
def nice_step(raw):
    """Round a positive number to 1, 2 or 5 x 10^k."""
    mag = 10 ** math.floor(math.log10(raw))
    nr = raw / mag
    return (1 if nr < 1.5 else 2 if nr < 3.5 else 5 if nr < 7.5 else 10) * mag


def activity_axis_scale(max_a):
    """Allowed range and step of the manual activity axis.

    step  = (data maximum / 50) rounded to 1, 2 or 5 x 10^k
    limit = 100 x data maximum, rounded down to a multiple of step
    dec   = decimals needed to show the step
    """
    step = nice_step(max_a / 50.0)
    dec = max(0, -math.floor(math.log10(step) + 1e-9))
    limit = round(math.floor(A_AXIS_FACTOR * max_a / step + 1e-9) * step, dec)
    return dict(max=max_a, step=step, dec=dec, limit=limit)


# ======================================================================
# Report text
# ======================================================================
def fmt_a(v):
    """Activity values of any magnitude."""
    if not np.isfinite(v):
        return str(v)
    a = abs(v)
    return f"{v:.3e}" if a != 0 and (a >= 1e5 or a < 1e-2) else f"{v:.3f}"


def stats_text(arr, name, unit, f):
    u = f" {unit}" if unit else ""
    return (f"Minimum {name} = {f(np.min(arr))}{u}\n"
            f"Maximum {name} = {f(np.max(arr))}{u}\n"
            f"Average {name} = {f(statistics.mean(arr))}{u}\n"
            f"Median {name} = {f(statistics.median(arr))}{u}\n")


def table_text(x, y):
    lines = [f"{'#':>4} {'Temperature':>12} {'Activity':>12}"]
    for i, (xi, yi) in enumerate(zip(x, y)):
        lines.append(f"{i:>4} {xi:>12.3f} {fmt_a(yi):>12}")
    return "\n".join(lines)


def summary_line(r):
    """Two-line fit summary shown under the plot title (matplotlib mathtext)."""
    a_b, q10 = r["popt"]
    a_e, q10_e = r["perr"]
    unit = _mt(r.get("unit", ""))
    return (f"Est. activity at {r['base_t']:.2f}{DEG}C = {fmt_a(a_b)} \u00B1 {fmt_a(a_e)}"
            f"{' ' + unit if unit else ''}\n"
            f"Q$_{{10}}$ = {q10:.3f} \u00B1 {q10_e:.3f}    R$^2$ = {r['r2']:.4f}")


def report_text(r, show_tab_delimited=False):
    x, y, xa, ya = r["x"], r["y"], r["xa"], r["ya"]
    a_b, q10 = r["popt"]
    a_e, q10_e = r["perr"]
    arr, unit = r["arr"], r.get("unit", "")
    us = f" {unit}" if unit else ""
    hdr = f"Activity ({unit})" if unit else "Activity"

    out = ["INPUT DATA", table_text(x, y), "", stats_text(x, "temperature", DEG + "C", lambda v: f"{v:.3f}"),
           stats_text(y, "activity", unit, fmt_a)]
    if r["use_range"]:
        out += [f"DATA USED FOR Q10 CALCULATION ({r['low']:g} - {r['high']:g} {DEG}C)", table_text(xa, ya), "",
                stats_text(xa, "temperature", DEG + "C", lambda v: f"{v:.3f}"),
                stats_text(ya, "activity", unit, fmt_a)]
    if show_tab_delimited:
        out += ["DATASET", f"Temp ({DEG}C)\t{hdr}"]
        out += [f"{a:.3f}\t{fmt_a(b)}" for a, b in zip(x, y)]
        out += ["", "FITTED CURVE", f"Temp ({DEG}C)\t{hdr}"]
        out += [f"{a:.3f}\t{fmt_a(b)}" for a, b in zip(r["fcx"], r["fcy"])]
        out += ["", "ARRHENIUS REGRESSION LINE", f"Temp ({DEG}C)\t1000/T (1/K)\tln(Activity)"]
        out += [f"{t:.3f}\t{1000 / (t + 273.15):.5f}\t{arrhenius_lna(arr, t):.5f}" for t in r["fcx"]]
        out.append("")
    out += [
        "RESULTS",
        "Model: A(T) = A_b * Q10^((T - T_b)/10)",
        f"Estimated activity at {r['base_t']:.2f} {DEG}C = {fmt_a(a_b)} \u00B1 {fmt_a(a_e)}{us}",
        f"Q10 (temperature coefficient) = {q10:.3f} \u00B1 {q10_e:.3f}",
        f"R\u00B2 = {r['r2']:.6f}",
        "",
        "ARRHENIUS REGRESSION (ln(Activity) vs. 1000/T, data used for Q10 calculation)",
        f"Activation energy Ea = {arr['Ea'] / 1000:.3f} \u00B1 {arr['se_Ea'] / 1000:.3f} kJ/mol",
        f"Slope d ln(A) / d(1000/T) = {arr['slope']:.4f} \u00B1 {arr['se_slope']:.4f} (10\u00B3 K)",
        f"Intercept = {arr['intercept']:.4f} \u00B1 {arr['se_int']:.4f}",
        f"R\u00B2 = {arr['r2']:.6f}",
    ]
    return "\n".join(out)


# ======================================================================
# Plotting (each function draws onto a matplotlib Axes)
#   activity     : x = T (°C),               y = A
#   ln(activity) : x = T (°C),               y = ln A
#   arrhenius    : x = 1000/(T + 273.15),    y = ln A
# Every drawer takes (ax, r, show_arr, axes).
# ======================================================================
def _mt(s):
    """Escape text for matplotlib mathtext."""
    return s.replace("$", r"\$")


def _inv(t):
    return 1000.0 / (np.asarray(t, float) + 273.15)


def _curve_x(r, axes):
    """Temperatures for the fitted curve: the default span, widened to a manual temperature axis."""
    a, b = r["span"]
    if axes and axes.get("t"):
        a, b = min(a, axes["t"][0]), max(b, axes["t"][1])
    return np.linspace(a, b, 240)


def _scatter(ax, r, xs_all, ys_all, xs_used, ys_used):
    """Data points: red = used in the fit, light green = outside the fit range."""
    if not r["use_range"]:
        ax.scatter(xs_used, ys_used, s=35, color=C_USED, label="Data", zorder=Z_POINT)
    else:
        ax.scatter(xs_all, ys_all, s=30, color=C_EXCL, label="Data", zorder=Z_POINT)
        ax.scatter(xs_used, ys_used, s=35, color=C_USED,
                   label="Data points used for Q$_{10}$ calculation", zorder=Z_POINT + 0.1)


def _fix_limits(ax, pairs):
    """Axis limits covering all (x, y) pairs plus a 5 % margin."""
    xs, ys = [], []
    for px, py in pairs:
        px, py = np.ravel(px).astype(float), np.ravel(py).astype(float)
        ok = np.isfinite(px) & np.isfinite(py)
        xs.append(px[ok]); ys.append(py[ok])
    xs, ys = np.concatenate(xs), np.concatenate(ys)
    if not len(xs):
        return
    x0, x1, y0, y1 = xs.min(), xs.max(), ys.min(), ys.max()
    dx, dy = (x1 - x0) or 1.0, (y1 - y0) or 1.0
    ax.set_xlim(x0 - 0.05 * dx, x1 + 0.05 * dx)
    ax.set_ylim(y0 - 0.05 * dy, y1 + 0.05 * dy)


def auto_activity_top(r):
    """Upper limit of the automatic activity axis (same rule as _fix_limits)."""
    ys = np.concatenate([r["y"], r["ref"]["fcy"]])
    ys = ys[np.isfinite(ys)]
    y0, y1 = ys.min(), ys.max()
    return y1 + 0.05 * ((y1 - y0) or 1.0)


def _fit_label(r):
    return f"Fit Q$_{{10}}$ = {r['popt'][1]:5.3f}"


def draw_activity(ax, r, show_arr=False, axes=None):
    cx = _curve_x(r, axes)
    _scatter(ax, r, r["x"], r["y"], r["xa"], r["ya"])
    ax.plot(cx, r["model"](cx, *r["popt"]), "--", color=C_FIT, zorder=Z_LINE, label=_fit_label(r))
    _fix_limits(ax, [(r["x"], r["y"]), (r["fcx"], r["ref"]["fcy"])])
    unit = _mt(r.get("unit", ""))
    ax.set_xlabel(f"Temperature ({DEG}C)")
    ax.set_ylabel(f"Activity ({unit})" if unit else "Activity")
    ax.set_title("Temperature vs. Activity")
    ax.legend()


def draw_lnact(ax, r, show_arr=False, axes=None):
    cx = _curve_x(r, axes)
    _scatter(ax, r, r["x"], np.log(r["y"]), r["xa"], np.log(r["ya"]))
    ax.plot(cx, np.log(r["model"](cx, *r["popt"])), "--", color=C_FIT, zorder=Z_LINE, label=_fit_label(r))
    _fix_limits(ax, [(r["x"], np.log(r["y"])), (r["fcx"], np.log(r["ref"]["fcy"]))])
    ax.set_xlabel(f"Temperature ({DEG}C)")
    ax.set_ylabel("ln(Activity)")
    ax.set_title("Temperature vs. ln(Activity)")
    ax.legend()


def draw_arrhenius(ax, r, show_arr=False, axes=None):
    """Arrhenius plot; with show_arr also the Arrhenius regression line (solid grey)."""
    cx = _curve_x(r, axes)
    _scatter(ax, r, _inv(r["x"]), np.log(r["y"]), _inv(r["xa"]), np.log(r["ya"]))
    ax.plot(_inv(cx), np.log(r["model"](cx, *r["popt"])), "--", color=C_FIT, zorder=Z_LINE, label=_fit_label(r))
    if show_arr:
        ax.plot(_inv(cx), arrhenius_lna(r["arr"], cx), "-", color=C_ARR, lw=1.5, zorder=Z_LINE,
                label=f"Arrhenius regression E$_a$ = {r['arr']['Ea'] / 1000:.2f} kJ/mol")
    pairs = [(_inv(r["x"]), np.log(r["y"])), (_inv(r["fcx"]), np.log(r["ref"]["fcy"]))]
    if show_arr:
        pairs.append((_inv(r["fcx"]), arrhenius_lna(r["ref"]["arr"], r["fcx"])))
    _fix_limits(ax, pairs)
    ax.set_xlabel(r"1000/Temperature (K$^{-1}$)")
    ax.set_ylabel("ln(Activity)")
    ax.set_title("Arrhenius Plot: 1000/T vs. ln(Activity)")
    ax.legend()


DRAWERS = [("Activity", draw_activity), ("ln(Activity)", draw_lnact), ("Arrhenius", draw_arrhenius)]

# x-axis value <-> temperature (°C); used for the range band and drag-to-select.
TO_X = {"Activity": lambda t: t, "ln(Activity)": lambda t: t,
        "Arrhenius": lambda t: 1000.0 / (t + 273.15)}
FROM_X = {"Activity": lambda v: v, "ln(Activity)": lambda v: v,
          "Arrhenius": lambda v: 1000.0 / v - 273.15}


def apply_axes(ax, name, axes):
    """Replace the automatic limits of one plot by the manual ranges that are set.
    axes = {"t": (Tmin, Tmax) or None, "a": (Amin, Amax) or None}."""
    if not axes:
        return
    if axes.get("t"):
        a, b = axes["t"]
        if name == "Arrhenius":
            ax.set_xlim(1000.0 / (b + 273.15), 1000.0 / (a + 273.15))
        else:
            ax.set_xlim(a, b)
    if axes.get("a"):
        lo, hi = axes["a"]
        if name == "Activity":
            ax.set_ylim(lo, hi)
        else:
            top = math.log(hi)
            bottom = math.log(lo) if lo > 0 else ax.get_ylim()[0]
            if top > bottom:
                ax.set_ylim(bottom, top)


def apply_header(fig):
    """(Re)place the title and the fit summary at the top of the figure."""
    hdr = getattr(fig, "_q10_header", None)
    if hdr is None:
        return
    title, summary, tsize, ssize = hdr
    for t in getattr(fig, "_q10_texts", []):
        try:
            t.remove()
        except Exception:
            pass
    fig._q10_texts = []
    w, h = fig.get_size_inches()
    y_in = 0.10
    if title:
        fig._q10_texts.append(fig.text(0.5, 1 - y_in / h, title, ha="center", va="top", fontsize=tsize))
        y_in += tsize * 1.25 / 72 + 0.25
    fig._q10_texts.append(fig.text(0.5, 1 - y_in / h, summary, ha="center", va="top", fontsize=ssize))
    y_in += 2 * ssize * 1.3 / 72 + 0.12
    try:
        fig.tight_layout(rect=[0, 0, 1, max(0.3, 1 - y_in / h)])
    except Exception:
        pass


def render_single(fig, name, drawer, r, title, show_arr=False, axes=None):
    """Draw one plot (for a tab) with the title and fit summary; return its Axes."""
    fig.clear()
    ax = fig.add_subplot(111)
    drawer(ax, r, show_arr, axes)
    apply_axes(ax, name, axes)
    if r["use_range"]:
        ax.axvspan(TO_X[name](r["low"]), TO_X[name](r["high"]), color=C_BAND, alpha=0.09, lw=0, zorder=0)
    fig._q10_header = (title, summary_line(r), TITLE_SIZE, 9)
    apply_header(fig)
    return ax


def render_combined(fig, r, title, show_arr=False, axes=None):
    """Figure for saving: title, summary, then the three plots stacked vertically."""
    fig.clear()
    for i, (name, drawer) in enumerate(DRAWERS, start=1):
        ax = fig.add_subplot(3, 1, i)
        drawer(ax, r, show_arr, axes)
        apply_axes(ax, name, axes)
        ax.set_box_aspect(3 / 4)
    fig._q10_header = (title, summary_line(r), TITLE_SIZE_COMBINED, 10)
    apply_header(fig)


def make_span_selector(ax, onselect):
    """Horizontal SpanSelector, compatible with matplotlib 3.3-3.4 (rectprops) and 3.5+ (props)."""
    style = dict(facecolor=C_BAND, alpha=0.16)
    try:
        return SpanSelector(ax, onselect, "horizontal", useblit=True, props=style)
    except TypeError:
        return SpanSelector(ax, onselect, "horizontal", useblit=True, rectprops=style)


# ======================================================================
# Help window
# ======================================================================
HELP_TEXT = [
    ("h2", "What Q10 means"),
    ("p", "Q10 is the factor by which a rate changes for a 10 °C rise in temperature. From rates "
          "R1 and R2 measured at temperatures T1 and T2:"),
    ("eq", "Q10 = (R2 / R1) ^ (10 / (T2 − T1))"),
    ("p", "Typical enzyme-catalysed reactions have Q10 of about 2 to 3. Q10 > 1 means the activity "
          "rises with temperature, Q10 = 1 means no temperature dependence, and Q10 < 1 means the "
          "activity falls as temperature rises."),
    ("h2", "The model that is fitted"),
    ("eq", "A(T) = A_b · Q10 ^ ((T − T_b) / 10)        (1)"),
    ("eq", "ln A(T) = ln A_b + ((T − T_b) / 10) · ln Q10        (2)"),
    ("p", "A_b is the activity at the base temperature T_b. The base temperature only decides at which "
          "temperature A_b is reported: Q10, its standard error, R² and the fitted curve are the same "
          "whatever base you choose."),
    ("h2", "The three plots"),
    ("p", "Activity: activity against temperature, with the fitted curve from model (1).\n"
          "ln(Activity): the logarithm of the activity against temperature; model (2) is a straight line here.\n"
          "Arrhenius: ln(activity) against 1000/T, where T is in kelvin. The optional solid grey line is the "
          "straight-line regression of ln(activity) on 1000/T, giving the activation energy "
          "Ea = −slope · 1000 · R."),
    ("p", "Red points are used in the fit, light-green points are outside the chosen temperature range, "
          "and the dashed line is the fitted Q10 curve."),
    ("h2", "How the fit is calculated"),
    ("p", "Non-linear least squares on the activity itself (not its logarithm), using SciPy's curve_fit. "
          "Starting values come from a straight-line regression of ln(activity) on temperature. The ± values "
          "are the square roots of the diagonal of s² · (JᵀJ)⁻¹ with s² = SSE / (n − 2). "
          "R² = 1 − SSE / SST."),
    ("h2", "Axis scales"),
    ("p", "The axis scales are set automatically from all data. To set them yourself, tick Temperature "
          "and/or Activity under Axis scales and enter Min and Max; press Enter, click the arrows or move "
          "to another field to apply. Untick to return to automatic scaling."),
    ("p", "Temperature: 0–100 °C in steps of 0.5 °C. The Arrhenius x axis is converted to 1000/T.\n"
          "Activity: 0 to 100 × the largest activity in the data. The step of the arrows is set from the "
          "data: (data maximum / 50) rounded to 1, 2 or 5 × 10ⁿ. The ln(Activity) and Arrhenius y axes are "
          "converted to ln; with Min = 0 their lower limit stays automatic."),
    ("h2", "How to use it"),
    ("p", "1. Click Clear, then enter one temperature (°C) and activity pair per line, or use Load CSV... "
          "(temperature in column 1, activity in column 2). Enter the activity unit if you like.\n"
          "2. To fit only part of the data, tick 'Use only data within a temperature range' and set Low and "
          "High, or drag across any plot after calculating. Double-click a plot or click Use all data to go "
          "back to all data.\n"
          "3. Choose the base temperature and click Calculate (Ctrl+Enter).\n"
          "4. Save figure... saves the three plots stacked in one PNG, PDF or SVG image; Save results... saves "
          "the results text."),
    ("h2", "Privacy"),
    ("p", "All calculations run on your computer. The program makes no network connections and stores "
          "nothing; files are written only when you click a Save button."),
]


class HelpWindow(tk.Toplevel):
    def __init__(self, master):
        super().__init__(master)
        self.title("How this calculator works")
        self.geometry(f"760x{min(720, master.winfo_screenheight() - 120)}")
        self.transient(master)
        base = tkfont.nametofont("TkDefaultFont")
        fam, size = base.actual("family"), max(10, int(base.actual("size")))
        bar = ttk.Frame(self, padding=(12, 8))
        bar.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(bar, text="How this calculator works", font=(fam, size + 4, "bold")).pack(side=tk.LEFT)
        ttk.Button(bar, text="Close", command=self.destroy).pack(side=tk.RIGHT)
        ttk.Separator(self).pack(side=tk.TOP, fill=tk.X)
        box = ttk.Frame(self)
        box.pack(fill=tk.BOTH, expand=True)
        txt = tk.Text(box, wrap="word", padx=18, pady=10, borderwidth=0, highlightthickness=0,
                      font=(fam, size + 1), spacing1=2, spacing3=2)
        sb = ttk.Scrollbar(box, orient="vertical", command=txt.yview)
        txt.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        txt.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        txt.tag_configure("h2", font=(fam, size + 4, "bold"), spacing1=14, spacing3=6)
        txt.tag_configure("p", spacing3=8)
        txt.tag_configure("eq", lmargin1=24, lmargin2=24, background="#F1F2F4", spacing1=4, spacing3=8)
        for kind, text in HELP_TEXT:
            txt.insert(tk.END, text + "\n", (kind,))
        txt.configure(state="disabled")
        self.bind("<Escape>", lambda e: self.destroy())
        self.focus_set()


# ======================================================================
# GUI
# ======================================================================
class Q10App(tk.Tk):
    RANGE_DELAY_MS = 400

    def __init__(self):
        super().__init__()
        self.title("Q10 calculator - temperature vs. activity")
        self.geometry(f"1280x{min(880, self.winfo_screenheight() - 80)}")
        self.minsize(1000, 600)

        self.result = None
        self.help_win = None
        self.plots = []          # one dict per plot tab
        self.a_scale = None      # activity-axis scale (activity_axis_scale())

        self.var_title = tk.StringVar()
        self.var_unit = tk.StringVar()
        self.var_range = tk.BooleanVar(value=False)
        self.var_low = tk.StringVar(value="15")
        self.var_high = tk.StringVar(value="35")
        self.var_base = tk.StringVar(value=BASE_CHOICES[3])   # 25°C
        self.var_custom = tk.StringVar(value="25")
        self.var_tab = tk.BooleanVar(value=False)
        self.var_arr = tk.BooleanVar(value=True)
        self.var_axt, self.var_axa = tk.BooleanVar(value=False), tk.BooleanVar(value=False)
        self.var_tmin, self.var_tmax = tk.StringVar(), tk.StringVar()
        self.var_amin, self.var_amax = tk.StringVar(), tk.StringVar()

        self._build_controls()
        self._build_output()
        self._toggle_states()

        self._range_job = None
        self._suppress_range = False
        for var in (self.var_low, self.var_high):
            var.trace_add("write", lambda *a: self._range_changed())
        self.var_unit.trace_add("write", lambda *a: self._unit_changed())
        self.bind("<Control-Return>", lambda e: self.run())

        self.load_example()
        self.run()
        self.status.set("Showing the demo calculation. Click Clear to start with your own data.")
        self.notice.pack(side=tk.TOP, fill=tk.X, pady=(0, 6), before=self.nb)

    # ---- layout ------------------------------------------------------
    def _build_controls(self):
        left = ttk.Frame(self, padding=8)
        left.pack(side=tk.LEFT, fill=tk.Y)

        head = ttk.Frame(left)
        head.pack(side=tk.TOP, fill=tk.X, pady=(0, 6))
        ttk.Label(head, text="Q10 calculator: temperature vs. activity",
                  font=("TkDefaultFont", 12, "bold")).pack(side=tk.LEFT)
        ttk.Button(head, text="How it works", command=self.show_help).pack(side=tk.RIGHT, padx=(8, 0))

        f = ttk.LabelFrame(left, text="Plot title", padding=6)
        f.pack(side=tk.TOP, fill=tk.X, pady=(0, 6))
        ttk.Entry(f, textvariable=self.var_title, width=46).pack(fill=tk.X)
        self.var_title.trace_add("write", lambda *a: self._redraw_plots())

        # Everything below the data box is packed to the BOTTOM first, so on a short
        # screen only the data box shrinks.
        self.status = tk.StringVar(value="Ready.")
        self.status_label = ttk.Label(left, textvariable=self.status, foreground="#555", wraplength=350)
        self.status_label.pack(side=tk.BOTTOM, anchor="w", pady=(8, 0))

        row = ttk.Frame(left)
        row.pack(side=tk.BOTTOM, fill=tk.X, pady=(6, 0))
        ttk.Button(row, text="Save figure...", command=self.save_figure).pack(side=tk.LEFT, expand=True, fill=tk.X)
        ttk.Button(row, text="Save results...", command=self.save_results).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=(4, 0))
        ttk.Button(left, text="Calculate (Ctrl+Enter)", command=self.run).pack(side=tk.BOTTOM, fill=tk.X, ipady=4)
        ttk.Checkbutton(left, text="Also show tab-delimited data in Results", variable=self.var_tab,
                        command=self._refresh_results).pack(side=tk.BOTTOM, anchor="w", pady=(0, 6))
        ttk.Checkbutton(left, text="Show Arrhenius regression line on the Arrhenius plot",
                        variable=self.var_arr, command=self._redraw_plots).pack(side=tk.BOTTOM, anchor="w", pady=(0, 2))

        # Axis scales: only redraw the plots; nothing is recalculated.
        f = ttk.LabelFrame(left, text="Axis scales", padding=6)
        f.pack(side=tk.BOTTOM, fill=tk.X, pady=(0, 6))
        ttk.Label(f, text="Tick to set; press Enter or move to another field to apply.",
                  foreground="#555").grid(row=0, column=0, columnspan=6, sticky="w", pady=(0, 4))
        ttk.Checkbutton(f, text="Temperature", variable=self.var_axt,
                        command=lambda: self._axis_toggled("t")).grid(row=1, column=0, sticky="w")
        ttk.Label(f, text="Min").grid(row=1, column=1, padx=(6, 2))
        self.sp_tmin = ttk.Spinbox(f, textvariable=self.var_tmin, width=8, from_=T_AXIS_RANGE[0],
                                   to=T_AXIS_RANGE[1], increment=T_AXIS_STEP, format="%.1f",
                                   command=self._apply_axes)
        self.sp_tmin.grid(row=1, column=2)
        ttk.Label(f, text="Max").grid(row=1, column=3, padx=(6, 2))
        self.sp_tmax = ttk.Spinbox(f, textvariable=self.var_tmax, width=8, from_=T_AXIS_RANGE[0],
                                   to=T_AXIS_RANGE[1], increment=T_AXIS_STEP, format="%.1f",
                                   command=self._apply_axes)
        self.sp_tmax.grid(row=1, column=4)
        ttk.Label(f, text=f"{DEG}C").grid(row=1, column=5, padx=(4, 0))
        ttk.Label(f, text=f"Range 0–100 {DEG}C, step {T_AXIS_STEP} {DEG}C", foreground="#555"
                  ).grid(row=2, column=1, columnspan=5, sticky="w", pady=(0, 6))
        ttk.Checkbutton(f, text="Activity", variable=self.var_axa,
                        command=lambda: self._axis_toggled("a")).grid(row=3, column=0, sticky="w")
        ttk.Label(f, text="Min").grid(row=3, column=1, padx=(6, 2))
        self.sp_amin = ttk.Spinbox(f, textvariable=self.var_amin, width=8, from_=0, to=1,
                                   increment=1, command=self._apply_axes)
        self.sp_amin.grid(row=3, column=2)
        ttk.Label(f, text="Max").grid(row=3, column=3, padx=(6, 2))
        self.sp_amax = ttk.Spinbox(f, textvariable=self.var_amax, width=8, from_=0, to=1,
                                   increment=1, command=self._apply_axes)
        self.sp_amax.grid(row=3, column=4)
        self.a_hint = tk.StringVar(value="Range and step are set from the data after calculating.")
        ttk.Label(f, textvariable=self.a_hint, foreground="#555", wraplength=300
                  ).grid(row=4, column=1, columnspan=5, sticky="w")
        for sp in (self.sp_tmin, self.sp_tmax, self.sp_amin, self.sp_amax):
            sp.bind("<Return>", lambda e: self._apply_axes())
            sp.bind("<KP_Enter>", lambda e: self._apply_axes())
            sp.bind("<FocusOut>", lambda e: self._apply_axes())

        # Base temperature
        f = ttk.LabelFrame(left, text="Base temperature", padding=6)
        f.pack(side=tk.BOTTOM, fill=tk.X, pady=(0, 6))
        cb = ttk.Combobox(f, textvariable=self.var_base, values=BASE_CHOICES, state="readonly", width=30)
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: (self._toggle_states(), self._rerun()))
        r = ttk.Frame(f)
        r.pack(fill=tk.X, pady=(4, 0))
        ttk.Label(r, text="Custom value").pack(side=tk.LEFT)
        self.ent_custom = ttk.Entry(r, textvariable=self.var_custom, width=8)
        self.ent_custom.pack(side=tk.LEFT, padx=4)
        ttk.Label(r, text=f"{DEG}C").pack(side=tk.LEFT)
        self.ent_custom.bind("<Return>", lambda e: self._rerun())
        self.ent_custom.bind("<FocusOut>", lambda e: self._rerun())

        # Temperature range for the fit
        f = ttk.LabelFrame(left, text="Temperature range for the fit", padding=6)
        f.pack(side=tk.BOTTOM, fill=tk.X, pady=(0, 6))
        ttk.Checkbutton(f, text="Use only data within a temperature range", variable=self.var_range,
                        command=lambda: (self._toggle_states(), self._rerun())).pack(anchor="w")
        r = ttk.Frame(f)
        r.pack(fill=tk.X, pady=(4, 0))
        ttk.Label(r, text="Low").pack(side=tk.LEFT)
        self.sp_low = ttk.Spinbox(r, textvariable=self.var_low, width=7, from_=-100, to=200, increment=0.5)
        self.sp_low.pack(side=tk.LEFT, padx=4)
        ttk.Label(r, text="High").pack(side=tk.LEFT)
        self.sp_high = ttk.Spinbox(r, textvariable=self.var_high, width=7, from_=-100, to=200, increment=0.5)
        self.sp_high.pack(side=tk.LEFT, padx=4)
        ttk.Label(r, text=f"{DEG}C").pack(side=tk.LEFT)
        ttk.Button(r, text="Use all data", command=self.use_all_data).pack(side=tk.RIGHT)

        # Data box (fills the remaining height)
        f = ttk.LabelFrame(left, text=f"Data: temperature ({DEG}C), activity — one pair per line", padding=6)
        f.pack(side=tk.TOP, fill=tk.BOTH, expand=True, pady=(0, 6))
        r = ttk.Frame(f)
        r.pack(side=tk.BOTTOM, fill=tk.X, pady=(4, 0))
        ttk.Label(r, text="Activity unit").pack(side=tk.LEFT)
        ttk.Entry(r, textvariable=self.var_unit, width=16).pack(side=tk.LEFT, padx=4, fill=tk.X, expand=True)
        r = ttk.Frame(f)
        r.pack(side=tk.BOTTOM, fill=tk.X, pady=(4, 0))
        ttk.Button(r, text="Load CSV...", command=self.load_csv).pack(side=tk.LEFT)
        ttk.Button(r, text="Load example", command=self.load_example).pack(side=tk.LEFT, padx=4)
        ttk.Button(r, text="Clear", command=self.clear).pack(side=tk.LEFT)
        box = ttk.Frame(f)
        box.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.txt_data = tk.Text(box, width=34, height=DATA_ROWS, font=("TkFixedFont", 10), undo=True)
        sb = ttk.Scrollbar(box, orient="vertical", command=self.txt_data.yview)
        self.txt_data.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.txt_data.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.txt_data.bind("<<Modified>>", self._data_modified)

    def _build_output(self):
        right = ttk.Frame(self, padding=(0, 8, 8, 8))
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.notice = ttk.Label(right, text=NOTICE_DEMO, wraplength=820, padding=8,
                                background="#EEF0F3", foreground="#1B1F24")
        self.nb = ttk.Notebook(right)
        self.nb.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        tab = ttk.Frame(self.nb)
        self.nb.add(tab, text="Results")
        self.txt_res = tk.Text(tab, wrap="none", font=("TkFixedFont", 10))
        sy = ttk.Scrollbar(tab, orient="vertical", command=self.txt_res.yview)
        sx = ttk.Scrollbar(tab, orient="horizontal", command=self.txt_res.xview)
        self.txt_res.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
        sy.pack(side=tk.RIGHT, fill=tk.Y)
        sx.pack(side=tk.BOTTOM, fill=tk.X)
        self.txt_res.pack(fill=tk.BOTH, expand=True)

        for name, drawer in DRAWERS:
            tab = ttk.Frame(self.nb)
            self.nb.add(tab, text=name)
            fig = Figure(figsize=(7.5, 6), dpi=100)
            canvas = FigureCanvasTkAgg(fig, master=tab)
            toolbar = NavigationToolbar2Tk(canvas, tab, pack_toolbar=False)
            toolbar.update()
            toolbar.pack(side=tk.BOTTOM, fill=tk.X)
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            p = dict(name=name, drawer=drawer, fig=fig, canvas=canvas, toolbar=toolbar, selector=None)
            canvas.mpl_connect("resize_event", lambda e, fg=fig, cv=canvas: (apply_header(fg), cv.draw_idle()))
            canvas.mpl_connect("button_press_event", lambda e: self.use_all_data() if e.dblclick else None)
            self.plots.append(p)
        self.nb.select(1)

    # ---- state helpers -------------------------------------------------
    def _set_state(self, widget, enabled):
        widget.state(["!disabled"] if enabled else ["disabled"])

    def _toggle_states(self):
        self._set_state(self.sp_low, self.var_range.get())
        self._set_state(self.sp_high, self.var_range.get())
        self._set_state(self.ent_custom, self.var_base.get() == CUSTOM_BASE)
        have = self.result is not None
        for sp in (self.sp_tmin, self.sp_tmax):
            self._set_state(sp, have and self.var_axt.get())
        for sp in (self.sp_amin, self.sp_amax):
            self._set_state(sp, have and self.var_axa.get() and self.a_scale is not None)

    def _error(self, msg, popup=False):
        self.status.set(msg)
        self.status_label.configure(foreground="#9A2A18")
        if popup:
            messagebox.showerror("Q10 calculator", msg, parent=self)

    def _ok(self, msg):
        self.status.set(msg)
        self.status_label.configure(foreground="#555")

    def _hide_notice(self):
        self.notice.pack_forget()

    def _data_modified(self, _event=None):
        if self.txt_data.edit_modified():
            self.txt_data.edit_modified(False)
            if getattr(self, "_loading", False):
                return
            self._hide_notice()

    # ---- data ----------------------------------------------------------
    def _set_data(self, text):
        self._loading = True
        self.txt_data.delete("1.0", tk.END)
        self.txt_data.insert("1.0", text)
        self.update_idletasks()
        self.txt_data.edit_modified(False)
        self._loading = False

    def load_example(self):
        self._set_data(EXAMPLE_DATA)
        self.var_title.set("Example")
        self.var_unit.set(EXAMPLE_UNIT)
        self._ok("Example data loaded.")

    def clear(self):
        self._hide_notice()
        self._set_data("")
        self.var_title.set("")
        self.var_unit.set("")
        self.result = None
        self.a_scale = None
        self.var_axt.set(False); self.var_axa.set(False)
        for v in (self.var_tmin, self.var_tmax, self.var_amin, self.var_amax):
            v.set("")
        self.a_hint.set("Range and step are set from the data after calculating.")
        self.txt_res.delete("1.0", tk.END)
        for p in self.plots:
            p["fig"].clear()
            p["selector"] = None
            p["canvas"].draw_idle()
        self._toggle_states()
        self._ok("Cleared.")

    def load_csv(self):
        path = filedialog.askopenfilename(parent=self, title="Load CSV",
                                          filetypes=[("CSV / text", "*.csv *.tsv *.txt"), ("All files", "*.*")])
        if not path:
            return
        rows = []
        try:
            with open(path, newline="", encoding="utf-8-sig") as fh:
                sample = fh.read(4096)
                fh.seek(0)
                try:
                    dialect = csv.Sniffer().sniff(sample, delimiters=",;\t ")
                except csv.Error:
                    dialect = csv.excel
                for cells in csv.reader(fh, dialect):
                    cells = [c.strip() for c in cells if c.strip()]
                    if len(cells) >= 2:
                        try:
                            rows.append((float(cells[0]), float(cells[1])))
                        except ValueError:
                            continue
        except OSError as e:
            self._error(f"Load failed: {e}", popup=True)
            return
        if not rows:
            self._error("Load failed: no numeric rows found (need two numeric columns).", popup=True)
            return
        self._hide_notice()
        self._set_data("\n".join(f"{a:g}, {b:g}" for a, b in rows))
        self._ok(f"Loaded {len(rows)} rows. Click Calculate.")

    # ---- calculation ---------------------------------------------------
    def run(self, quiet=False):
        try:
            temps, acts = parse_pairs(self.txt_data.get("1.0", tk.END))
            use_range = self.var_range.get()
            low = high = 0.0
            if use_range:
                low, high = float(self.var_low.get()), float(self.var_high.get())
            custom = 0.0
            if self.var_base.get() == CUSTOM_BASE:
                try:
                    custom = float(self.var_custom.get())
                except ValueError:
                    raise ValueError("Enter a custom base temperature.") from None
            r = analyze(temps, acts, use_range, low, high, self.var_base.get(), custom)
        except (ValueError, RuntimeError) as e:
            self._error(f"Cannot calculate: {e}", popup=not quiet)
            return
        r["unit"] = self.var_unit.get().strip()
        self.result = r
        self._configure_a_scale(float(np.max(r["y"])))
        self._toggle_states()
        self._refresh_results()
        self._redraw_plots()
        a_b, q10 = r["popt"]
        self._ok(f"Q10 = {q10:.3f}, activity at {r['base_t']:.2f}{DEG}C = {fmt_a(a_b)}, R\u00B2 = {r['r2']:.4f}")

    def _rerun(self):
        if self.result is not None:
            self.run(quiet=True)

    def _range_changed(self):
        if self._suppress_range or self.result is None or not self.var_range.get():
            return
        if self._range_job:
            self.after_cancel(self._range_job)
        self._range_job = self.after(self.RANGE_DELAY_MS, self._rerun)

    def use_all_data(self):
        if self.var_range.get():
            self.var_range.set(False)
            self._toggle_states()
            self._rerun()

    def _span_selected(self, name, toolbar, xmin, xmax):
        if toolbar.mode or self.result is None or abs(xmax - xmin) < 1e-9:
            return
        t1, t2 = FROM_X[name](xmin), FROM_X[name](xmax)
        lo, hi = min(t1, t2), max(t1, t2)
        self._suppress_range = True
        self.var_low.set(f"{lo:.1f}")
        self.var_high.set(f"{hi:.1f}")
        self._suppress_range = False
        self.var_range.set(True)
        self._toggle_states()
        self.run(quiet=True)

    def _unit_changed(self):
        if self.result is not None:
            self.result["unit"] = self.var_unit.get().strip()
            self._refresh_results()
            self._redraw_plots()

    # ---- axis scales ---------------------------------------------------
    def _configure_a_scale(self, max_a):
        s = activity_axis_scale(max_a)
        self.a_scale = s
        fmt = f"%.{s['dec']}f"
        for sp in (self.sp_amin, self.sp_amax):
            sp.configure(from_=0, to=s["limit"], increment=s["step"], format=fmt)
        self.a_hint.set(f"Range 0–{s['limit']:.{s['dec']}f} (100 \u00D7 data max {max_a:.6g}), "
                        f"step {s['step']:.{s['dec']}f} (set from the data)")
        if self.var_axa.get():
            lo, hi = self._float(self.var_amin), self._float(self.var_amax)
            lo = min(max(lo if lo is not None else 0.0, 0.0), s["limit"])
            hi = min(max(hi if hi is not None else 0.0, 0.0), s["limit"])
            if not hi > lo:
                lo, hi = self._default_a()
            self.var_amin.set(f"{lo:.{s['dec']}f}")
            self.var_amax.set(f"{hi:.{s['dec']}f}")

    def _default_a(self):
        s = self.a_scale
        top = math.ceil(auto_activity_top(self.result) / s["step"] - 1e-9) * s["step"]
        return 0.0, min(round(top, s["dec"]), s["limit"])

    def _default_t(self):
        x = self.result["x"]
        fcx = self.result["fcx"]
        x0, x1 = min(x.min(), fcx.min()), max(x.max(), fcx.max())
        d = (x1 - x0) or 1.0
        lo = math.floor((x0 - 0.05 * d) / T_AXIS_STEP) * T_AXIS_STEP
        hi = math.ceil((x1 + 0.05 * d) / T_AXIS_STEP) * T_AXIS_STEP
        return max(lo, T_AXIS_RANGE[0]), min(hi, T_AXIS_RANGE[1])

    def _axis_toggled(self, which):
        if self.result is not None:
            if which == "t" and self.var_axt.get() and not (self.var_tmin.get() and self.var_tmax.get()):
                lo, hi = self._default_t()
                self.var_tmin.set(f"{lo:.1f}"); self.var_tmax.set(f"{hi:.1f}")
            if which == "a" and self.var_axa.get() and not (self.var_amin.get() and self.var_amax.get()):
                lo, hi = self._default_a()
                d = self.a_scale["dec"]
                self.var_amin.set(f"{lo:.{d}f}"); self.var_amax.set(f"{hi:.{d}f}")
        self._toggle_states()
        self._apply_axes()

    @staticmethod
    def _float(var):
        try:
            return float(var.get())
        except ValueError:
            return None

    def _read_axes(self):
        """Return (axes, message). Invalid ranges are left automatic and reported."""
        axes, msgs = {"t": None, "a": None}, []
        if self.var_axt.get():
            a, b = self._float(self.var_tmin), self._float(self.var_tmax)
            if a is not None and b is not None and T_AXIS_RANGE[0] <= a < b <= T_AXIS_RANGE[1]:
                axes["t"] = (a, b)
            else:
                msgs.append(f"Temperature axis: enter 0 \u2264 Min < Max \u2264 100 {DEG}C.")
        if self.var_axa.get() and self.a_scale:
            a, b = self._float(self.var_amin), self._float(self.var_amax)
            lim = self.a_scale["limit"]
            if a is not None and b is not None and 0 <= a < b <= lim + 1e-12:
                axes["a"] = (a, b)
            else:
                msgs.append(f"Activity axis: enter 0 \u2264 Min < Max \u2264 {lim:.{self.a_scale['dec']}f}.")
        return axes, " ".join(msgs)

    def _apply_axes(self):
        if self.result is None:
            return
        _, msg = self._read_axes()
        if msg:
            self._error(msg + " That axis stays automatic.")
        self._redraw_plots()

    # ---- output --------------------------------------------------------
    def _refresh_results(self):
        if self.result is None:
            return
        self.txt_res.delete("1.0", tk.END)
        self.txt_res.insert("1.0", report_text(self.result, self.var_tab.get()))

    def _redraw_plots(self):
        if self.result is None:
            return
        axes, _ = self._read_axes()
        title = self.var_title.get().strip()
        for p in self.plots:
            ax = render_single(p["fig"], p["name"], p["drawer"], self.result, title,
                               self.var_arr.get(), axes)
            p["selector"] = make_span_selector(
                ax, lambda a, b, n=p["name"], tb=p["toolbar"]: self._span_selected(n, tb, a, b))
            p["canvas"].draw_idle()

    def _file_base(self):
        t = re.sub(r"[^\w\- ]+", "", self.var_title.get()).strip()
        return re.sub(r"\s+", "_", t) or "q10_activity"

    def save_figure(self):
        if self.result is None:
            self._error("Calculate first.", popup=True)
            return
        path = filedialog.asksaveasfilename(
            parent=self, title="Save figure", defaultextension=".png",
            initialfile=self._file_base() + "_q10_plots.png",
            filetypes=[("PNG image", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
        if not path:
            return
        axes, _ = self._read_axes()
        fig = Figure(figsize=(8, 21), dpi=100)
        render_combined(fig, self.result, self.var_title.get().strip(), self.var_arr.get(), axes)
        try:
            fig.savefig(path, dpi=200, facecolor="white")
        except Exception as e:
            self._error(f"Save failed: {e}", popup=True)
            return
        self._ok(f"Saved {path}")

    def save_results(self):
        if self.result is None:
            self._error("Calculate first.", popup=True)
            return
        path = filedialog.asksaveasfilename(
            parent=self, title="Save results", defaultextension=".txt",
            initialfile=self._file_base() + "_q10_results.txt",
            filetypes=[("Text", "*.txt"), ("All files", "*.*")])
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(report_text(self.result, self.var_tab.get()) + "\n")
        except OSError as e:
            self._error(f"Save failed: {e}", popup=True)
            return
        self._ok(f"Saved {path}")

    def show_help(self):
        if self.help_win is not None and self.help_win.winfo_exists():
            self.help_win.lift()
            return
        self.help_win = HelpWindow(self)


if __name__ == "__main__":
    Q10App().mainloop()

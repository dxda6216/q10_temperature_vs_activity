# Q10 calculator: temperature vs. activity

Fit a Q<sub>10</sub> temperature coefficient to activity measured at several temperatures, and see the fit on Activity, ln(Activity) and Arrhenius plots.

**Web app:** https://dxda6216.github.io/q10_temperature_vs_activity/

The calculator comes in two versions that give the same results:

| File | Version | How to run |
|---|---|---|
| [`index.html`](index.html) | Web app (single HTML file) | Open the [web app](https://dxda6216.github.io/q10_temperature_vs_activity/), or download `index.html` and open it in a browser. No installation or server needed. |
| [`q10_temperature_vs_activity.py`](q10_temperature_vs_activity.py) | Desktop app (Python / Tkinter) | `python q10_temperature_vs_activity.py` |

## Features

- Paste temperature–activity pairs or load a CSV file.
- Non-linear least-squares fit of Q<sub>10</sub> and the activity at a base temperature of your choice, with standard errors and R².
- Three plots: Activity vs. temperature, ln(Activity) vs. temperature, and an Arrhenius plot (ln(Activity) vs. 1000/T) with an optional Arrhenius regression line and activation energy E<sub>a</sub>.
- Fit only part of the data: set a temperature range, or drag across any plot with the mouse. Double-click a plot to go back to all data.
- Manual axis scales for temperature and activity.
- Save the three plots as one image (web: PNG or SVG; desktop: PNG, PDF or SVG) and the results as a text report.
- Starts with a demo calculation on example data.

## Model

The activity A at temperature T (°C) is modelled as

```
A(T) = A_b · Q10 ^ ((T − T_b) / 10)
```

where A<sub>b</sub> is the activity at the base temperature T<sub>b</sub>. Taking logarithms gives a straight line in T:

```
ln A(T) = ln A_b + ((T − T_b) / 10) · ln Q10
```

Q<sub>10</sub> > 1 means the activity rises with temperature, Q<sub>10</sub> = 1 means no temperature dependence, and Q<sub>10</sub> < 1 means the activity falls as temperature rises. Typical enzyme-catalysed reactions have Q<sub>10</sub> of about 2 to 3.

The base temperature (0, 4, 20, 25, 30 or 37 °C, the minimum, maximum, average or median of the data, or a custom value) only decides at which temperature A<sub>b</sub> is reported. Q<sub>10</sub>, its standard error, R² and the fitted curve are the same for every choice.

## How the fit is calculated

- **Fit:** non-linear least squares on the activity itself (not its logarithm), minimising SSE = Σ (A<sub>i</sub> − Â<sub>i</sub>)². The desktop app uses SciPy's `curve_fit`; the web app uses the Levenberg–Marquardt method and reaches the same minimum.
- **Starting values:** a straight-line regression of ln(activity) on temperature.
- **Standard errors (±):** square roots of the diagonal of s² (JᵀJ)⁻¹ with s² = SSE / (n − 2), the same definition as `curve_fit` with `absolute_sigma=False`.
- **R²** = 1 − SSE / SST.
- **Arrhenius regression:** a straight line of ln(activity) against 1000/T<sub>K</sub> (T<sub>K</sub> in kelvin), fitted to the same data as Q<sub>10</sub>. The activation energy is E<sub>a</sub> = −slope × 1000 × R, with R = 8.314 J mol⁻¹ K⁻¹.

## Input data

Enter one pair per line: temperature in °C, then activity. Values can be separated by commas, spaces, tabs or semicolons. A header line and lines starting with `#` are skipped.

```
# temperature (°C), activity (U/mg)
10, 14.8
15, 23.4
20, 32.9
25, 49.1
30, 75.6
35, 111.9
40, 166.4
```

CSV files should have temperature in the first column and activity in the second. At least 3 data points are needed, and all activities must be positive. You can also enter the activity unit, which is used in the plot labels and the report.

## Temperature range for the fit

Tick **Use only data within a temperature range** and set Low and High, or drag across any plot after calculating. Points used in the fit are shown in red and points outside the range in green. Double-click a plot or click **Use all data** to return to all data. The automatic axis scales stay those of the all-data fit, so they do not jump when the range changes.

## Axis scales

The axes are scaled automatically. To set them yourself, tick **Temperature** and/or **Activity** under **Axis scales** (below the plot) and enter Min and Max. Press Enter, click the arrows or move to another field to apply; untick to return to automatic scaling.

| Axis | Range | Step of the arrows |
|---|---|---|
| Temperature | 0–100 °C | 0.5 °C |
| Activity | 0 to 100 × the largest activity in the data | Set from the data: data maximum / 50, rounded to 1, 2 or 5 × 10ⁿ |

The Arrhenius x axis is converted from the temperature range to 1000/T. The ln(Activity) and Arrhenius y axes are converted from the activity range to ln; with Min = 0 their lower limit stays automatic, because ln(0) is undefined.

## Desktop app requirements

- Python 3.9 or later
- numpy, scipy, matplotlib (Tkinter is included with most Python installations)

```
pip install numpy scipy matplotlib
python q10_temperature_vs_activity.py
```

The comment block at the top of `q10_temperature_vs_activity.py` describes the background, calculation method and usage in more detail.

## Privacy

All calculations run locally, in your browser or on your computer. The data you enter is not sent anywhere or stored. Files are written only when you use the save buttons.

## License

This project is released under the [MIT License](LICENSE).

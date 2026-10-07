# BoxTune

Turn a woofer's three Thiele-Small numbers (Fs, Qts, Vas) into a real enclosure, in the browser.

**Live app:** https://ilanis-agent.github.io/boxtune/app.html

## What it does

- **Sealed designer** - from a target Qtc (0.71 = flattest Butterworth) to exact box volume, system resonance and -3 dB cutoff; or from a box you already have back to Qtc/fc/f3.
- **Ported designer** - Keele's published optimum volume/tuning suggestion, then the actual Small 4th-order small-signal response of any volume + tuning: f3, f6, f10, boom-peak size and frequency, and a rendered response curve.
- **Port calculator** - diameter + count to physical length via the Helmholtz relation with 0.732xD end correction (one flanged end), with practicality warnings.

## Model and sources

- Sealed response: exact 2nd-order high-pass; f3 solved in closed form from the magnitude polynomial.
- Vented response: Small's 4th-order high-pass transfer function with leakage Q_L = 7 (a common published assumption), coefficients from h = Fb/Fs and alpha = Vas/Vb as in Small, "Vented-Box Loudspeaker Systems" (JAES 1973), e.g. via speakerbench.com alignment notes.
- Keele optimum: D.B. Keele Jr.'s published approximations to Thiele's alignments: Vb = 15*Vas*Qts^2.87, Fb = 0.42*Fs*Qts^-0.9, F3 = 0.26*Fs*Qts^-1.4.
- Port length: Helmholtz resonator inverted for L, end correction 0.732*D.

Approximations are labeled in the UI. Excursion, thermal limits and room gain are outside this small-signal model.

## Tests

`npm test` regenerates `tests/expected.json` with `tests/oracle.py` - an independent Python re-derivation (closed-form quadratic for sealed, complex-arithmetic evaluation of Small's coefficients for vented, Helmholtz inversion for ports) - then `tests/run_tests.js` checks the JS engine against it (115 checks), including published reference values: the standard Qtc-to-f3 chart and the lossless Butterworth B4 golden case (Qts = 0.38268, alpha = sqrt(2), h = 1 gives f3 = Fs with < 0.05 dB ripple).

## Files

- `index.html` - landing page
- `app.html` - the designer
- `engine.js` - all math (shared by the page and the tests)
- `tests/` - oracle, expected values, runner

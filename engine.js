// BoxTune engine: sealed and vented loudspeaker enclosure design from Thiele-Small parameters.
// Model: standard small-signal theory (Thiele/Small). Vented response uses Small's 4th-order
// high-pass transfer function with box leakage Q_L (default 7, a common published assumption).
// All approximations are labeled where shown.

var C_SOUND = 343.0; // speed of sound, m/s at 20 C (published standard value)

function bad(v){ return !(typeof v === 'number' && isFinite(v)); }

// ---- sealed ----
// 2nd-order high-pass with system Q Qtc and corner fc:
// |H(f)|^2 = u^4 / ((1-u^2)^2 + u^2/Qtc^2), u = f/fc  (standard 2nd-order HP magnitude)
function sealedF3Ratio(qtc) {
  // exact -3 dB point: solve |H|^2 = 1/2 -> y^2 + y(2 - 1/Qtc^2) - 1 = 0 with y = u^2
  var q2 = 1 / (qtc * qtc);
  var y = ((q2 - 2) + Math.sqrt((2 - q2) * (2 - q2) + 4)) / 2;
  return Math.sqrt(y); // f3 / fc
}
function sealedFromQtc(fs, qts, vas, qtc) {
  if ([fs, qts, vas, qtc].some(bad) || fs <= 0 || qts <= 0 || vas <= 0) return { error: 'Fs, Qts and Vas must be positive numbers.' };
  if (qtc <= qts) return { error: 'Target Qtc must be higher than the driver Qts (' + qts + '). A box only raises Q.' };
  var alpha = (qtc / qts) * (qtc / qts) - 1;
  var vb = vas / alpha;
  var fc = fs * (qtc / qts);
  return { Vb: vb, fc: fc, f3: fc * sealedF3Ratio(qtc), alpha: alpha, Qtc: qtc };
}
function sealedFromVb(fs, qts, vas, vb) {
  if ([fs, qts, vas, vb].some(bad) || fs <= 0 || qts <= 0 || vas <= 0 || vb <= 0) return { error: 'Fs, Qts, Vas and box volume must be positive numbers.' };
  var alpha = vas / vb;
  var qtc = qts * Math.sqrt(alpha + 1);
  var fc = fs * Math.sqrt(alpha + 1);
  return { Vb: vb, fc: fc, f3: fc * sealedF3Ratio(qtc), alpha: alpha, Qtc: qtc };
}

// ---- vented: Keele's published approximations to Small's optimum alignments ----
// (D.B. Keele Jr., from A.N. Thiele's alignment tables - widely republished)
function keeleOptimum(fs, qts, vas) {
  if ([fs, qts, vas].some(bad) || fs <= 0 || qts <= 0 || vas <= 0) return { error: 'Fs, Qts and Vas must be positive numbers.' };
  return {
    Vb: 15 * vas * Math.pow(qts, 2.87),
    Fb: 0.42 * fs * Math.pow(qts, -0.9),
    F3: 0.26 * fs * Math.pow(qts, -1.4)
  };
}
function keeleGivenVb(fs, qts, vas, vb) {
  if ([fs, qts, vas, vb].some(bad) || vb <= 0) return { error: 'Volumes must be positive numbers.' };
  return { F3: Math.sqrt(vas / vb) * fs, Fb: Math.pow(vas / vb, 0.32) * fs };
}

// ---- vented: Small 4th-order transfer function ----
// G(u) = u^4 / (u^4 + a1 u^3 + a2 u^2 + a3 u + 1), u = f/fb (normalized to box tuning)
// h = fb/fs, alpha = Vas/Vb, QL = box leakage Q (default 7):
//   a1 = (QL + h*Qts)/(sqrt(h)*QL*Qts)
//   a2 = (h + (alpha+1+h^2)*QL*Qts)/(h*QL*Qts)
//   a3 = (h*QL + Qts)/(sqrt(h)*QL*Qts)
function ventedCoeffs(fs, qts, vas, vb, fb, ql) {
  ql = ql || 7;
  var h = fb / fs, alpha = vas / vb, sh = Math.sqrt(h);
  return {
    a1: (ql + h * qts) / (sh * ql * qts),
    a2: (h + (alpha + 1 + h * h) * ql * qts) / (h * ql * qts),
    a3: (h * ql + qts) / (sh * ql * qts),
    h: h, alpha: alpha, ql: ql
  };
}
function ventedMagSq(u, c) {
  // |G(u)|^2 = u^8 / |D(u)|^2, D(u) = u^4 + a1(ju)^3-> re: u^4 - a2 u^2 + 1, im: -a1 u^3 + a3 u
  var re = u * u * u * u - c.a2 * u * u + 1;
  var im = -c.a1 * u * u * u + c.a3 * u;
  var u8 = Math.pow(u, 8);
  return u8 / (re * re + im * im);
}
function ventedAnalyze(fs, qts, vas, vb, fb, ql) {
  ql = ql || 7;
  if ([fs, qts, vas, vb, fb, ql].some(bad) || fs <= 0 || qts <= 0 || vas <= 0 || vb <= 0 || fb <= 0 || ql <= 0)
    return { error: 'Fs, Qts, Vas, box volume, tuning and QL must be positive numbers.' };
  var c = ventedCoeffs(fs, qts, vas, vb, fb, ql);
  // response curve, log grid 12 Hz..400 Hz
  var curve = [];
  var fmin = 12, fmax = 400, npts = 240;
  for (var i = 0; i < npts; i++) {
    var f = fmin * Math.pow(fmax / fmin, i / (npts - 1));
    var m = Math.sqrt(ventedMagSq(f / fb, c));
    curve.push({ f: f, dB: 20 * Math.log10(m) });
  }
  // peak (boom) over u in [0.4, 3]
  var peak = 0, peakU = 1;
  for (var j = 0; j <= 600; j++) {
    var u = 0.4 * Math.pow(3 / 0.4, j / 600);
    var ms = ventedMagSq(u, c);
    if (ms > peak) { peak = ms; peakU = u; }
  }
  var peakDb = 10 * Math.log10(peak);
  // f3/f6/f10: scan downward from u=2.5, first crossing of target magnitude
  function cut(target) {
    var prevU = 2.5, prevM = Math.sqrt(ventedMagSq(2.5, c));
    for (var k = 1; k <= 4000; k++) {
      var u2 = 2.5 * Math.pow(0.02 / 2.5, k / 4000);
      var m2 = Math.sqrt(ventedMagSq(u2, c));
      if (m2 < target && prevM >= target) {
        // bisect between prevU and u2
        var lo = u2, hi = prevU;
        for (var b = 0; b < 60; b++) {
          var mid = (lo + hi) / 2;
          if (Math.sqrt(ventedMagSq(mid, c)) >= target) hi = mid; else lo = mid;
        }
        return ((lo + hi) / 2) * fb;
      }
      prevU = u2; prevM = m2;
    }
    return null;
  }
  var f3 = cut(Math.SQRT1_2), f6 = cut(0.5), f10 = cut(Math.pow(10, -0.5));
  return { f3: f3, f6: f6, f10: f10, peakDb: peakDb, peakF: peakU * fb, curve: curve, coeffs: c };
}

// ---- port (vent) length from the Helmholtz resonator relation ----
// fb = (c/2pi) sqrt(Atot / (Vb * Leff));  Leff = L + endCorr * Dv  (endCorr 0.732, one flanged end)
function portLength(vbLiters, fb, diaMm, nPorts) {
  nPorts = nPorts || 1;
  if ([vbLiters, fb, diaMm, nPorts].some(bad) || vbLiters <= 0 || fb <= 0 || diaMm <= 0 || nPorts < 1)
    return { error: 'Box volume, tuning, port diameter and count must be positive numbers.' };
  var d = diaMm / 1000, a = Math.PI * d * d / 4, v = vbLiters / 1000;
  var leff = C_SOUND * C_SOUND * nPorts * a / (4 * Math.PI * Math.PI * fb * fb * v);
  var corr = 0.732 * d;
  var L = leff - corr;
  if (L <= 0) return { error: 'That port tunes the box with no physical length at all (it would need to be ' + (leff * 1000).toFixed(0) + ' mm deep including end correction). Use a smaller diameter or fewer ports.' };
  return { lengthMm: L * 1000, leffMm: leff * 1000, endCorrMm: corr * 1000, areaCm2: nPorts * a * 1e4 };
}

var engine = {
  sealedFromQtc: sealedFromQtc, sealedFromVb: sealedFromVb, sealedF3Ratio: sealedF3Ratio,
  keeleOptimum: keeleOptimum, keeleGivenVb: keeleGivenVb,
  ventedCoeffs: ventedCoeffs, ventedMagSq: ventedMagSq, ventedAnalyze: ventedAnalyze,
  portLength: portLength
};
if (typeof module !== 'undefined') module.exports = engine;

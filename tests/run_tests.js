// BoxTune test suite: engine vs independent Python oracle (tests/expected.json)
// plus published-reference checks (sealed f3 table, Butterworth B4 golden case).
const e = require('../engine.js');
const exp = require('./expected.json');
let pass = 0, fail = 0, fails = [];
function ok(cond, msg) { if (cond) pass++; else { fail++; fails.push(msg); } }
function near(a, b, tol, msg) { ok(typeof a === 'number' && Math.abs(a - b) <= tol, `${msg}: got ${a}, want ${b} +/- ${tol}`); }

// 1. published sealed f3/fc table (usenclosure / Dickason chart) - within published rounding
for (const t of exp.sealed_table) {
  const r = e.sealedF3Ratio(t.qtc);
  near(r, t.ratio, 1e-12, `sealed table qtc=${t.qtc} engine==oracle`);
  ok(Math.abs(r - t.published) <= 0.012, `sealed table qtc=${t.qtc} vs published ${t.published} (got ${r.toFixed(4)})`);
}

// 2. sealed designers vs oracle
const D = {
  budget10: { fs: 28.9, qts: 0.38, vas: 57.5 },
  compact65: { fs: 45, qts: 0.32, vas: 14 },
  pa12: { fs: 40, qts: 0.24, vas: 95 },
  b4textbook: { fs: 30, qts: 0.38268, vas: 42.47 }
};
for (const [n, d] of Object.entries(D)) {
  const a = e.sealedFromQtc(d.fs, d.qts, d.vas, 0.7071), o = exp.sealed_qtc[n];
  near(a.Vb, o.Vb, 0.02, `${n} sealedQtc Vb`); near(a.fc, o.fc, 0.02, `${n} sealedQtc fc`); near(a.f3, o.f3, 0.02, `${n} sealedQtc f3`);
  const b = e.sealedFromVb(d.fs, d.qts, d.vas, 20), p = exp.sealed_vb[n];
  near(b.Qtc, p.Qtc, 1e-9, `${n} sealedVb Qtc`); near(b.f3, p.f3, 0.02, `${n} sealedVb f3`);
  // identity: Vb from Qtc round-trips
  const rt = e.sealedFromVb(d.fs, d.qts, d.vas, a.Vb);
  near(rt.Qtc, 0.7071, 1e-9, `${n} sealed roundtrip Qtc`);
}

// 3. Keele approximations vs oracle
for (const [n, d] of Object.entries(D)) {
  const k = e.keeleOptimum(d.fs, d.qts, d.vas), o = exp.keele[n];
  near(k.Vb, o.Vb, 0.02, `${n} keele Vb`); near(k.Fb, o.Fb, 0.01, `${n} keele Fb`); near(k.F3, o.F3, 0.01, `${n} keele F3`);
  const g = e.keeleGivenVb(d.fs, d.qts, d.vas, 30), p = exp.keele_vb[n];
  near(g.F3, p.F3, 0.02, `${n} keeleVb F3`); near(g.Fb, p.Fb, 0.02, `${n} keeleVb Fb`);
}

// 4. vented analysis vs oracle (independent complex-arithmetic path)
const VC = {
  b4_lossless: [30, 0.38268, 42.47, 42.47 / 1.41421356, 30, 1e9],
  b4_ql7: [30, 0.38268, 42.47, 42.47 / 1.41421356, 30, 7],
  budget10_keele: [28.9, 0.38, 57.5, 53.7, 29, 7],
  boomy_small: [28.9, 0.38, 57.5, 18, 40, 7]
};
for (const [n, a] of Object.entries(VC)) {
  const r = e.ventedAnalyze(...a), o = exp.vented[n];
  near(r.f3, o.f3, 0.05, `${n} f3`); near(r.f6, o.f6, 0.05, `${n} f6`); near(r.f10, o.f10, 0.05, `${n} f10`);
  near(r.peakDb, o.peakDb, 0.01, `${n} peakDb`); near(r.peakF, o.peakF, 0.5, `${n} peakF`);
  // curve equality on 12 grid samples
  let maxDev = 0;
  o.samples.forEach((p, i) => { const q = r.curve[i * 20]; maxDev = Math.max(maxDev, Math.abs(p.f - q.f), Math.abs(p.dB - q.dB)); });
  ok(maxDev < 1e-9, `${n} curve sample match (maxDev ${maxDev})`);
}
// golden properties: lossless B4 is maximally flat with f3 = fs = fb
const b4 = e.ventedAnalyze(30, 0.38268, 42.47, 42.47 / 1.41421356, 30, 1e9);
near(b4.f3, 30, 0.05, 'B4 golden f3=fs');
ok(Math.abs(b4.peakDb) < 0.05, `B4 golden flat, ripple ${b4.peakDb.toFixed(4)} dB`);
// boomy small box must show a real peak above +0.5 dB
ok(e.ventedAnalyze(28.9, 0.38, 57.5, 18, 40, 7).peakDb > 0.5, 'boomy small box peak detected');

// 5. port lengths vs oracle + one hand-checked physics value
for (const [k, o] of Object.entries(exp.ports)) {
  const m = k.match(/([\d.]+)L_([\d.]+)Hz_([\d.]+)mm_x(\d+)/);
  const r = e.portLength(+m[1], +m[2], +m[3], +m[4]);
  near(r.lengthMm, o.lengthMm, 0.5, `port ${k} length`);
  near(r.areaCm2, o.areaCm2, 0.05, `port ${k} area`);
}
// hand check: 4" port (101.6mm), 56.6 L, 30 Hz -> Leff = c^2 A/(4 pi^2 fb^2 V)
const hand = e.portLength(56.6, 30, 101.6, 1);
const Leff = 343 * 343 * (Math.PI * 0.1016 * 0.1016 / 4) / (4 * Math.PI * Math.PI * 900 * 0.0566);
near(hand.leffMm, Leff * 1000, 0.01, 'port hand Leff');
near(hand.lengthMm, (Leff - 0.732 * 0.1016) * 1000, 0.01, 'port hand length');
// two ports of same total area => same tuning => 2x area needs 2x leff
const p1 = e.portLength(50, 32, 70, 1), p2 = e.portLength(50, 32, 70, 2);
near(p2.leffMm, 2 * p1.leffMm, 0.01, 'port N scaling');

// 6. error paths
ok(e.sealedFromQtc(30, 0.38, 57, 0.3).error, 'qtc<=qts rejected');
ok(e.sealedFromQtc(-1, 0.38, 57, 0.7).error, 'negative fs rejected');
ok(e.sealedFromVb(30, 0.38, 57, NaN).error, 'NaN vb rejected');
ok(e.keeleOptimum(30, 0, 57).error, 'zero qts rejected');
ok(e.ventedAnalyze(30, 0.38, 57, -5, 30).error, 'negative vb rejected');
ok(e.ventedAnalyze(30, 0.38, 57, 30, 0).error, 'zero fb rejected');
ok(e.portLength(50, 30, -10).error, 'negative port dia rejected');
ok(e.portLength(1, 200, 10, 1).error, 'absurd port (negative length) rejected');
ok(e.portLength(50, 30, 'abc').error, 'string port dia rejected');

console.log(`${pass} passed, ${fail} failed`);
if (fail) { console.log(fails.slice(0, 20).join('\n')); process.exit(1); }

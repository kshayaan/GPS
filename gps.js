/*
 * gps.js — the mechanism: how a GPS receiver works out where it is.
 *
 * Each satellite sends its position and the exact time it sent the signal.
 * The receiver notes when the signal arrived, by its own clock. Then:
 *
 *   distance = speed of light x (arrival time - send time)
 *
 * If the receiver's clock is wrong by dt, every distance is wrong by the same
 * amount, b = c x dt. These "wrong-by-the-same-amount" distances are called
 * pseudoranges:
 *
 *   rho_i = |satellite_i - receiver| + b
 *
 * The unknowns are the receiver's position (x, y, and z in 3D) and b.
 * With one more satellite than position unknowns, we can solve for all of
 * them. We solve by Gauss-Newton least squares: guess, measure how wrong the
 * guess is, correct, repeat.
 *
 * Units: metres for distances. Clock errors are stored as b (metres) and
 * converted to microseconds only for display. No UI code lives here.
 */
(function (root) {
  'use strict';

  var C = 299792458;          // speed of light, m/s
  var C_PER_US = C / 1e6;      // metres light travels in one microsecond (~299.79)

  function makeRng(seed) {
    var a = seed >>> 0;
    return function () {
      a = (a + 0x6d2b79f5) >>> 0;
      var t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // Normally distributed random number (Box-Muller), for measurement noise.
  function gaussian(rng) {
    var u = 1 - rng(), v = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  function dist(a, b) {
    var s = 0;
    for (var i = 0; i < a.length; i++) s += (a[i] - b[i]) * (a[i] - b[i]);
    return Math.sqrt(s);
  }

  // Solve A x = g for a small square system (Gaussian elimination with pivoting).
  function solveLinear(A, g) {
    var n = g.length, M = A.map(function (row, i) { return row.concat([g[i]]); });
    for (var c = 0; c < n; c++) {
      var piv = c;
      for (var r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
      if (Math.abs(M[piv][c]) < 1e-12) return null;
      var t = M[c]; M[c] = M[piv]; M[piv] = t;
      for (r = 0; r < n; r++) {
        if (r === c) continue;
        var f = M[r][c] / M[c][c];
        for (var k = c; k <= n; k++) M[r][k] -= f * M[c][k];
      }
    }
    return M.map(function (row, i) { return row[n] / row[i]; });
  }

  function invert(A) {
    var n = A.length, cols = [];
    for (var j = 0; j < n; j++) {
      var e = []; for (var i = 0; i < n; i++) e.push(i === j ? 1 : 0);
      var col = solveLinear(A, e);
      if (!col) return null;
      cols.push(col);
    }
    return A.map(function (_, i) { return cols.map(function (col) { return col[i]; }); });
  }

  /*
   * What the receiver measures: one pseudorange per satellite.
   *   truth      true receiver position
   *   clockUs    receiver clock error in microseconds
   *   extra      optional extra metres per satellite (atmosphere, a bounced signal)
   */
  function measure(sats, truth, clockUs, extra) {
    var b = (clockUs || 0) * C_PER_US;
    return sats.map(function (s, i) { return dist(s, truth) + b + ((extra && extra[i]) || 0); });
  }

  /*
   * Gauss-Newton fix.
   *   sats        satellite positions, each [x, y] or [x, y, z]
   *   rho         pseudoranges (metres)
   *   solveClock  true: treat the clock error as an unknown (the real GPS trick)
   *               false: trust the receiver's clock (b = 0)
   *   guess       starting position (default: centre of the satellites)
   */
  function solve(sats, rho, options) {
    var o = options || {};
    var dims = sats[0].length;
    var solveClock = o.solveClock !== false;
    var unknowns = dims + (solveClock ? 1 : 0);
    if (sats.length < unknowns) {
      return { ok: false, reason: 'too-few', unknowns: unknowns, measurements: sats.length };
    }
    var x = o.guess ? o.guess.slice() : sats[0].map(function (_, d) {
      return sats.reduce(function (s, p) { return s + p[d]; }, 0) / sats.length;
    });
    var b = 0, iterations = 0;
    for (; iterations < 60; iterations++) {
      var A = [], g = [];
      for (var i = 0; i < unknowns; i++) { A.push(new Array(unknowns).fill(0)); g.push(0); }
      sats.forEach(function (s, k) {
        var r = Math.max(1e-9, dist(x, s));
        var predicted = r + (solveClock ? b : 0);
        var e = rho[k] - predicted;
        // How the predicted range changes if we nudge each unknown.
        var J = x.map(function (xi, d) { return (xi - s[d]) / r; });
        if (solveClock) J.push(1);
        for (var a = 0; a < unknowns; a++) {
          g[a] += J[a] * e;
          for (var c2 = 0; c2 < unknowns; c2++) A[a][c2] += J[a] * J[c2];
        }
      });
      var step = solveLinear(A, g);
      if (!step) return { ok: false, reason: 'bad-geometry', unknowns: unknowns, measurements: sats.length };
      for (var d = 0; d < dims; d++) x[d] += step[d];
      if (solveClock) b += step[dims];
      var size = 0; for (d = 0; d < unknowns; d++) size += step[d] * step[d];
      if (Math.sqrt(size) < 1e-7) { iterations++; break; }
    }
    var residuals = sats.map(function (s, k) { return rho[k] - (dist(x, s) + b); });
    var rms = Math.sqrt(residuals.reduce(function (s, r) { return s + r * r; }, 0) / residuals.length);
    return { ok: true, pos: x, clockMetres: b, clockUs: b / C_PER_US, residuals: residuals, rms: rms,
      iterations: iterations, unknowns: unknowns, measurements: sats.length, spare: sats.length - unknowns };
  }

  // The 0, 1 or 2 points where two circles cross (2D). Used to show
  // "two satellites leave you in one of two places".
  function circleCrossings(c1, r1, c2, r2) {
    var dx = c2[0] - c1[0], dy = c2[1] - c1[1], d = Math.hypot(dx, dy);
    if (d === 0 || d > r1 + r2 || d < Math.abs(r1 - r2)) return [];
    var a = (r1 * r1 - r2 * r2 + d * d) / (2 * d);
    var h = Math.sqrt(Math.max(0, r1 * r1 - a * a));
    var mx = c1[0] + a * dx / d, my = c1[1] + a * dy / d;
    return [[mx + h * dy / d, my - h * dx / d], [mx - h * dy / d, my + h * dx / d]];
  }

  /*
   * Dilution of precision: how much the satellites' layout magnifies
   * measurement errors. 1-2 is excellent; above 6 is poor.
   */
  function dop(sats, pos, solveClock) {
    var useClock = solveClock !== false;
    var G = sats.map(function (s) {
      var r = dist(pos, s);
      var row = pos.map(function (p, d) { return (p - s[d]) / r; });
      if (useClock) row.push(1);
      return row;
    });
    var n = G[0].length;
    var GtG = [];
    for (var i = 0; i < n; i++) {
      GtG.push([]);
      for (var j = 0; j < n; j++) GtG[i].push(G.reduce(function (s, row) { return s + row[i] * row[j]; }, 0));
    }
    var Q = invert(GtG);
    if (!Q) return { ok: false };
    var dims = pos.length, p = 0, t = 0;
    for (i = 0; i < n; i++) { t += Q[i][i]; if (i < dims) p += Q[i][i]; }
    return { ok: true, position: Math.sqrt(p), total: Math.sqrt(t) };
  }

  /*
   * A cloud of fixes, each with fresh random noise on every distance.
   * Shows how far the blue dot wanders for a given satellite layout.
   */
  function noiseCloud(sats, truth, sigmaMetres, trials, seed) {
    var rng = makeRng(seed || 1), pts = [];
    for (var k = 0; k < trials; k++) {
      var noise = sats.map(function () { return gaussian(rng) * sigmaMetres; });
      var f = solve(sats, measure(sats, truth, 0, noise), { guess: truth.slice() });
      if (f.ok) pts.push(f.pos);
    }
    var errs = pts.map(function (p) { return dist(p, truth); }).sort(function (a, b) { return a - b; });
    return { points: pts, median: errs[Math.floor(errs.length / 2)] || 0, p95: errs[Math.floor(errs.length * 0.95)] || 0 };
  }

  /*
   * Catching a bad signal (for example one bounced off a building).
   * With spare satellites the distances can't all agree, so the leftover
   * error (rms) is large. Leave each satellite out in turn; the one whose
   * removal makes the rest agree best is the suspect.
   */
  function findBadSignal(sats, rho, thresholdMetres) {
    var full = solve(sats, rho);
    if (!full.ok) return { canCheck: false, full: full };
    if (full.spare < 1) return { canCheck: false, full: full, reason: 'no-spare' };
    var threshold = thresholdMetres || 5;
    if (full.rms < threshold) return { canCheck: true, flagged: false, full: full };
    if (full.spare < 2) return { canCheck: true, flagged: true, suspect: null, full: full, reason: 'cannot-tell-which' };
    var best = null;
    sats.forEach(function (_, k) {
      var s2 = sats.filter(function (_, j) { return j !== k; });
      var r2 = rho.filter(function (_, j) { return j !== k; });
      var f = solve(s2, r2);
      if (f.ok && (!best || f.rms < best.fix.rms)) best = { index: k, fix: f };
    });
    return { canCheck: true, flagged: true, suspect: best ? best.index : null, without: best ? best.fix : null, full: full };
  }

  /*
   * A signal in a street of tall buildings (a "street canyon").
   * The street runs east-west and is `streetWidth` metres wide. One row of buildings
   * (height southHeight) stands between you and the satellite, which is due south at
   * `elevationDeg` above the horizon; the row across the street has height northHeight.
   *
   * Direct signal: the straight line from you towards the satellite. It is blocked if it
   * passes the south buildings' wall lower than their roof.
   *
   * Bounced signal: a copy that reflects off the north buildings' wall. Mirror trick: a
   * bounce off a flat wall travels exactly as far as a straight line to your mirror image
   * behind that wall. Your mirror image is streetWidth further north, so the bounced
   * signal travels streetWidth x cos(elevation) further than the direct one would.
   *
   * GPS signals are radio waves, not visible light, but they travel at the speed of light
   * and bounce off walls the same way.
   */
  function streetCanyon(o) {
    var w = o.streetWidth, e = o.elevationDeg * Math.PI / 180;
    var hS = o.southHeight, hN = o.northHeight, h0 = o.phoneHeight || 1.5;
    var dir = [0, -Math.cos(e), Math.sin(e)];          // from you towards the satellite
    var you = [0, 0, h0];
    var rise = Math.tan(e);                             // metres up per metre across
    var directZ = h0 + (w / 2) * rise;                  // direct ray's height at the south wall
    var blockPoint = [0, -w / 2, directZ];
    var directBlocked = directZ < hS;
    var bouncePoint = [0, w / 2, h0 + (w / 2) * rise];  // where the bounce hits the north wall
    var hitsWall = bouncePoint[2] <= hN;                // the north wall is tall enough to reflect it
    var incomingZ = h0 + 1.5 * w * rise;                // bounced ray's height as it passes the south wall
    var incomingClear = incomingZ >= hS;
    var bounceOk = hitsWall && incomingClear;
    var extra = w * Math.cos(e);
    var heard = !directBlocked ? 'direct' : bounceOk ? 'bounced' : 'none';
    return { dir: dir, you: you, directZ: directZ, blockPoint: blockPoint, directBlocked: directBlocked,
      bouncePoint: bouncePoint, hitsWall: hitsWall, incomingZ: incomingZ, incomingClear: incomingClear,
      bounceOk: bounceOk, extra: extra, heard: heard, extraUsed: heard === 'bounced' ? extra : 0 };
  }

  var api = { C: C, C_PER_US: C_PER_US, makeRng: makeRng, gaussian: gaussian, dist: dist, measure: measure, solve: solve,
    circleCrossings: circleCrossings, dop: dop, noiseCloud: noiseCloud, findBadSignal: findBadSignal, streetCanyon: streetCanyon };
  root.GPS = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);

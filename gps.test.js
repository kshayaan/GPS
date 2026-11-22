// Run with: node gps.test.js
const assert = require('assert');
const G = require('./gps.js');

let passed = 0;
function test(name, fn) {
  try { fn(); passed++; console.log('  ok  ' + name); }
  catch (e) { console.log('  FAIL ' + name + '\n       ' + e.message); process.exitCode = 1; }
}

const sats2d = [[250, 250], [2750, 350], [1500, 1850], [350, 1700]];
const you = [1250, 900];

test('with a perfect clock and no noise, 3 satellites find you exactly', () => {
  const f = G.solve(sats2d.slice(0, 3), G.measure(sats2d.slice(0, 3), you, 0));
  assert.ok(f.ok);
  assert.ok(G.dist(f.pos, you) < 1e-6, 'off by ' + G.dist(f.pos, you));
});

test('the clock trick recovers a 1 microsecond clock error and the true position', () => {
  const rho = G.measure(sats2d.slice(0, 3), you, 1);
  const f = G.solve(sats2d.slice(0, 3), rho, { solveClock: true });
  assert.ok(G.dist(f.pos, you) < 1e-6);
  assert.ok(Math.abs(f.clockUs - 1) < 1e-9, 'clock ' + f.clockUs);
});

test('trusting a 1 microsecond wrong clock puts you over 100 m off', () => {
  const f = G.solve(sats2d.slice(0, 3), G.measure(sats2d.slice(0, 3), you, 1), { solveClock: false });
  assert.ok(f.ok);
  assert.ok(G.dist(f.pos, you) > 100, 'only ' + G.dist(f.pos, you));
});

test('too few satellites: refuses instead of guessing', () => {
  const f = G.solve(sats2d.slice(0, 2), G.measure(sats2d.slice(0, 2), you, 1), { solveClock: true });
  assert.strictEqual(f.ok, false);
  assert.strictEqual(f.reason, 'too-few');
  assert.strictEqual(f.unknowns, 3);
});

test('two circles cross at two points, both the right distance from each satellite', () => {
  const r1 = G.dist(sats2d[0], you), r2 = G.dist(sats2d[1], you);
  const pts = G.circleCrossings(sats2d[0], r1, sats2d[1], r2);
  assert.strictEqual(pts.length, 2);
  pts.forEach((p) => {
    assert.ok(Math.abs(G.dist(p, sats2d[0]) - r1) < 1e-6);
    assert.ok(Math.abs(G.dist(p, sats2d[1]) - r2) < 1e-6);
  });
  assert.ok(pts.some((p) => G.dist(p, you) < 1e-6), 'one crossing is you');
});

test('real scale, 3D: 4 satellites 20,200 km up recover position and clock', () => {
  const R = 6371000, orbit = R + 20200000;
  const receiver = [0, 0, R];
  const dirs = [[0.3, 0.2, 1], [-0.6, 0.4, 0.8], [0.5, -0.7, 0.6], [-0.2, -0.6, 0.9], [0.8, 0.5, 0.4]];
  const satPos = dirs.map((d) => { const n = Math.hypot(...d); return d.map((v) => v / n * orbit); });
  const rho = G.measure(satPos.slice(0, 4), receiver, 1000); // clock off by 1 ms
  const f = G.solve(satPos.slice(0, 4), rho, { guess: [0, 0, 0] });
  assert.ok(f.ok);
  assert.ok(G.dist(f.pos, receiver) < 1e-3, 'off by ' + G.dist(f.pos, receiver));
  assert.ok(Math.abs(f.clockUs - 1000) < 1e-6);
  assert.strictEqual(G.solve(satPos.slice(0, 3), rho.slice(0, 3)).ok, false, '3 satellites are not enough in 3D');
});

test('bunched satellites have worse geometry (higher DOP) than spread-out ones', () => {
  const spread = G.dop(sats2d, you, true).position;
  const bunched = G.dop([[200, 200], [500, 150], [350, 420], [150, 450]], you, true).position;
  assert.ok(bunched > 2 * spread, 'spread ' + spread + ', bunched ' + bunched);
});

test('the noise cloud is wider when satellites are bunched', () => {
  const a = G.noiseCloud(sats2d, you, 5, 300, 9).median;
  const b = G.noiseCloud([[200, 200], [500, 150], [350, 420], [150, 450]], you, 5, 300, 9).median;
  assert.ok(b > 2 * a, 'spread ' + a + ', bunched ' + b);
});

test('a bounced signal is caught and the right satellite is blamed (5 satellites)', () => {
  const sats = sats2d.concat([[2700, 1800]]);
  const extra = [0, 0, 0, 180, 0];
  const r = G.findBadSignal(sats, G.measure(sats, you, 0.5, extra));
  assert.ok(r.flagged);
  assert.strictEqual(r.suspect, 3);
  assert.ok(G.dist(r.without.pos, you) < 1e-6);
});

test('with no spare satellite a bounced signal cannot be detected', () => {
  const sats = sats2d.slice(0, 3);
  const r = G.findBadSignal(sats, G.measure(sats, you, 0, [0, 150, 0]));
  assert.strictEqual(r.canCheck, false);
  assert.ok(G.dist(r.full.pos, you) > 20, 'the fix is quietly wrong');
});

test('same seed, same noise', () => {
  const a = G.noiseCloud(sats2d, you, 5, 50, 3).median, b = G.noiseCloud(sats2d, you, 5, 50, 3).median;
  assert.strictEqual(a, b);
});

console.log('\n' + passed + ' tests passed');

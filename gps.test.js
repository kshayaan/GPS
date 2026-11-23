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

test('street canyon: an open street lets the direct signal through', () => {
  const r = G.streetCanyon({ streetWidth: 30, elevationDeg: 50, southHeight: 15, northHeight: 15 });
  assert.strictEqual(r.directBlocked, false);
  assert.strictEqual(r.heard, 'direct');
  assert.strictEqual(r.extraUsed, 0);
});

test('street canyon: a tall building blocks the direct signal but the bounce gets through', () => {
  const r = G.streetCanyon({ streetWidth: 50, elevationDeg: 30, southHeight: 40, northHeight: 70 });
  assert.strictEqual(r.directBlocked, true);
  assert.strictEqual(r.heard, 'bounced');
  assert.ok(Math.abs(r.extra - 50 * Math.cos(Math.PI / 6)) < 1e-9);
});

test('street canyon: the bounce path really is longer by streetWidth x cos(elevation)', () => {
  // Check the mirror-trick formula against walking the actual path: satellite -> wall -> you.
  const o = { streetWidth: 40, elevationDeg: 25, southHeight: 10, northHeight: 80 };
  const r = G.streetCanyon(o);
  const far = 1e7, sat = r.you.map((v, i) => v + r.dir[i] * far);   // a satellite very far away
  const direct = G.dist(sat, r.you);
  const bounced = G.dist(sat, r.bouncePoint) + G.dist(r.bouncePoint, r.you);
  assert.ok(Math.abs((bounced - direct) - r.extra) < 1e-3, 'walked ' + (bounced - direct) + ', formula ' + r.extra);
});

test('street canyon: a very tall building blocks the bounce too, so the signal is lost', () => {
  const r = G.streetCanyon({ streetWidth: 20, elevationDeg: 20, southHeight: 80, northHeight: 80 });
  assert.strictEqual(r.directBlocked, true);
  assert.strictEqual(r.incomingClear, false);
  assert.strictEqual(r.heard, 'none');
});

test('street canyon: a wall across the street that is too low cannot reflect the signal', () => {
  const r = G.streetCanyon({ streetWidth: 50, elevationDeg: 30, southHeight: 40, northHeight: 10 });
  assert.strictEqual(r.directBlocked, true);
  assert.strictEqual(r.hitsWall, false);
  assert.strictEqual(r.heard, 'none');
});

test('same seed, same noise', () => {
  const a = G.noiseCloud(sats2d, you, 5, 50, 3).median, b = G.noiseCloud(sats2d, you, 5, 50, 3).median;
  assert.strictEqual(a, b);
});

console.log('\n' + passed + ' tests passed');

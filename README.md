# The Stopwatch in Your Pocket: how GPS finds you

**The one idea:** your phone finds you with a stopwatch. It times how late each satellite's signal arrives, and each delay becomes a distance; the distances cross at one spot. The twist is that the stopwatch is cheap (one millionth of a second off makes every distance 300 m wrong), so GPS treats the stopwatch's error as one more unknown and solves for it too. This project lets you step through that on a map, break it, and fix it.

**Live site:** https://kshayaan.github.io/GPS/ (turn on GitHub Pages first: Settings → Pages → Deploy from a branch → main / root)

## Run it

Plain HTML, CSS and JavaScript. No libraries, no build step.

```
python3 -m http.server
# then open http://localhost:8000
```

## Test it

```
node gps.test.js
```

11 tests cover the solver: a perfect fix with a perfect clock; recovering a 1 µs clock error exactly; how far off you are when the clock is trusted; refusing to answer with too few satellites; the two crossing points of two circles; a real-scale 3D fix with satellites 20,200 km up and a clock 1 ms wrong; bunched satellites giving worse accuracy; and catching a bounced signal and blaming the right satellite.

## Files

| File | What it does |
|---|---|
| `gps.js` | The mechanism: pseudoranges, the Gauss–Newton position-and-clock solver (2D or 3D), circle crossings, dilution of precision, noise clouds, bad-signal detection. No UI code. |
| `index.html` | The page: intro cards, the 6-step story on a city map, free play with challenges, and three "go deeper" experiments. It only calls `gps.js` and draws what it returns. |
| `gps.test.js` | Tests for `gps.js`, run with Node. |

## How it works

1. Each satellite sends its position and the time it sent the signal.
2. The phone notes the arrival time by its own clock. Distance = speed of light × delay.
3. If the phone's clock is wrong by `dt`, every distance is wrong by the same `b = c × dt`. These measurements are called pseudoranges: `rho_i = |satellite_i − you| + b`.
4. The unknowns are your position (x, y, and height in 3D) and `b`. With as many satellites as unknowns, we solve by Gauss–Newton least squares: guess, compare predicted and measured distances, correct, repeat.
5. On the flat map that is 3 unknowns (x, y, clock), so 3 satellites. In the real 3D world it is 4.

## How the page is ordered

1. **What you already know:** open Maps and a blue dot appears.
2. **The two parts:** satellites that only send, and a GPS chip in your phone that only listens.
3. **How it works:** satellites shout the time; the phone's stopwatch turns each delay into a distance; three distances meet at you.
4. **Problem 1, inside the story:** the phone's stopwatch is cheap, so GPS solves for its error too.
5. **Why pinpointing you is hard:** the other problems, as experiments.

## The other problems (deeper experiments)

- **Geometry:** the same 5 m measurement errors give a fix within about 4 m when satellites are spread out, and about 100 m when they are bunched in one corner (dilution of precision).
- **A bounced signal:** with 5 satellites the distances can't all agree, and leaving each one out in turn finds the liar. With only 3 there is no spare, so the phone is quietly wrong.
- **Real scale:** 4 satellites 20,200 km up, a phone clock 1 ms wrong (300 km of error in every distance), and the solver still places you within a centimetre.

## What is computed, drawn, and simplified

Every distance, delay, circle, position, clock error, error cloud and table value is computed in the browser by `gps.js` from what the visitor presses. The streets, river and park, and the intro-card pictures, are hand-drawn decoration.

Simplified: the main map is flat, so 3 satellites are enough; satellites sit on the map's edges, not in orbit, and don't move; we leave out atmospheric delay corrections and each satellite's own clock correction; the real-scale section uses a fixed, made-up set of satellite directions.

## Where it falls short

_Write this yourselves before submitting: be honest about what you didn't get to (for example moving satellites, atmospheric delay, or a full 3D view)._

## References

1. GPS.gov, official US government information about the Global Positioning System.
2. Gauss–Newton method for non-linear least squares.
3. Dilution of precision (DOP), and receiver autonomous integrity monitoring (RAIM).

## How we built it

_Write this yourselves: who did what, which tools you used (including any AI tools) and how, and why you made the main choices._

## Licence

MIT. See [LICENSE](LICENSE).

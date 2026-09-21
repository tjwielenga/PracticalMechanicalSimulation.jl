import assert from "node:assert/strict";
import test from "node:test";
import { niceAxis, numberLabel } from "../src/plot.js";

test("uses clean multiples and includes zero when data crosses zero", () => {
  const axis = niceAxis([-0.3, 8.1]);
  assert.equal(axis.step, 2);
  assert.deepEqual(axis.ticks, [-2, 0, 2, 4, 6, 8, 10]);
});

test("uses the standard 1-2-5 progression at smaller scales", () => {
  const axis = niceAxis([-0.03, 0.081]);
  assert.equal(axis.step, 0.02);
  assert.deepEqual(axis.ticks,
    [-0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1]);
});

test("does not extend an axis to zero when positive data is well away", () => {
  assert.deepEqual(niceAxis([12, 20]).ticks, [12, 14, 16, 18, 20]);
});

test("formats tick labels using the selected step", () => {
  assert.equal(numberLabel(-2, 2), "-2");
  assert.equal(numberLabel(0, 2), "0");
  assert.equal(numberLabel(0.05, 0.05), "0.05");
  assert.equal(numberLabel(200000, 50000), "2e+5");
});

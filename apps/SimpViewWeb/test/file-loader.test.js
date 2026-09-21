import assert from "node:assert/strict";
import test from "node:test";
import { resampleSeries } from "../src/file-loader.js";

test("resamples stored signals onto a denser animation timeline", () => {
  const values = resampleSeries([0, 1, 2], [0, 10, 30],
    [0, 0.25, 1, 1.5, 2]);
  assert.deepEqual(values, [0, 2.5, 10, 20, 30]);
});

test("holds endpoint signal values outside the stored interval", () => {
  assert.deepEqual(resampleSeries([1, 2], [4, 8], [0, 1.5, 3]),
    [4, 6, 8]);
});

test("rejects signal values without corresponding source times", () => {
  assert.throws(() => resampleSeries([0, 1], [2], [0.5]),
    /matching source times/);
});

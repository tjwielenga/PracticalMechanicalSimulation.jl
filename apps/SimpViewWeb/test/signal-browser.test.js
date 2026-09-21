import assert from "node:assert/strict";
import test from "node:test";
import {
  optionTree,
  popupPlacement,
  signalOptions,
} from "../src/signal-browser.js";

test("builds assembly and component hierarchy from signal names", () => {
  const options = signalOptions([
    { name: "van.steering.pitman.theta" },
    { name: "van.steering.cross_link.V_x" },
    { name: "van.chassis.R_z" },
  ], true);
  const tree = optionTree(options);

  assert.equal(tree.leaves[0].name, "time (s)");
  const van = tree.groups.get("van");
  assert.ok(van);
  assert.equal(van.groups.get("chassis").leaves[0].label, "R_z");
  const steering = van.groups.get("steering");
  assert.equal(steering.groups.get("pitman").leaves[0].label, "theta");
  assert.equal(steering.groups.get("cross_link").leaves[0].label, "V_x");
});

test("gives time a distinct X-axis option", () => {
  const options = signalOptions([{ name: "body.R_x" }], true);
  assert.deepEqual(options.map((option) => option.id), ["time", "signal:0"]);
  assert.equal(options[0].index, null);
});

test("labels a static history axis as samples", () => {
  const options = signalOptions([], true, "static history sample");
  assert.equal(options[0].id, "time");
  assert.equal(options[0].name, "static history sample");
});

test("opens a tall signal menu toward the larger viewport area", () => {
  const nearBottom = popupPlacement({ left: 300, top: 650, bottom: 685 },
    1200, 800);
  assert.equal(nearBottom.above, true);
  assert.equal(nearBottom.top, null);
  assert.equal(nearBottom.bottom, 155);
  assert.equal(nearBottom.maxHeight, 432);

  const nearTop = popupPlacement({ left: 30, top: 40, bottom: 75 },
    600, 800);
  assert.equal(nearTop.above, false);
  assert.equal(nearTop.top, 80);
  assert.equal(nearTop.bottom, null);
});

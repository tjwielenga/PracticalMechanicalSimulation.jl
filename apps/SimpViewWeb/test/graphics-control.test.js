import assert from "node:assert/strict";
import test from "node:test";
import {
  defaultGraphicsNodeOpen,
  effectiveGraphicsScale,
  followTargetPath,
  graphicsPathKey,
  graphicsPathPrefixes,
  graphicsTree,
} from "../src/graphics-control.js";

test("builds graphics categories and assembly hierarchy", () => {
  const tree = graphicsTree([
    ["Model", "van", "Bodies", "chassis", "Geometry", "body"],
    ["Model", "van", "Bodies", "chassis", "Inertia"],
    ["Model", "van", "left_front", "Forces", "spring", "Reaction"],
    ["Model", "van", "left_front", "Joints", "upper_arm"],
    ["Model", "van", "rear", "left", "shackle", "Joints", "pivot"],
  ]);

  assert.equal(tree.name, "All graphics");
  const van = tree.children.get("Model").children.get("van");
  assert.ok(van.children.get("Bodies").children.get("chassis")
    .children.get("Geometry"));
  assert.ok(van.children.get("left_front").children.get("Forces")
    .children.get("spring"));
  assert.ok(van.children.get("left_front").children.get("Joints")
    .children.get("upper_arm"));
  assert.ok(van.children.get("rear").children.get("left")
    .children.get("shackle").children.get("Joints")
    .children.get("pivot"));
});

test("starts with every graphics hierarchy branch collapsed", () => {
  assert.equal(defaultGraphicsNodeOpen(["Follow"]), false);
  assert.equal(defaultGraphicsNodeOpen(["Model"]), false);
  assert.equal(defaultGraphicsNodeOpen(["Model", "van"]), false);
  assert.equal(defaultGraphicsNodeOpen(
    ["Model", "van", "rear", "left", "shackle"]), false);
  assert.equal(defaultGraphicsNodeOpen(
    ["Model", "van", "rear", "left", "shackle", "Joints"]), false);
  assert.equal(defaultGraphicsNodeOpen(
    ["Model", "van", "Bodies", "chassis"]), false);
});

test("places follow bodies in the assembly hierarchy", () => {
  assert.deepEqual(followTargetPath("van.front_left.spindle"),
    ["Follow", "van", "front_left", "spindle"]);
});

test("uses ground as the fixed follow reference", () => {
  assert.deepEqual(followTargetPath("ground"), ["Follow", "ground"]);
});

test("multiplies local graphics scales through the selected hierarchy", () => {
  const scales = new Map([
    [graphicsPathKey([]), 2],
    [graphicsPathKey(["Forces"]), 10],
    [graphicsPathKey(["Forces", "Reactions"]), 0.5],
  ]);
  const path = ["Forces", "Reactions", "van", "spring"];

  assert.deepEqual(graphicsPathPrefixes(path).slice(0, 3),
    [[], ["Forces"], ["Forces", "Reactions"]]);
  assert.equal(effectiveGraphicsScale(path, scales), 10);
  scales.delete(graphicsPathKey(["Forces"]));
  assert.equal(effectiveGraphicsScale(path, scales), 1);
});

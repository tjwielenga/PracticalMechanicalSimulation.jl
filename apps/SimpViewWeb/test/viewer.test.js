import assert from "node:assert/strict";
import test from "node:test";
import * as THREE from "three";
import {
  applyBodyFollow,
  defaultCameraConvention,
  interpolatedSample,
  modelGraphicPath,
  modeTrackBaseline,
  sceneFollowTargets,
} from "../src/viewer.js";

test("organizes stored graphics around model elements", () => {
  const bodies = ["van.chassis", "van.front.spindle"];
  assert.deepEqual(modelGraphicPath(
    ["Geometry", "van", "chassis", "graphics", "body"], bodies),
  ["Model", "van", "Bodies", "chassis", "Geometry", "body"]);
  assert.deepEqual(modelGraphicPath(
    ["Markers", "van", "front", "spindle", "cm"], bodies),
  ["Model", "van", "front", "Bodies", "spindle", "Markers", "cm"]);
  assert.deepEqual(modelGraphicPath(
    ["Forces", "Reactions", "van", "front", "spring"], bodies),
  ["Model", "van", "front", "Forces", "spring", "Reaction"]);
  assert.deepEqual(modelGraphicPath(
    ["Torques", "Applied", "steering"], bodies),
  ["Model", "Forces", "steering", "Applied torque"]);
  assert.deepEqual(modelGraphicPath(
    ["Bodies", "van", "chassis", "Geometry", "body"], bodies),
  ["Model", "van", "Bodies", "chassis", "Geometry", "body"]);
  assert.deepEqual(modelGraphicPath(
    ["Forces", "Applied", "gravity", "1"], bodies),
  ["Model", "Forces", "gravity", "Applied", "1"]);
  assert.deepEqual(modelGraphicPath(
    ["Torques", "Reactions", "van", "steering", "steering_gear",
      "coordinate_1"], bodies),
  ["Model", "van", "steering", "Forces", "steering_gear",
    "Reaction torque", "coordinate_1"]);
  assert.deepEqual(modelGraphicPath(
    ["Geometry", "Connectors", "van", "left_front_shock"], bodies),
  ["Model", "van", "Forces", "left_front_shock", "Default graphic"]);
  assert.deepEqual(modelGraphicPath(
    ["Geometry", "Belt spans", "drive", "upper_span"], bodies),
  ["Model", "drive", "Forces", "upper_span", "Default graphic"]);
  assert.deepEqual(modelGraphicPath(
    ["Joints", "van", "front", "upper_ball"], bodies),
  ["Model", "van", "front", "Joints", "upper_ball",
    "Default graphic"]);
});

test("uses y up while looking at the x-y plane for planar models", () => {
  assert.deepEqual(defaultCameraConvention("planar"), {
    up: [0, 1, 0],
    position: [0, 0, 1],
  });
});

test("interpolates display values between stored samples", () => {
  assert.equal(interpolatedSample([0, 4], 0.25), 1);
  assert.deepEqual(interpolatedSample([[0, 2, 4], [4, 6, 8]], 0.25),
    [1, 3, 5]);
});

test("finds the unperturbed modal pose from opposite phase samples", () => {
  const baseline = modeTrackBaseline({
    position: [[3, -1, 2], [1, 0, 2], [-1, 1, 2], [1, 0, 2],
      [3, -1, 2]],
    quaternion: [[0, 0, Math.sin(0.1), Math.cos(0.1)], [0, 0, 0, 1],
      [0, 0, -Math.sin(0.1), Math.cos(0.1)], [0, 0, 0, 1],
      [0, 0, Math.sin(0.1), Math.cos(0.1)]],
    scale: [[1.2, 0.8, 1], [1, 1, 1], [0.8, 1.2, 1], [1, 1, 1],
      [1.2, 0.8, 1]],
  });
  assert.deepEqual(baseline.position, [1, 0, 2]);
  assert.deepEqual(baseline.scale, [1, 1, 1]);
  assert.ok(Math.abs(baseline.quaternion[2]) < 1e-12);
  assert.ok(Math.abs(baseline.quaternion[3] - 1) < 1e-12);
});

test("retains the z-up oblique view for spatial models", () => {
  assert.deepEqual(defaultCameraConvention("spatial"), {
    up: [0, 0, 1],
    position: [4, -6, 3.5],
  });
});

test("uses explicit body-following tracks", () => {
  const targets = [{ name: "van.chassis", track: "body:van.chassis" }];
  assert.equal(sceneFollowTargets({ follow_targets: targets }), targets);
});

test("recognizes body tracks in existing viewer documents", () => {
  const scene = {
    tracks: [{ id: "surface:van.chassis.inertia_ellipsoid" }],
    instances: [{
      category: "inertia",
      track: "surface:van.chassis.inertia_ellipsoid",
    }],
  };
  assert.deepEqual(sceneFollowTargets(scene), [{
    name: "van.chassis",
    track: "surface:van.chassis.inertia_ellipsoid",
  }]);
});

test("translation-only follow preserves camera orientation and pan", () => {
  const camera = { position: new THREE.Vector3(3, 1, 2),
    up: new THREE.Vector3(0, 0, 1) };
  const orbitTarget = new THREE.Vector3(0, 1, 0);
  applyBodyFollow(camera, orbitTarget,
    new THREE.Vector3(0, 0, 0), new THREE.Quaternion(),
    new THREE.Vector3(5, 0, 0),
    new THREE.Quaternion().setFromAxisAngle(
      new THREE.Vector3(0, 0, 1), Math.PI / 2), false);
  assert.deepEqual(camera.position.toArray(), [8, 1, 2]);
  assert.deepEqual(orbitTarget.toArray(), [5, 1, 0]);
  assert.deepEqual(camera.up.toArray(), [0, 0, 1]);
});

test("rotating follow keeps the camera pose fixed relative to the body", () => {
  const camera = { position: new THREE.Vector3(3, 1, 2),
    up: new THREE.Vector3(0, 0, 1) };
  const orbitTarget = new THREE.Vector3(0, 1, 0);
  const rotation = new THREE.Quaternion().setFromAxisAngle(
    new THREE.Vector3(1, 0, 0), Math.PI / 2);
  const oldPosition = new THREE.Vector3(0, 0, 0);
  const newPosition = new THREE.Vector3(5, 0, 0);
  applyBodyFollow(camera, orbitTarget, oldPosition, new THREE.Quaternion(),
    newPosition, rotation, true);
  const cameraInBody = camera.position.clone().sub(newPosition)
    .applyQuaternion(rotation.clone().invert());
  const targetInBody = orbitTarget.clone().sub(newPosition)
    .applyQuaternion(rotation.clone().invert());
  const upInBody = camera.up.clone()
    .applyQuaternion(rotation.clone().invert());
  assert.ok(cameraInBody.distanceTo(new THREE.Vector3(3, 1, 2)) < 1e-12);
  assert.ok(targetInBody.distanceTo(new THREE.Vector3(0, 1, 0)) < 1e-12);
  assert.ok(upInBody.distanceTo(new THREE.Vector3(0, 0, 1)) < 1e-12);
});

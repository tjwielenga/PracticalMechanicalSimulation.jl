import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import {
  effectiveGraphicsScale,
  graphicsPathKey,
  graphicsPathPrefixes,
} from "./graphics-control.js";

const ENGINEERING_COLORS = {
  blue3: "#0000cd",
  darkorange2: "#ee7600",
  gold2: "#eeb422",
  gray25: "#404040",
  gray35: "#595959",
  gray60: "#999999",
  gray65: "#a6a6a6",
  gray70: "#b3b3b3",
  green3: "#00cd00",
  red3: "#cd0000",
};

function colorValue(name, fallback = "#6d7782") {
  if (!name) return fallback;
  const normalized = String(name).toLowerCase();
  const gray = /^gray(\d{1,3})$/.exec(normalized);
  if (gray && Number(gray[1]) <= 100) {
    const channel = Math.round(255 * Number(gray[1]) / 100)
      .toString(16).padStart(2, "0");
    return `#${channel}${channel}${channel}`;
  }
  return ENGINEERING_COLORS[normalized] ?? name;
}

function material(color, opacity = 1) {
  return new THREE.MeshStandardMaterial({
    color: colorValue(color),
    opacity,
    transparent: opacity < 1,
    roughness: 0.64,
    metalness: 0.04,
    side: THREE.DoubleSide,
  });
}

function lineMaterial(color, opacity = 1) {
  return new THREE.LineBasicMaterial({
    color: colorValue(color),
    opacity,
    transparent: opacity < 1,
  });
}

function vector(values) {
  return new THREE.Vector3(values[0], values[1], values[2]);
}

export function defaultCameraConvention(dimension) {
  return dimension === "planar"
    ? { up: [0, 1, 0], position: [0, 0, 1] }
    : { up: [0, 0, 1], position: [4, -6, 3.5] };
}

function interpolatedValue(first, second, fraction) {
  if (Array.isArray(first)) {
    return first.map((value, index) =>
      interpolatedValue(value, second[index], fraction));
  }
  return first + fraction * (second - first);
}

export function interpolatedSample(values, index) {
  const position = Math.max(0, Math.min(Number(index), values.length - 1));
  const lower = Math.floor(position);
  const upper = Math.min(lower + 1, values.length - 1);
  const fraction = position - lower;
  return fraction === 0 || lower === upper ? values[lower] :
    interpolatedValue(values[lower], values[upper], fraction);
}

const sample = interpolatedSample;

function normalizedQuaternion(values) {
  return new THREE.Quaternion(...values).normalize();
}

function interpolatedQuaternion(values, index) {
  const position = Math.max(0, Math.min(Number(index), values.length - 1));
  const lower = Math.floor(position);
  const upper = Math.min(lower + 1, values.length - 1);
  const fraction = position - lower;
  const first = normalizedQuaternion(values[lower]);
  if (fraction === 0 || lower === upper) return first.toArray();
  const second = normalizedQuaternion(values[upper]);
  if (first.dot(second) < 0) {
    second.set(-second.x, -second.y, -second.z, -second.w);
  }
  return first.slerp(second, fraction).normalize().toArray();
}

export function modeTrackBaseline(track) {
  const last = Math.max((track.position?.length ?? 1) - 1, 0);
  const opposite = Math.round(last / 2);
  const midpoint = (first, second) => first.map((value, index) =>
    0.5 * (value + second[index]));
  const firstQuaternion = normalizedQuaternion(sample(track.quaternion, 0));
  const oppositeQuaternion = normalizedQuaternion(
    sample(track.quaternion, opposite));
  if (firstQuaternion.dot(oppositeQuaternion) < 0) {
    oppositeQuaternion.set(-oppositeQuaternion.x, -oppositeQuaternion.y,
      -oppositeQuaternion.z, -oppositeQuaternion.w);
  }
  return {
    position: midpoint(sample(track.position, 0),
      sample(track.position, opposite)),
    quaternion: firstQuaternion.clone().slerp(oppositeQuaternion, 0.5)
      .toArray(),
    scale: midpoint(sample(track.scale, 0), sample(track.scale, opposite)),
  };
}

function scaledQuaternion(baseline, current, scale) {
  const base = normalizedQuaternion(baseline);
  const value = normalizedQuaternion(current);
  if (base.dot(value) < 0) {
    value.set(-value.x, -value.y, -value.z, -value.w);
  }
  const delta = base.clone().invert().multiply(value).normalize();
  const sine = Math.hypot(delta.x, delta.y, delta.z);
  if (sine < 1e-12) return base;
  const angle = 2 * Math.atan2(sine, Math.max(-1, Math.min(1, delta.w)));
  const axis = new THREE.Vector3(delta.x, delta.y, delta.z)
    .multiplyScalar(1 / sine);
  return base.multiply(new THREE.Quaternion().setFromAxisAngle(axis,
    scale * angle)).normalize();
}

function setCylinder(anchor, mesh, pointA, pointB, radius) {
  const first = vector(pointA);
  const second = vector(pointB);
  const direction = second.clone().sub(first);
  const length = Math.max(direction.length(), 1e-8);
  anchor.position.copy(first.add(second).multiplyScalar(0.5));
  mesh.scale.set(radius, length, radius);
  anchor.quaternion.setFromUnitVectors(
    new THREE.Vector3(0, 1, 0),
    direction.normalize(),
  );
}

function defaultGraphicPath(group, name, subgroup = null) {
  const path = [group];
  if (subgroup) path.push(subgroup);
  if (name) path.push(...String(name).split("."));
  return path;
}

export function modelGraphicPath(path, bodyNames = []) {
  const parts = path.map(String);
  if (parts.length === 0 || parts[0] === "Model") return parts;
  const bodyPath = (ownerParts, group, remainder = []) => {
    if (ownerParts.length === 1 && ownerParts[0].toLowerCase() === "ground") {
      return ["Model", "ground", group, ...remainder];
    }
    return ["Model", ...ownerParts.slice(0, -1), "Bodies",
      ownerParts.at(-1), group, ...remainder];
  };
  const elementPath = (group, nameParts, suffix = []) => [
    "Model", ...nameParts.slice(0, -1), group, nameParts.at(-1), ...suffix,
  ];
  const group = parts[0];
  if (group === "Bodies") {
    const categoryIndex = parts.findIndex((part, index) => index > 0 &&
      ["Geometry", "Inertia", "Markers", "Frames", "Deformation"]
        .includes(part));
    if (categoryIndex > 1) {
      return bodyPath(parts.slice(1, categoryIndex), parts[categoryIndex],
        parts.slice(categoryIndex + 1));
    }
  }
  if (["Joints", "Measurements"].includes(group)) {
    const name = parts.slice(1);
    const suffix = group === "Joints" && name.at(-1) !== "Default graphic"
      ? ["Default graphic"] : [];
    return elementPath(group, name, suffix);
  }
  if (group === "Geometry" &&
      ["Connectors", "Belt spans"].includes(parts[1])) {
    return elementPath("Forces", parts.slice(2), ["Default graphic"]);
  }
  if (["Geometry", "Inertia", "Markers", "Frames"].includes(group)) {
    const candidates = [...bodyNames, "ground"].sort(
      (first, second) => second.length - first.length);
    const remainderText = parts.slice(1).join(".");
    const owner = candidates.find((name) => remainderText === name ||
      remainderText.startsWith(`${name}.`));
    if (owner) {
      const ownerParts = owner.split(".");
      const remainder = parts.slice(1 + ownerParts.length);
      for (let index = remainder.length - 1; index >= 0; index -= 1) {
        if (remainder[index] === "graphics") remainder.splice(index, 1);
      }
      return bodyPath(ownerParts, group, remainder);
    }
    return elementPath(group, parts.slice(1));
  }
  if (["Forces", "Torques", "Bushings"].includes(group)) {
    let kind = group === "Bushings" ? "Default graphic" : null;
    let name = parts.slice(1);
    if (["Applied", "Reactions"].includes(name[0])) {
      const oldKind = name.shift();
      kind = group === "Torques"
        ? `${oldKind === "Applied" ? "Applied" : "Reaction"} torque`
        : oldKind === "Applied" ? "Applied" : "Reaction";
    } else {
      const knownKind = ["Applied", "Reaction", "Applied torque",
        "Reaction torque", "Bushing", "Connector", "Belt span",
        "Default graphic"];
      if (knownKind.includes(name.at(-1))) kind = name.pop();
    }
    if (["Bushing", "Connector", "Belt span"].includes(kind)) {
      kind = "Default graphic";
    }
    const generated = name.length > 1 &&
      (/^\d+$/.test(name.at(-1)) || /^coordinate_\d+$/.test(name.at(-1)))
      ? [name.pop()] : [];
    return elementPath("Forces", name,
      [...(kind ? [kind] : []), ...generated]);
  }
  return parts;
}

function categoryObject(object, category) {
  object.userData.viewerCategory = category;
  object.userData.sampleVisible = true;
  return object;
}

function meshObject(asset, instance) {
  const surfaceMaterial = material(instance.color, instance.opacity);
  if (asset.kind === "cylinder") {
    return new THREE.Mesh(new THREE.CylinderGeometry(1, 1, 1, 20),
      surfaceMaterial);
  }
  if (asset.kind === "sphere") {
    return new THREE.Mesh(new THREE.SphereGeometry(0.5, 24, 16),
      surfaceMaterial);
  }
  if (asset.kind === "box") {
    return new THREE.Mesh(new THREE.BoxGeometry(1, 1, 1), surfaceMaterial);
  }
  if (asset.kind === "frustum") {
    return new THREE.Mesh(new THREE.CylinderGeometry(asset.radius_b,
      asset.radius_a, 1, asset.segments ?? 24), surfaceMaterial);
  }
  if (asset.kind === "surface") {
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.Float32BufferAttribute(
      asset.vertices.flat(), 3));
    geometry.setIndex(asset.faces.flat().map((index) => index - 1));
    geometry.computeVertexNormals();
    return new THREE.Mesh(geometry, surfaceMaterial);
  }
  if (asset.kind === "lines") {
    const positions = [];
    for (const edge of asset.edges) {
      positions.push(...asset.vertices[edge[0] - 1],
        ...asset.vertices[edge[1] - 1]);
    }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position",
      new THREE.Float32BufferAttribute(positions, 3));
    return new THREE.LineSegments(geometry,
      lineMaterial(instance.color, instance.opacity));
  }
  if (asset.kind === "arrow") {
    const arrow = new THREE.Group();
    const shaft = new THREE.Mesh(new THREE.CylinderGeometry(0.055, 0.055,
      0.68, 12), surfaceMaterial);
    shaft.position.y = 0.34;
    const head = new THREE.Mesh(new THREE.ConeGeometry(0.14, 0.32, 16),
      surfaceMaterial);
    head.position.y = 0.84;
    arrow.add(shaft, head);
    return arrow;
  }
  if (asset.kind === "torque") {
    const torque = new THREE.Group();
    const shaft = new THREE.Mesh(new THREE.BoxGeometry(0.11, 0.58, 0.11),
      surfaceMaterial);
    shaft.position.y = 0.29;
    const first = new THREE.Mesh(new THREE.ConeGeometry(0.16, 0.22, 16),
      surfaceMaterial);
    first.position.y = 0.68;
    const second = new THREE.Mesh(new THREE.ConeGeometry(0.11, 0.18, 16),
      surfaceMaterial);
    second.position.y = 0.88;
    torque.add(shaft, first, second);
    return torque;
  }
  throw new Error(`Unsupported viewer mesh kind: ${asset.kind}`);
}

function arrayMaximum(items, field) {
  let maximum = 0;
  for (const item of items ?? []) {
    for (const row of item[field] ?? []) {
      maximum = Math.max(maximum, vector(row).length());
    }
  }
  return maximum;
}

export function sceneFollowTargets(scene) {
  const explicit = scene.follow_targets ?? [];
  if (explicit.length > 0) return explicit;
  const trackIds = new Set((scene.tracks ?? []).map((track) => track.id));
  const inferred = new Map();
  for (const instance of scene.instances ?? []) {
    if (instance.category !== "inertia" || !trackIds.has(instance.track)) {
      continue;
    }
    const path = instance.path ?? [];
    if (path[0] === "Inertia" && path.length >= 2) {
      inferred.set(path.slice(1).join("."), instance.track);
      continue;
    }
    const prefix = "surface:";
    const suffix = ".inertia_ellipsoid";
    if (instance.track.startsWith(prefix) && instance.track.endsWith(suffix)) {
      inferred.set(instance.track.slice(prefix.length, -suffix.length),
        instance.track);
    }
  }
  return [...inferred].map(([name, track]) => ({ name, track }));
}

export function applyBodyFollow(camera, orbitTarget, previousPosition,
    previousQuaternion, position, quaternion, rotateWithBody) {
  if (!rotateWithBody) {
    const translation = position.clone().sub(previousPosition);
    camera.position.add(translation);
    orbitTarget.add(translation);
    return;
  }
  const rotation = quaternion.clone().multiply(
    previousQuaternion.clone().invert()).normalize();
  camera.position.sub(previousPosition).applyQuaternion(rotation).add(position);
  orbitTarget.sub(previousPosition).applyQuaternion(rotation).add(position);
  camera.up.applyQuaternion(rotation).normalize();
}

export class MechanismScene {
  constructor(canvas) {
    this.canvas = canvas;
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color("#f4f6f8");
    this.camera = new THREE.PerspectiveCamera(38, 1, 0.001, 1e6);
    const initialCamera = defaultCameraConvention("spatial");
    this.camera.up.fromArray(initialCamera.up);
    this.camera.position.fromArray(initialCamera.position);
    this.controls = this.createControls();
    this.root = new THREE.Group();
    this.scene.add(this.root);
    this.updaters = [];
    this.graphicObjects = [];
    this.graphicScales = new Map();
    this.graphicVisibility = new Map();
    this.deformationScales = new Map();
    this.deformationGroups = new Set();
    this.modelBodyNames = [];
    this.followTargets = new Map();
    this.followTargetName = null;
    this.followPosition = new THREE.Vector3();
    this.followQuaternion = new THREE.Quaternion();
    this.followRotation = false;
    this.dimension = "spatial";
    this.modal = false;
    this.modeScale = 1;
    this.sampleIndex = 0;

    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x75808b, 2.2));
    const light = new THREE.DirectionalLight(0xffffff, 2.6);
    light.position.set(4, -5, 8);
    this.scene.add(light);

    this.resizeObserver = new ResizeObserver(() => this.resize());
    this.resizeObserver.observe(canvas.parentElement);
    this.animate = this.animate.bind(this);
    requestAnimationFrame(this.animate);
  }

  createControls() {
    const controls = new OrbitControls(this.camera, this.canvas);
    controls.enableDamping = true;
    controls.dampingFactor = 0.08;
    controls.screenSpacePanning = true;
    return controls;
  }

  animate() {
    this.controls.update();
    this.renderer.render(this.scene, this.camera);
    requestAnimationFrame(this.animate);
  }

  resize() {
    const bounds = this.canvas.parentElement.getBoundingClientRect();
    if (bounds.width < 2 || bounds.height < 2) return;
    this.renderer.setSize(bounds.width, bounds.height, false);
    this.camera.aspect = bounds.width / bounds.height;
    this.camera.updateProjectionMatrix();
  }

  clear() {
    this.root.traverse((object) => {
      object.geometry?.dispose?.();
      if (Array.isArray(object.material)) {
        object.material.forEach((item) => item.dispose?.());
      } else {
        object.material?.dispose?.();
      }
    });
    this.scene.remove(this.root);
    this.root = new THREE.Group();
    this.scene.add(this.root);
    this.updaters = [];
    this.graphicObjects = [];
    this.deformationGroups.clear();
    this.modelBodyNames = [];
    this.followTargets.clear();
    this.followTargetName = null;
  }

  registerGraphic(object, path, includeInFit = true) {
    const normalized = modelGraphicPath(path, this.modelBodyNames);
    object.userData.viewerPath = normalized;
    object.userData.viewerBaseScale = object.scale.clone();
    object.userData.viewerIncludeInFit = includeInFit;
    this.graphicObjects.push(object);
    return object;
  }

  graphicPaths() {
    return this.graphicObjects.map((object) => object.userData.viewerPath);
  }

  deformationPaths() {
    return [...this.deformationGroups]
      .sort((first, second) => first.localeCompare(second, undefined, {
        numeric: true,
        sensitivity: "base",
      }))
      .map((name) => {
        const parts = String(name).split(".");
        return ["Model", ...parts.slice(0, -1), "Bodies", parts.at(-1),
          "Deformation"];
      });
  }

  registerFollowTarget(name, object) {
    this.followTargets.set(String(name), object);
  }

  followTargetNames() {
    return [...this.followTargets.keys()].sort((first, second) =>
      first.localeCompare(second, undefined, {
        numeric: true,
        sensitivity: "base",
      }));
  }

  getFollowTarget() {
    return this.followTargetName;
  }

  getFollowRotation() {
    return this.followRotation;
  }

  setFollowRotation(enabled) {
    this.followRotation = Boolean(enabled);
  }

  setFollowTarget(name) {
    if (name === null) {
      this.followTargetName = null;
      return;
    }
    const normalized = String(name);
    const target = this.followTargets.get(normalized);
    if (!target) throw new Error(`Unknown follow body: ${normalized}`);
    this.followTargetName = normalized;
    target.getWorldPosition(this.followPosition);
    target.getWorldQuaternion(this.followQuaternion);
  }

  addCylinder(trajectory, category = "geometry", color = "#687785",
    opacity = 1, path = null) {
    const anchor = categoryObject(new THREE.Group(), category);
    const mesh = new THREE.Mesh(
      new THREE.CylinderGeometry(1, 1, 1, 20),
      material(trajectory.color ?? color, trajectory.opacity ?? opacity),
    );
    anchor.add(mesh);
    this.root.add(anchor);
    this.registerGraphic(anchor, path ?? defaultGraphicPath(
      category === "joints" ? "Joints" : "Geometry", trajectory.name));
    this.updaters.push((index) => setCylinder(anchor, mesh,
      sample(trajectory.point_a, index), sample(trajectory.point_b, index),
      trajectory.radius ?? 0.01));
  }

  addFrustum(trajectory, path = null) {
    const anchor = categoryObject(new THREE.Group(), "geometry");
    const mesh = new THREE.Mesh(
      new THREE.CylinderGeometry(trajectory.radius_b, trajectory.radius_a,
        1, trajectory.kind === "gear" ? 32 : 24),
      material(trajectory.color, trajectory.opacity),
    );
    anchor.add(mesh);
    this.root.add(anchor);
    this.registerGraphic(anchor, path ?? defaultGraphicPath(
      "Geometry", trajectory.name));
    this.updaters.push((index) => setCylinder(anchor, mesh,
      sample(trajectory.point_a, index), sample(trajectory.point_b, index), 1));
  }

  addSurface(surface) {
    if (surface.encoding === "rigid_pose") {
      this.addRigidSurface(surface);
      return;
    }
    for (const patch of surface.patches) {
      const geometry = new THREE.BufferGeometry();
      const positions = new Float32Array(surface.vertices[0].length * 3);
      geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
      geometry.setIndex(patch.faces.flat().map((index) => index - 1));
      geometry.computeVertexNormals();
      const mesh = categoryObject(new THREE.Mesh(geometry,
        material(patch.color, patch.opacity)),
      surface.category === "inertia" ? "inertia" : "geometry");
      this.root.add(mesh);
      this.registerGraphic(mesh, defaultGraphicPath(
        surface.category === "inertia" ? "Inertia" : "Geometry",
        surface.name, patch.name), surface.include_in_fit ?? true);
      this.updaters.push((index) => {
        const vertices = sample(surface.vertices, index);
        geometry.attributes.position.array.set(vertices.flat());
        geometry.attributes.position.needsUpdate = true;
        geometry.computeVertexNormals();
        geometry.computeBoundingSphere();
      });
    }
    if (surface.edges.length > 0) {
      const positions = new Float32Array(surface.edges.length * 6);
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
      const lines = categoryObject(new THREE.LineSegments(geometry,
        lineMaterial(surface.edge_color)),
      surface.category === "inertia" ? "inertia" : "geometry");
      this.root.add(lines);
      this.registerGraphic(lines, defaultGraphicPath(
        surface.category === "inertia" ? "Inertia" : "Geometry",
        surface.name, "Edges"), surface.include_in_fit ?? true);
      this.updaters.push((index) => {
        const vertices = sample(surface.vertices, index);
        const values = [];
        for (const edge of surface.edges) {
          values.push(...vertices[edge[0] - 1], ...vertices[edge[1] - 1]);
        }
        geometry.attributes.position.array.set(values);
        geometry.attributes.position.needsUpdate = true;
      });
    }
  }

  addRigidSurface(surface) {
    const category = surface.category === "inertia" ? "inertia" : "geometry";
    const group = categoryObject(new THREE.Group(), category);
    this.root.add(group);
    this.registerGraphic(group, defaultGraphicPath(
      category === "inertia" ? "Inertia" : "Geometry", surface.name),
      surface.include_in_fit ?? true);
    for (const patch of surface.patches) {
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute("position", new THREE.Float32BufferAttribute(
        surface.vertices.flat(), 3));
      geometry.setIndex(patch.faces.flat().map((index) => index - 1));
      geometry.computeVertexNormals();
      group.add(new THREE.Mesh(geometry, material(patch.color, patch.opacity)));
    }
    if (surface.edges.length > 0) {
      const positions = [];
      for (const edge of surface.edges) {
        positions.push(...surface.vertices[edge[0] - 1],
          ...surface.vertices[edge[1] - 1]);
      }
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute("position",
        new THREE.Float32BufferAttribute(positions, 3));
      group.add(new THREE.LineSegments(geometry,
        lineMaterial(surface.edge_color)));
    }
    this.updaters.push((index) => {
      group.position.copy(vector(sample(surface.position, index)));
      group.quaternion.fromArray(
        interpolatedQuaternion(surface.quaternion, index));
    });
  }

  addMarker(marker) {
    const geometry = marker.shape === "box"
      ? new THREE.BoxGeometry(1, 1, 1)
      : new THREE.SphereGeometry(0.5, 20, 12);
    const mesh = categoryObject(new THREE.Mesh(geometry,
      material(marker.color, marker.opacity)),
    marker.category === "marker" ? "markers" : "geometry");
    mesh.scale.set(...marker.size);
    this.root.add(mesh);
    this.registerGraphic(mesh, defaultGraphicPath(
      marker.category === "marker" ? "Markers" : "Geometry", marker.name));
    this.updaters.push((index) => {
      mesh.position.copy(vector(sample(marker.center, index)));
      mesh.rotation.set(0, 0, sample(marker.angle, index));
    });
  }

  addJoint(joint) {
    const sphere = categoryObject(new THREE.Mesh(
      new THREE.SphereGeometry(joint.diameter / 2, 18, 12),
      material("#b4bbc2")), "joints");
    this.root.add(sphere);
    this.registerGraphic(sphere, defaultGraphicPath("Joints", joint.name));
    this.updaters.push((index) => sphere.position.copy(
      vector(sample(joint.position, index))));
  }

  addBody(body, color) {
    this.addCylinder(body, "geometry", color, 1,
      defaultGraphicPath("Geometry", body.name));
    const ellipsoid = categoryObject(new THREE.Mesh(
      new THREE.SphereGeometry(0.5, 24, 16),
      material(color, 0.35)), "inertia");
    ellipsoid.scale.set(...body.ellipsoid_axes);
    this.root.add(ellipsoid);
    this.registerGraphic(ellipsoid, defaultGraphicPath("Inertia", body.name));
    this.registerFollowTarget(body.name, ellipsoid);
    this.updaters.push((index) => {
      ellipsoid.position.copy(vector(sample(body.center, index)));
      ellipsoid.rotation.set(0, 0, sample(body.angle, index));
    });
  }

  addPlanarWheel(item, color, opacity) {
    const trajectory = {
      point_a: item.center.map((center) =>
        [center[0], center[1], center[2] - item.half_width]),
      point_b: item.center.map((center) =>
        [center[0], center[1], center[2] + item.half_width]),
      radius: item.radius,
    };
    this.addCylinder(trajectory, "geometry", color, opacity,
      defaultGraphicPath("Geometry", item.name));
  }

  addFrame(frame) {
    const isJoint = frame.category === "perp";
    const category = frame.category === "marker" ? "markers"
      : isJoint ? "joints" : "frames";
    const axes = isJoint
      ? [[frame.x_direction, "#9a9a9a"], [frame.y_direction, "#9a9a9a"]]
      : [[frame.x_direction, "#bd3131"], [frame.y_direction, "#27914d"],
        [frame.z_direction, "#315dbd"]];
    for (const [directions, color] of axes) {
      const arrow = categoryObject(new THREE.ArrowHelper(
        new THREE.Vector3(1, 0, 0), new THREE.Vector3(), frame.axis_length,
        color, 0.22 * frame.axis_length, 0.08 * frame.axis_length,
      ), category);
      this.root.add(arrow);
      this.registerGraphic(arrow, defaultGraphicPath(
        isJoint ? "Joints" : category === "markers" ? "Markers" : "Frames",
        frame.name));
      this.updaters.push((index) => {
        arrow.position.copy(vector(sample(frame.origin, index)));
        arrow.setDirection(vector(sample(directions, index)).normalize());
      });
    }
  }

  addForce(arrow, appearance, lengthScale) {
    const category = arrow.category === "torque" ? "torques" : "loads";
    const color = arrow.category === "reaction"
      ? appearance.reaction_color : appearance.applied_color;
    const object = categoryObject(new THREE.ArrowHelper(
      new THREE.Vector3(0, 0, 1), new THREE.Vector3(), 1,
      colorValue(color), 0.2, 0.08,
    ), category);
    this.root.add(object);
    this.registerGraphic(object, defaultGraphicPath("Forces", arrow.name,
      arrow.category === "reaction" ? "Reactions" : "Applied"));
    this.updaters.push((index) => {
      const value = vector(sample(arrow.force, index));
      const length = value.length() * lengthScale;
      object.userData.sampleVisible = length > 1e-12;
      object.position.copy(vector(sample(arrow.position, index)));
      if (length > 1e-12) {
        object.setDirection(value.normalize());
        object.setLength(length, Math.min(0.22 * length, 0.08),
          Math.min(0.09 * length, 0.035));
      }
    });
  }

  addTorque(arrow, appearance, lengthScale) {
    const color = arrow.category === "reaction"
      ? appearance.reaction_color : appearance.applied_color;
    const object = categoryObject(new THREE.ArrowHelper(
      new THREE.Vector3(0, 0, 1), new THREE.Vector3(), 1,
      colorValue(color), 0.2, 0.08,
    ), "torques");
    this.root.add(object);
    this.registerGraphic(object, defaultGraphicPath("Torques", arrow.name,
      arrow.category === "reaction" ? "Reactions" : "Applied"));
    this.updaters.push((index) => {
      const magnitude = sample(arrow.torque, index);
      const axis = vector(sample(arrow.axis, index));
      const length = Math.abs(magnitude) * lengthScale;
      object.userData.sampleVisible = length > 1e-12;
      object.position.copy(vector(sample(arrow.position, index)));
      if (length > 1e-12) {
        if (magnitude < 0) axis.negate();
        object.setDirection(axis.normalize());
        object.setLength(length, Math.min(0.22 * length, 0.08),
          Math.min(0.09 * length, 0.035));
      }
    });
  }

  loadSceneGraph(scene) {
    const meshes = scene.meshes ?? {};
    const followTargets = sceneFollowTargets(scene);
    this.modelBodyNames = followTargets.map((target) => String(target.name));
    const trackGroups = new Map();
    const trackLookup = new Map((scene.tracks ?? [])
      .map((track) => [track.id, track]));
    const baselines = new Map((scene.tracks ?? []).map((track) => [track.id,
      this.modal && (track.position?.length ?? 0) > 1
        ? modeTrackBaseline(track) : null]));
    const trackPose = (track, index) => {
      const position = vector(sample(track.position, index));
      const quaternion = new THREE.Quaternion().fromArray(
        interpolatedQuaternion(track.quaternion, index));
      const scale = vector(sample(track.scale, index));
      const baseline = baselines.get(track.id);
      if (!baseline) return { position, quaternion, scale };
      const basePosition = vector(baseline.position);
      const baseScale = vector(baseline.scale);
      return {
        position: basePosition.clone().add(
          position.sub(basePosition).multiplyScalar(this.modeScale)),
        quaternion: scaledQuaternion(baseline.quaternion,
          quaternion.toArray(), this.modeScale),
        scale: baseScale.clone().add(
          scale.sub(baseScale).multiplyScalar(this.modeScale)),
      };
    };
    const deformationReferencePose = (deformation, index) => {
      const referenceTrack = trackLookup.get(deformation.reference_track);
      if (!referenceTrack) return null;
      const parent = trackPose(referenceTrack, index);
      const localPosition = vector(deformation.local_position ?? [0, 0, 0])
        .multiply(parent.scale).applyQuaternion(parent.quaternion);
      const localQuaternion = new THREE.Quaternion().fromArray(
        deformation.local_quaternion ?? [0, 0, 0, 1]);
      return {
        position: parent.position.clone().add(localPosition),
        quaternion: parent.quaternion.clone().multiply(localQuaternion),
        scale: parent.scale.clone().multiply(
          vector(deformation.local_scale ?? [1, 1, 1])),
      };
    };
    for (const track of scene.tracks ?? []) {
      const group = new THREE.Group();
      if (track.deformation?.group) {
        this.deformationGroups.add(String(track.deformation.group));
      }
      this.root.add(group);
      trackGroups.set(track.id, group);
      this.updaters.push((index) => {
        const pose = trackPose(track, index);
        const deformation = track.deformation;
        const reference = deformation
          ? deformationReferencePose(deformation, index) : null;
        if (reference) {
          const factor = this.getDeformationScale(deformation.group);
          group.position.copy(reference.position.clone().add(
            pose.position.sub(reference.position).multiplyScalar(factor)));
          group.quaternion.copy(scaledQuaternion(
            reference.quaternion.toArray(), pose.quaternion.toArray(), factor));
          group.scale.copy(reference.scale.clone().add(
            pose.scale.sub(reference.scale).multiplyScalar(factor)));
        } else {
          group.position.copy(pose.position);
          group.quaternion.copy(pose.quaternion);
          group.scale.copy(pose.scale);
        }
      });
    }
    for (const target of followTargets) {
      const group = trackGroups.get(target.track);
      if (group) this.registerFollowTarget(target.name, group);
    }
    for (const instance of scene.instances ?? []) {
      const parent = trackGroups.get(instance.track);
      const asset = meshes[instance.mesh];
      if (!parent || !asset) continue;
      const object = categoryObject(meshObject(asset, instance),
        instance.category);
      object.position.fromArray(instance.local_position ?? [0, 0, 0]);
      object.quaternion.fromArray(instance.local_quaternion ?? [0, 0, 0, 1]);
      object.scale.fromArray(instance.local_scale ?? [1, 1, 1]);
      parent.add(object);
      this.registerGraphic(object, instance.path ?? defaultGraphicPath(
        instance.category === "loads" ? "Forces"
          : instance.category === "torques" ? "Torques"
            : instance.category === "joints" ? "Joints"
              : instance.category === "markers" ? "Markers"
                : instance.category === "frames" ? "Frames"
                  : instance.category === "inertia" ? "Inertia" : "Geometry",
        instance.name), instance.include_in_fit ?? true);
    }
    for (const surface of scene.dynamic_surfaces ?? []) {
      this.addSurface(surface);
    }
  }

  load(choice, appearance = {}, dimension = "spatial", resetView = false,
      modal = false) {
    this.clear();
    this.modal = modal;
    if (resetView || dimension !== this.dimension) {
      const convention = defaultCameraConvention(dimension);
      this.dimension = dimension;
      const nextUp = vector(convention.up);
      const upChanged = !this.camera.up.equals(nextUp);
      this.camera.up.fromArray(convention.up);
      this.camera.position.fromArray(convention.position);
      if (upChanged) {
        this.controls.dispose();
        this.controls = this.createControls();
      }
      this.controls.target.set(0, 0, 0);
      this.controls.update();
    }
    this.scene.background = new THREE.Color(colorValue(appearance.background,
      "#f4f6f8"));
    let scene = choice.scene;
    if (scene.encoding === "mesh_instances") {
      this.loadSceneGraph(scene);
      scene = {};
    }
    const palette = appearance.body_palette ??
      ["steelblue", "darkorange", "seagreen", "orchid"];
    if (scene.bodies) {
      this.modelBodyNames = scene.bodies.map((body) => String(body.name));
    }
    for (const [index, item] of (scene.bodies ?? []).entries()) {
      this.addBody(item, palette[index % palette.length]);
    }
    for (const [index, item] of (scene.gears ?? []).entries()) {
      this.addPlanarWheel(item, palette[index % palette.length],
        item.internal ? 0.22 : 0.72);
    }
    for (const [index, item] of (scene.pulleys ?? []).entries()) {
      this.addPlanarWheel(item, palette[index % palette.length], 0.68);
    }
    for (const item of scene.graphic_cylinders ?? []) this.addCylinder(item);
    for (const item of scene.graphic_frustums ?? []) this.addFrustum(item);
    for (const item of scene.graphic_surfaces ?? []) this.addSurface(item);
    for (const item of scene.graphic_markers ?? []) this.addMarker(item);
    for (const item of scene.xy_frames ?? []) this.addFrame(item);
    for (const item of scene.joints ?? []) this.addJoint(item);
    for (const item of scene.guides ?? []) {
      this.addCylinder(item, "joints", "#aeb6be", 1,
        defaultGraphicPath("Joints", item.name));
    }
    for (const item of scene.connectors ?? []) {
      const bushing = item.kind === "bushing";
      this.addCylinder(item, "geometry", "#25282b", 1,
        defaultGraphicPath(bushing ? "Bushings" : "Geometry", item.name,
          bushing ? null : "Connectors"));
    }
    for (const item of scene.belt_spans ?? []) {
      this.addCylinder({ ...item, point_a: item.point_1,
        point_b: item.point_2, radius: 0.008 }, "geometry", "#25282b", 1,
      defaultGraphicPath("Geometry", item.name, "Belt spans"));
    }
    for (const item of scene.spheres ?? []) {
      this.addJoint({ position: item.center, diameter: 2 * item.radius });
    }
    this.update(0);
    const box = this.fitBox();
    const span = box.isEmpty() ? 1 : box.getSize(new THREE.Vector3()).length();
    const forceMaximum = arrayMaximum(scene.force_arrows, "force");
    const torqueMaximum = Math.max(...(scene.torque_arrows ?? [])
      .flatMap((item) => item.torque.map(Math.abs)), 0);
    const forceScale = forceMaximum > 0 ? 0.2 * span / forceMaximum : 1;
    const torqueScale = torqueMaximum > 0 ? 0.16 * span / torqueMaximum : 1;
    for (const item of scene.force_arrows ?? []) {
      this.addForce(item, appearance, forceScale);
    }
    for (const item of scene.torque_arrows ?? []) {
      this.addTorque(item, appearance, torqueScale);
    }
    this.update(0);
    this.applyGraphicSettings();
    this.fit();
  }

  update(index) {
    this.sampleIndex = index;
    for (const update of this.updaters) update(index);
    const followed = this.followTargets.get(this.followTargetName);
    if (followed) {
      const position = followed.getWorldPosition(new THREE.Vector3());
      const quaternion = followed.getWorldQuaternion(new THREE.Quaternion());
      applyBodyFollow(this.camera, this.controls.target,
        this.followPosition, this.followQuaternion,
        position, quaternion, this.followRotation);
      this.followPosition.copy(position);
      this.followQuaternion.copy(quaternion);
    }
    this.applyGraphicSettings();
  }

  setModeScale(scale) {
    this.modeScale = Number.isFinite(scale) && scale > 0 ? scale : 1;
    if (this.modal) this.update(this.sampleIndex);
  }

  getGraphicScale(path) {
    return this.graphicScales.get(graphicsPathKey(path)) ?? 1;
  }

  getEffectiveGraphicScale(path) {
    return effectiveGraphicsScale(path, this.graphicScales);
  }

  setGraphicScale(path, scale) {
    this.graphicScales.set(graphicsPathKey(path), scale);
    this.applyGraphicSettings();
  }

  resetGraphicScale(path) {
    this.graphicScales.delete(graphicsPathKey(path));
    this.applyGraphicSettings();
  }

  resetAllGraphicScales() {
    this.graphicScales.clear();
    this.deformationScales.clear();
    this.update(this.sampleIndex);
  }

  getDeformationScale(group) {
    return this.deformationScales.get(String(group)) ?? 1;
  }

  setDeformationScale(group, scale) {
    this.deformationScales.set(String(group), scale);
    this.update(this.sampleIndex);
  }

  resetDeformationScale(group) {
    this.deformationScales.delete(String(group));
    this.update(this.sampleIndex);
  }

  getGraphicVisible(path) {
    return this.graphicVisibility.get(graphicsPathKey(path)) ?? true;
  }

  setGraphicVisible(path, visible) {
    this.graphicVisibility.set(graphicsPathKey(path), visible);
    this.applyGraphicSettings();
  }

  effectiveGraphicVisible(path) {
    return graphicsPathPrefixes(path).every((prefix) =>
      this.getGraphicVisible(prefix));
  }

  resetGraphicSettings() {
    this.graphicScales.clear();
    this.graphicVisibility.clear();
    this.deformationScales.clear();
    this.applyGraphicSettings();
  }

  applyGraphicSettings() {
    for (const object of this.graphicObjects) {
      const path = object.userData.viewerPath;
      const base = object.userData.viewerBaseScale;
      object.scale.copy(base).multiplyScalar(
        this.getEffectiveGraphicScale(path));
      object.visible = this.effectiveGraphicVisible(path) &&
        object.userData.sampleVisible !== false;
    }
  }

  fitBox() {
    const box = new THREE.Box3();
    this.root.updateWorldMatrix(true, true);
    for (const object of this.graphicObjects) {
      if (!object.visible || object.userData.viewerIncludeInFit === false) {
        continue;
      }
      box.union(new THREE.Box3().setFromObject(object));
    }
    return box;
  }

  fit() {
    const box = this.fitBox();
    if (box.isEmpty()) return;
    const center = box.getCenter(new THREE.Vector3());
    const size = box.getSize(new THREE.Vector3());
    const radius = Math.max(size.length() * 0.55, 0.1);
    const direction = this.camera.position.clone()
      .sub(this.controls.target).normalize();
    this.controls.target.copy(center);
    this.camera.position.copy(center).add(direction.multiplyScalar(2.7 * radius));
    this.camera.near = Math.max(radius / 10000, 1e-5);
    this.camera.far = Math.max(radius * 1000, 100);
    this.camera.updateProjectionMatrix();
    this.controls.update();
  }
}

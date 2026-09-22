const PATH_SEPARATOR = "\u001f";
const FOLLOW_GROUP = "Follow";
const ROOT_ORDER = [FOLLOW_GROUP, "Model"];
const MODEL_GROUP_ORDER = [
  "Bodies", "Joints", "Forces", "Measurements", "Geometry", "Inertia",
  "Markers", "Frames", "Deformation",
];

export function graphicsPathKey(path) {
  return path.join(PATH_SEPARATOR);
}

export function graphicsPathPrefixes(path) {
  const prefixes = [[]];
  for (let length = 1; length <= path.length; length += 1) {
    prefixes.push(path.slice(0, length));
  }
  return prefixes;
}

export function followTargetPath(name) {
  return [FOLLOW_GROUP, ...String(name).split(".")];
}

export function effectiveGraphicsScale(path, scales) {
  return graphicsPathPrefixes(path).reduce((value, prefix) =>
    value * (scales.get(graphicsPathKey(prefix)) ?? 1), 1);
}

export function defaultGraphicsNodeOpen(path) {
  return false;
}

export function graphicsTree(paths) {
  const root = { name: "All graphics", path: [], children: new Map() };
  for (const path of paths) {
    let node = root;
    for (const name of path) {
      if (!node.children.has(name)) {
        const childPath = [...node.path, name];
        node.children.set(name, {
          name,
          path: childPath,
          children: new Map(),
        });
      }
      node = node.children.get(name);
    }
  }
  return root;
}

function orderedChildren(node) {
  return [...node.children.values()].sort((first, second) => {
    if (node.path.length === 0) {
      const firstIndex = ROOT_ORDER.indexOf(first.name);
      const secondIndex = ROOT_ORDER.indexOf(second.name);
      if (firstIndex >= 0 || secondIndex >= 0) {
        return (firstIndex < 0 ? ROOT_ORDER.length : firstIndex) -
          (secondIndex < 0 ? ROOT_ORDER.length : secondIndex);
      }
    }
    if ((node.path.length === 1 && node.path[0] === FOLLOW_GROUP) ||
        node.path[0] === "Model") {
      const firstGround = first.name.toLowerCase() === "ground";
      const secondGround = second.name.toLowerCase() === "ground";
      if (firstGround !== secondGround) return firstGround ? -1 : 1;
    }
    if (node.path[0] === "Model") {
      const firstIndex = MODEL_GROUP_ORDER.indexOf(first.name);
      const secondIndex = MODEL_GROUP_ORDER.indexOf(second.name);
      if ((firstIndex >= 0) !== (secondIndex >= 0)) {
        return firstIndex >= 0 ? 1 : -1;
      }
      if (firstIndex >= 0 && secondIndex >= 0) {
        return firstIndex - secondIndex;
      }
    }
    return first.name.localeCompare(second.name, undefined, {
      numeric: true,
      sensitivity: "base",
    });
  });
}

function formatScale(value) {
  if (value >= 100 || value < 0.01) return `${value.toPrecision(3)}×`;
  if (value >= 10 || value < 0.1) return `${value.toFixed(2)}×`;
  return `${value.toFixed(3).replace(/0+$/, "").replace(/\.$/, "")}×`;
}

export class GraphicsControl {
  constructor(root, callbacks) {
    this.root = root;
    this.callbacks = callbacks;
    this.treeElement = root.querySelector(".graphics-tree");
    this.rotationToggle = root.querySelector(".follow-rotation-toggle");
    this.selectedName = root.querySelector(".graphics-selected-name");
    this.scaleSlider = root.querySelector(".graphics-scale-slider");
    this.localScale = root.querySelector(".graphics-local-scale");
    this.effectiveScale = root.querySelector(".graphics-effective-scale");
    this.localScaleLabel = root.querySelector(".graphics-local-scale-label");
    this.effectiveScaleLabel = root.querySelector(
      ".graphics-effective-scale-label");
    this.resetButton = root.querySelector(".graphics-reset");
    this.resetAllButton = root.querySelector(".graphics-reset-all");
    this.selectedPath = [];
    this.paths = [];
    this.followPaths = [];
    this.followTargets = new Map();
    this.deformationPaths = [];
    this.deformationKeys = new Set();
    this.visibilityInputs = [];
    this.followInputs = [];
    this.openPaths = null;

    this.scaleSlider.addEventListener("input", () => {
      const scale = 10 ** Number(this.scaleSlider.value);
      if (this.isDeformation(this.selectedPath)) {
        this.callbacks.setDeformationScale(
          this.deformationGroup(this.selectedPath), scale);
      } else {
        this.callbacks.setScale(this.selectedPath, scale);
      }
      this.refreshValues();
    });
    this.scaleSlider.addEventListener("dblclick", () => this.resetSelected());
    this.resetButton.addEventListener("click", () => this.resetSelected());
    this.resetAllButton.addEventListener("click", () => {
      this.callbacks.resetAllScales();
      this.refreshValues();
    });
    this.rotationToggle.addEventListener("change", () => {
      this.callbacks.setFollowRotation(this.rotationToggle.checked);
    });
  }

  setPaths(paths, followTargets = [], deformationPaths = []) {
    const unique = new Map();
    for (const path of paths) unique.set(graphicsPathKey(path), path);
    this.paths = [...unique.values()];
    this.deformationPaths = deformationPaths.map((path) => [...path]);
    this.deformationKeys = new Set(this.deformationPaths.map(graphicsPathKey));
    const available = new Set([""]);
    for (const path of this.paths) {
      for (const prefix of graphicsPathPrefixes(path)) {
        available.add(graphicsPathKey(prefix));
      }
    }
    for (const path of this.deformationPaths) {
      for (const prefix of graphicsPathPrefixes(path)) {
        available.add(graphicsPathKey(prefix));
      }
    }
    this.followTargets.clear();
    const groundPath = [FOLLOW_GROUP, "Ground"];
    this.followPaths = [groundPath];
    this.followTargets.set(graphicsPathKey(groundPath), null);
    for (const target of followTargets) {
      const path = followTargetPath(target);
      this.followPaths.push(path);
      this.followTargets.set(graphicsPathKey(path), String(target));
    }
    if (!available.has(graphicsPathKey(this.selectedPath))) {
      this.selectedPath = [];
    }
    this.render();
    this.refreshFollow();
    this.refreshValues();
  }

  resetSelected() {
    if (this.isDeformation(this.selectedPath)) {
      this.callbacks.resetDeformationScale(
        this.deformationGroup(this.selectedPath));
    } else {
      this.callbacks.resetScale(this.selectedPath);
    }
    this.refreshValues();
  }

  isDeformation(path) {
    return this.deformationKeys.has(graphicsPathKey(path));
  }

  deformationGroup(path) {
    const bodies = path.indexOf("Bodies");
    if (path[0] === "Model" && bodies > 0 && bodies + 1 < path.length) {
      return [...path.slice(1, bodies), path[bodies + 1]].join(".");
    }
    return path.slice(1, -1).join(".");
  }

  select(path) {
    this.selectedPath = [...path];
    this.render();
    this.refreshValues();
  }

  follow(target) {
    this.callbacks.setFollowTarget(target);
    this.refreshFollow();
  }

  renderNode(node, container, rootNode = false) {
    const row = document.createElement("div");
    row.className = "graphics-tree-row";
    if (graphicsPathKey(node.path) === graphicsPathKey(this.selectedPath)) {
      row.classList.add("selected");
    }
    const nodeKey = graphicsPathKey(node.path);
    const followNode = node.path[0] === FOLLOW_GROUP;
    const deformationNode = this.isDeformation(node.path);
    const hasFollowTarget = this.followTargets.has(nodeKey);
    const control = document.createElement(hasFollowTarget ||
      (!followNode && !deformationNode) ? "input" : "span");
    if (hasFollowTarget) {
      control.type = "radio";
      control.name = "follow-body";
      control.title = node.name === "Ground" ? "Follow ground" :
        `Follow ${this.followTargets.get(nodeKey)}`;
      control.addEventListener("change", () => {
        if (control.checked) this.follow(this.followTargets.get(nodeKey));
      });
      this.followInputs.push({
        input: control,
        row,
        target: this.followTargets.get(nodeKey),
      });
    } else if (followNode || deformationNode) {
      control.className = "graphics-tree-spacer";
    } else {
      control.type = "checkbox";
      control.checked = this.callbacks.getVisible(node.path);
      control.title = `Show ${node.name}`;
      control.addEventListener("change", () => {
        this.callbacks.setVisible(node.path, control.checked);
        this.refreshVisibility();
        this.refreshValues();
      });
      this.visibilityInputs.push({ checkbox: control, path: node.path });
    }
    const label = document.createElement("button");
    label.type = "button";
    label.textContent = node.name;
    label.title = node.path.join(".") || "All graphics";
    if (hasFollowTarget) {
      label.addEventListener("click", (event) => {
        event.preventDefault();
        event.stopPropagation();
        this.follow(this.followTargets.get(nodeKey));
      });
    } else if (followNode) {
      label.addEventListener("click", (event) => {
        event.preventDefault();
        event.stopPropagation();
      });
    } else {
      label.addEventListener("click", (event) => {
        event.preventDefault();
        event.stopPropagation();
        this.select(node.path);
      });
    }
    row.append(control, label);

    if (rootNode || node.children.size === 0) {
      container.append(row);
    } else {
      const details = document.createElement("details");
      details.graphicsPathKey = nodeKey;
      details.open = (this.openPaths === null
        ? defaultGraphicsNodeOpen(node.path)
        : this.openPaths.has(nodeKey));
      const summary = document.createElement("summary");
      summary.append(row);
      details.append(summary);
      const children = document.createElement("div");
      children.className = "graphics-tree-children";
      for (const child of orderedChildren(node)) {
        this.renderNode(child, children);
      }
      details.append(children);
      container.append(details);
    }
    if (rootNode) {
      const children = document.createElement("div");
      children.className = "graphics-tree-root-children";
      for (const child of orderedChildren(node)) {
        this.renderNode(child, children);
      }
      container.append(children);
    }
  }

  render() {
    const existingDetails = this.treeElement.querySelectorAll("details");
    if (existingDetails.length > 0) {
      this.openPaths = new Set([...existingDetails]
        .filter((details) => details.open)
        .map((details) => details.graphicsPathKey));
    }
    const scrollTop = this.treeElement.scrollTop;
    this.visibilityInputs = [];
    this.followInputs = [];
    this.treeElement.replaceChildren();
    const tree = graphicsTree([
      ...this.paths, ...this.followPaths, ...this.deformationPaths]);
    const children = document.createElement("div");
    children.className = "graphics-tree-root-children";
    for (const child of orderedChildren(tree)) {
      this.renderNode(child, children);
    }
    this.treeElement.append(children);
    this.treeElement.scrollTop = scrollTop;
  }

  refreshVisibility() {
    for (const { checkbox, path } of this.visibilityInputs) {
      checkbox.checked = this.callbacks.getVisible(path);
    }
  }

  refreshFollow() {
    const followed = this.callbacks.getFollowTarget();
    this.rotationToggle.checked = this.callbacks.getFollowRotation();
    this.rotationToggle.disabled = followed === null;
    for (const { input, row, target } of this.followInputs) {
      const selected = followed === target;
      input.checked = selected;
      row.classList.toggle("following", selected);
    }
  }

  refreshValues() {
    const enabled = this.paths.length > 0;
    this.scaleSlider.disabled = !enabled;
    this.resetButton.disabled = !enabled;
    this.resetAllButton.disabled = !enabled;
    this.selectedName.textContent = this.selectedPath.join(" › ") ||
      "All graphics";
    const deformation = this.isDeformation(this.selectedPath);
    const local = deformation
      ? this.callbacks.getDeformationScale(
        this.deformationGroup(this.selectedPath))
      : this.callbacks.getScale(this.selectedPath);
    const effective = deformation ? local :
      this.callbacks.getEffectiveScale(this.selectedPath);
    this.localScaleLabel.firstChild.textContent = deformation
      ? "Amplification " : "Local ";
    this.effectiveScaleLabel.hidden = deformation;
    this.scaleSlider.value = String(Math.log10(local));
    this.localScale.value = formatScale(local);
    this.effectiveScale.value = formatScale(effective);
  }
}

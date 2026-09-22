function signalPath(name) {
  return String(name).split(".").filter(Boolean);
}

export function signalOptions(signals, includeTime = false,
    timeLabel = "time (s)") {
  const options = signals.map((signal, index) => ({
    id: `signal:${index}`,
    index,
    name: signal.name,
    path: signalPath(signal.name),
  }));
  if (includeTime) {
    options.unshift({
      id: "time",
      index: null,
      name: timeLabel,
      path: [timeLabel],
    });
  }
  return options;
}

export function optionTree(options) {
  const root = { groups: new Map(), leaves: [] };
  for (const option of options) {
    let node = root;
    for (const part of option.path.slice(0, -1)) {
      if (!node.groups.has(part)) {
        node.groups.set(part, { groups: new Map(), leaves: [] });
      }
      node = node.groups.get(part);
    }
    node.leaves.push({ ...option, label: option.path.at(-1) });
  }
  return root;
}

function appendTree(parent, node, onSelect, expand = false) {
  const leaves = [...node.leaves]
    .sort((first, second) => first.label.localeCompare(second.label));
  const appendLeaf = (leaf) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "signal-leaf";
    button.dataset.optionId = leaf.id;
    button.textContent = leaf.label;
    button.title = leaf.name;
    button.addEventListener("click", () => onSelect(leaf));
    parent.append(button);
  };
  for (const leaf of leaves.filter((item) => item.id === "time")) {
    appendLeaf(leaf);
  }
  for (const [name, child] of [...node.groups.entries()]
    .sort(([first], [second]) => first.localeCompare(second))) {
    const details = document.createElement("details");
    details.open = expand;
    const summary = document.createElement("summary");
    summary.textContent = name;
    details.append(summary);
    const contents = document.createElement("div");
    contents.className = "signal-tree-children";
    appendTree(contents, child, onSelect, expand);
    details.append(contents);
    parent.append(details);
  }
  for (const leaf of leaves.filter((item) => item.id !== "time")) {
    appendLeaf(leaf);
  }
}

export function popupPlacement(bounds, viewportWidth, viewportHeight) {
  const margin = 8;
  const gap = 5;
  const width = Math.min(480, Math.max(80, viewportWidth - 2 * margin));
  const left = Math.max(margin, Math.min(bounds.left + 24,
    viewportWidth - margin - width));
  const availableAbove = Math.max(0, bounds.top - gap - margin);
  const availableBelow = Math.max(0,
    viewportHeight - bounds.bottom - gap - margin);
  const above = availableAbove >= availableBelow;
  const available = above ? availableAbove : availableBelow;
  return {
    above,
    left,
    width,
    maxHeight: Math.min(432, available),
    top: above ? null : bounds.bottom + gap,
    bottom: above ? viewportHeight - bounds.top + gap : null,
  };
}

export class SignalBrowser {
  constructor(root, { includeTime = false, onSelect = () => {} } = {}) {
    this.root = root;
    this.includeTime = includeTime;
    this.onSelect = onSelect;
    this.button = root.querySelector(".signal-browser-button");
    this.popup = root.querySelector(".signal-browser-popup");
    this.search = root.querySelector(".signal-browser-search");
    this.tree = root.querySelector(".signal-tree");
    this.options = [];
    this.selected = null;
    this.positionPopup = this.positionPopup.bind(this);

    this.button.addEventListener("click", () => this.toggle());
    this.search.addEventListener("input", () => this.render());
    this.search.addEventListener("keydown", (event) => {
      if (event.key === "Escape") this.close();
    });
    document.addEventListener("pointerdown", (event) => {
      if (!this.root.contains(event.target)) this.close();
    });
  }

  setSignals(signals, preferredId = null, timeLabel = "time (s)") {
    this.options = signalOptions(signals, this.includeTime, timeLabel);
    const selected = this.options.find((option) => option.id === preferredId) ??
      this.options[0] ?? null;
    this.select(selected, false);
    this.button.disabled = this.options.length === 0;
    this.render();
  }

  select(option, notify = true) {
    if (!option) {
      this.selected = null;
      this.button.textContent = "No variables";
      return;
    }
    this.selected = option;
    this.button.textContent = option.name;
    this.button.title = option.name;
    this.close();
    this.renderSelection();
    if (notify) this.onSelect(option);
  }

  renderSelection() {
    for (const leaf of this.tree.querySelectorAll(".signal-leaf")) {
      leaf.classList.toggle("selected",
        leaf.dataset.optionId === this.selected?.id);
    }
  }

  render() {
    const query = this.search.value.trim().toLowerCase();
    const visible = query
      ? this.options.filter((option) => option.name.toLowerCase().includes(query))
      : this.options;
    this.tree.replaceChildren();
    if (visible.length === 0) {
      const empty = document.createElement("div");
      empty.className = "signal-tree-empty";
      empty.textContent = "No matching variables";
      this.tree.append(empty);
      return;
    }
    appendTree(this.tree, optionTree(visible), (option) => this.select(option),
      query.length > 0);
    this.renderSelection();
  }

  toggle() {
    if (this.popup.hidden) this.open();
    else this.close();
  }

  open() {
    this.popup.hidden = false;
    this.positionPopup();
    window.addEventListener("resize", this.positionPopup);
    this.button.setAttribute("aria-expanded", "true");
    this.search.focus();
    this.search.select();
  }

  close() {
    this.popup.hidden = true;
    window.removeEventListener("resize", this.positionPopup);
    this.button.setAttribute("aria-expanded", "false");
  }

  positionPopup() {
    if (this.popup.hidden) return;
    const placement = popupPlacement(this.root.getBoundingClientRect(),
      window.innerWidth, window.innerHeight);
    this.popup.style.left = `${placement.left}px`;
    this.popup.style.width = `${placement.width}px`;
    this.popup.style.maxHeight = `${placement.maxHeight}px`;
    this.popup.style.top = placement.top === null ? "auto" :
      `${placement.top}px`;
    this.popup.style.bottom = placement.bottom === null ? "auto" :
      `${placement.bottom}px`;
  }
}

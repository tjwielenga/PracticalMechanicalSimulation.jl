let virtualFileNumber = 0;
const PREVIEW_SERVICE = "http://127.0.0.1:8123/api/preview";

async function readModelPreview(file) {
  let response;
  try {
    response = await fetch(PREVIEW_SERVICE, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ name: file.name, source: await file.text() }),
    });
  } catch {
    throw new Error(
      "The SimpView model service is not running. Start bin/simpview-server " +
      "in another terminal, then reopen the model.",
    );
  }
  const value = await response.json().catch(() => ({}));
  if (!response.ok) {
    throw new Error(value.error ?? "The model preview could not be created.");
  }
  return value;
}

function scalar(dataset) {
  const value = dataset.value;
  if (ArrayBuffer.isView(value) || Array.isArray(value)) return value[0];
  return value;
}

function attribute(object, name, fallback = null) {
  try {
    return object.get_attribute(name, true);
  } catch {
    return fallback;
  }
}

function strings(dataset) {
  return Array.from(dataset.value ?? [], String);
}

function rows(dataset, columns = null) {
  const shape = dataset.shape ?? [];
  // HDF5.jl reverses array dimensions at the file boundary so C readers see
  // Julia's column-major values naturally. A Julia samples-by-components
  // matrix therefore appears here as components-by-samples.
  if (shape.length !== 2 || (columns !== null && shape[0] !== columns)) {
    throw new Error(`Unexpected graphics array shape for ${dataset.path}.`);
  }
  const [columnCount, rowCount] = shape;
  const values = dataset.value;
  return Array.from({ length: rowCount }, (_, row) =>
    Array.from({ length: columnCount }, (_, column) =>
      Number(values[column * rowCount + row])));
}

export function resampleSeries(sourceTimes, sourceValues, targetTimes) {
  if (sourceTimes.length !== sourceValues.length || sourceTimes.length === 0) {
    throw new Error("Cannot resample a signal without matching source times.");
  }
  let lower = 0;
  return targetTimes.map((target) => {
    if (target <= sourceTimes[0]) return sourceValues[0];
    if (target >= sourceTimes.at(-1)) return sourceValues.at(-1);
    while (lower + 1 < sourceTimes.length &&
        sourceTimes[lower + 1] < target) lower += 1;
    const upper = lower + 1;
    const duration = sourceTimes[upper] - sourceTimes[lower];
    if (!(duration > 0)) return sourceValues[upper];
    const fraction = (target - sourceTimes[lower]) / duration;
    return sourceValues[lower] +
      fraction * (sourceValues[upper] - sourceValues[lower]);
  });
}

function sampledVertices(dataset) {
  const shape = dataset.shape ?? [];
  if (shape.length !== 3 || shape[0] !== 3) {
    throw new Error(`Unexpected surface history shape for ${dataset.path}.`);
  }
  const [, vertexCount, sampleCount] = shape;
  const values = dataset.value;
  return Array.from({ length: sampleCount }, (_, sample) =>
    Array.from({ length: vertexCount }, (_, vertex) =>
      Array.from({ length: 3 }, (_, coordinate) =>
        Number(values[sample + sampleCount *
          (vertex + vertexCount * coordinate)]))));
}

function readMeshGroup(group) {
  const mesh = { kind: String(scalar(group.get("kind"))) };
  if (group.keys().includes("vertices")) mesh.vertices = rows(group.get("vertices"), 3);
  if (group.keys().includes("triangles")) mesh.faces = rows(group.get("triangles"), 3);
  if (group.keys().includes("edges")) mesh.edges = rows(group.get("edges"), 2);
  for (const name of ["radius_a", "radius_b", "segments"]) {
    if (group.keys().includes(name)) mesh[name] = Number(scalar(group.get(name)));
  }
  return mesh;
}

function readScene(choice) {
  const meshGroup = choice.get("meshes");
  const meshIds = strings(meshGroup.get("id"));
  const meshes = Object.fromEntries(meshIds.map((id, index) =>
    [id, readMeshGroup(meshGroup.get(String(index + 1)))]));

  const trackGroup = choice.get("tracks");
  const trackIds = strings(trackGroup.get("id"));
  const packedFields = Object.fromEntries(
    [["position", 3], ["quaternion", 4], ["scale", 3]].map(
      ([name, columns]) => [name, {
        rows: rows(trackGroup.get(name), columns),
        offsets: Array.from(trackGroup.get(`${name}_offset`).value, Number),
        counts: Array.from(trackGroup.get(`${name}_count`).value, Number),
      }]),
  );
  const tracks = trackIds.map((id, index) => ({
    id,
    ...Object.fromEntries(Object.entries(packedFields).map(([name, field]) =>
      [name, field.rows.slice(field.offsets[index],
        field.offsets[index] + field.counts[index])])),
  }));

  const instanceGroup = choice.get("instances");
  const names = strings(instanceGroup.get("name"));
  const meshNames = strings(instanceGroup.get("mesh"));
  const trackNames = strings(instanceGroup.get("track"));
  const categories = strings(instanceGroup.get("category"));
  const colors = strings(instanceGroup.get("color"));
  const opacities = Array.from(instanceGroup.get("opacity").value, Number);
  const included = Array.from(instanceGroup.get("include_in_fit").value,
    Boolean);
  const positions = rows(instanceGroup.get("local_position"), 3);
  const quaternions = rows(instanceGroup.get("local_quaternion"), 4);
  const scales = rows(instanceGroup.get("local_scale"), 3);
  const paths = strings(instanceGroup.get("path"));
  const instances = names.map((name, index) => ({
    name,
    mesh: meshNames[index],
    track: trackNames[index],
    category: categories[index],
    color: colors[index],
    opacity: opacities[index],
    include_in_fit: included[index],
    local_position: positions[index],
    local_quaternion: quaternions[index],
    local_scale: scales[index],
    path: paths[index] ? paths[index].split("\n") : [],
  }));

  const followGroup = choice.get("follow_targets");
  const followNames = strings(followGroup.get("name"));
  const followTracks = strings(followGroup.get("track"));
  const surfaceGroup = choice.get("dynamic_surfaces");
  const surfaceCount = Number(scalar(surfaceGroup.get("count")));
  const dynamicSurfaces = Array.from({ length: surfaceCount }, (_, index) => {
    const item = surfaceGroup.get(String(index + 1));
    const patches = item.get("patches");
    const patchCount = Number(scalar(patches.get("count")));
    return {
      name: String(scalar(item.get("name"))),
      category: String(scalar(item.get("category"))),
      edge_color: String(scalar(item.get("edge_color"))),
      edge_width: Number(scalar(item.get("edge_width"))),
      include_in_fit: Boolean(scalar(item.get("include_in_fit"))),
      vertices: sampledVertices(item.get("vertices")),
      edges: rows(item.get("edges"), 2),
      patches: Array.from({ length: patchCount }, (_, patchIndex) => {
        const patch = patches.get(String(patchIndex + 1));
        return {
          name: String(scalar(patch.get("name"))),
          color: String(scalar(patch.get("color"))),
          opacity: Number(scalar(patch.get("opacity"))),
          faces: rows(patch.get("triangles"), 3),
        };
      }),
    };
  });
  return {
    encoding: "mesh_instances",
    meshes,
    tracks,
    instances,
    follow_targets: followNames.map((name, index) =>
      ({ name, track: followTracks[index] })),
    dynamic_surfaces: dynamicSurfaces,
  };
}

function readChoice(choice, label) {
  const signals = choice.get("signals");
  const signalNames = strings(signals.get("name"));
  const signalValues = signalNames.length > 0
    ? rows(signals.get("values")) : [];
  return {
    label,
    times: Array.from(choice.get("time").value, Number),
    scene: readScene(choice),
    signals: signalNames.map((name, column) => ({
      name,
      values: signalValues.map((row) => row[column]),
    })),
    bookmarks: [],
  };
}

function readAppearance(group) {
  return {
    background: String(scalar(group.get("background"))),
    body_palette: strings(group.get("body_palette")),
    reaction_color: String(scalar(group.get("reaction_color"))),
    applied_color: String(scalar(group.get("applied_color"))),
    show_reactions: Boolean(scalar(group.get("show_reactions"))),
    show_applied_loads: Boolean(scalar(group.get("show_applied_loads"))),
    show_torques: Boolean(scalar(group.get("show_torques"))),
    show_ground_loads: Boolean(scalar(group.get("show_ground_loads"))),
  };
}

async function readNativeSimp(file) {
  const { default: h5wasm } = await import("h5wasm");
  const { FS } = await h5wasm.ready;
  virtualFileNumber += 1;
  const virtualPath = `/simpview-${virtualFileNumber}.simp`;
  FS.writeFile(virtualPath, new Uint8Array(await file.arrayBuffer()));
  let stored = null;
  try {
    stored = new h5wasm.File(virtualPath, "r");
    if (!stored.keys().includes("graphics")) {
      throw new Error(
        "This .simp file does not contain native graphics data.",
      );
    }
    const graphics = stored.get("graphics");
    if (String(attribute(graphics, "format", "")) !== "SimpGraphics") {
      throw new Error("The .simp graphics section has an unsupported format.");
    }
    const choiceGroup = graphics.get("choices");
    const labels = strings(choiceGroup.get("label"));
    const status = attribute(stored, "status", "complete");
    const statusMessage = attribute(stored, "status_message", "");
    const failureTime = attribute(stored, "failure_time", null);
    const statusText = String(status);
    const baseTitle = String(scalar(graphics.get("title")));
    const failureLocation = failureTime === null ? "" :
      ` at t=${Number(failureTime).toPrecision(6)} s`;
    const documentValue = {
      format: "SimpView",
      version: 4,
      title: statusText === "complete" ? baseTitle :
        `${baseTitle} — ${statusText}${failureLocation}` +
        `${statusMessage ? `: ${statusMessage}` : ""}`,
      status: statusText,
      status_message: String(statusMessage),
      failure_time: failureTime === null ? null : Number(failureTime),
      dimension: String(scalar(graphics.get("dimension"))),
      choice_name: String(scalar(graphics.get("choice_name"))),
      appearance: readAppearance(graphics.get("appearance")),
      choices: labels.map((label, index) =>
        readChoice(choiceGroup.get(String(index + 1)), label)),
    };
    const analysisMode = stored.keys().includes("diagnostics")
      ? String(attribute(stored.get("diagnostics"), "analysis_mode", "")) : "";
    for (const choice of documentValue.choices) {
      if (analysisMode === "static" || choice.label.startsWith("Static")) {
        choice.time_label = "static history sample";
      }
    }
    if (documentValue.choices.length === 1 &&
        (documentValue.choices[0].signals?.length ?? 0) === 0 &&
        stored.keys().includes("simulation") &&
        stored.keys().includes("variables") &&
        stored.keys().includes("results")) {
      const times = Array.from(stored.get("simulation/time").value);
      const components = Array.from(stored.get("variables/component").value);
      const names = Array.from(stored.get("variables/name").value);
      const values = stored.get("results/values").value;
      const sampleCount = times.length;
      const displayTimes = documentValue.choices[0].times;
      documentValue.choices[0].signals = names.map((name, variable) => {
        const angular = name === "theta" || name.startsWith("theta_") ||
          name.startsWith("psi_") || name.toLowerCase().includes("angle");
        const factor = angular ? 180 / Math.PI : 1;
        const storedValues = Array.from(values.slice(variable * sampleCount,
          (variable + 1) * sampleCount), (value) => factor * value);
        return {
          name: `${components[variable]}.${name}${angular ? " (deg)" : ""}`,
          values: resampleSeries(times, storedValues, displayTimes),
        };
      });
    }
    return documentValue;
  } finally {
    stored?.close();
    FS.unlink(virtualPath);
  }
}

export async function readViewerFile(file) {
  const lowerName = file.name.toLowerCase();
  if (lowerName.endsWith(".simp")) {
    return readNativeSimp(file);
  }
  if (lowerName.endsWith(".toml") || lowerName.endsWith(".lua")) {
    return readModelPreview(file);
  }
  throw new Error("SimpView opens .simp, .toml, and .lua files.");
}

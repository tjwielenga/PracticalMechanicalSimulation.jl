import "./style.css";
import { interpolatedSample, MechanismScene } from "./viewer.js";
import { SignalPlot } from "./plot.js";
import { SignalBrowser } from "./signal-browser.js";
import { GraphicsControl } from "./graphics-control.js";
import { prepareDocument } from "./schema.js";
import { readViewerFile } from "./file-loader.js";
import { saveResultFile } from "./file-save.js";

const elements = {
  openButton: document.querySelector("#open-button"),
  fileInput: document.querySelector("#file-input"),
  title: document.querySelector("#document-title"),
  fileName: document.querySelector("#file-name"),
  choiceControl: document.querySelector("#choice-control"),
  choiceName: document.querySelector("#choice-name"),
  choiceSelect: document.querySelector("#choice-select"),
  fitButton: document.querySelector("#fit-button"),
  runToolbar: document.querySelector("#run-toolbar"),
  analysisSelect: document.querySelector("#analysis-select"),
  endTimeInput: document.querySelector("#end-time-input"),
  frameRateInput: document.querySelector("#frame-rate-input"),
  runButton: document.querySelector("#run-button"),
  downloadButton: document.querySelector("#download-button"),
  runStatus: document.querySelector("#run-status"),
  dropMessage: document.querySelector("#drop-message"),
  xValueReadout: document.querySelector("#x-value-readout"),
  yValueReadout: document.querySelector("#y-value-readout"),
  sampleReadout: document.querySelector("#sample-readout"),
  playButton: document.querySelector("#play-button"),
  modeScaleControl: document.querySelector("#mode-scale-control"),
  modeScaleInput: document.querySelector("#mode-scale-input"),
  modeScaleReadout: document.querySelector("#mode-scale-readout"),
  modeScaleReset: document.querySelector("#mode-scale-reset"),
  timeSlider: document.querySelector("#time-slider"),
  timeReadout: document.querySelector("#time-readout"),
};

const mechanism = new MechanismScene(document.querySelector("#scene-canvas"));
const graphicsControl = new GraphicsControl(
  document.querySelector("#graphics-control"), {
    getScale: (path) => mechanism.getGraphicScale(path),
    getEffectiveScale: (path) => mechanism.getEffectiveGraphicScale(path),
    setScale: (path, scale) => mechanism.setGraphicScale(path, scale),
    resetScale: (path) => mechanism.resetGraphicScale(path),
    resetAllScales: () => mechanism.resetAllGraphicScales(),
    getDeformationScale: (group) => mechanism.getDeformationScale(group),
    setDeformationScale: (group, scale) =>
      mechanism.setDeformationScale(group, scale),
    resetDeformationScale: (group) =>
      mechanism.resetDeformationScale(group),
    getVisible: (path) => mechanism.getGraphicVisible(path),
    setVisible: (path, visible) => mechanism.setGraphicVisible(path, visible),
    getFollowTarget: () => mechanism.getFollowTarget(),
    setFollowTarget: (name) => mechanism.setFollowTarget(name),
    getFollowRotation: () => mechanism.getFollowRotation(),
    setFollowRotation: (enabled) => mechanism.setFollowRotation(enabled),
  });
const plot = new SignalPlot(document.querySelector("#plot-canvas"));
const xBrowser = new SignalBrowser(document.querySelector("#x-signal-browser"), {
  includeTime: true,
  onSelect: updatePlotSignals,
});
const yBrowser = new SignalBrowser(document.querySelector("#y-signal-browser"), {
  onSelect: updatePlotSignals,
});

let viewerDocument = null;
let choice = null;
let sampleIndex = 0;
let playing = false;
let animationStart = 0;
let resultStart = 0;
let openedModel = null;
let activeRun = null;
let runPollTimer = null;
let staticRun = null;
let dynamicRun = null;
let stateRun = null;
let runBusy = false;
let modeScale = 1;

const RUN_SERVICE = "http://127.0.0.1:8123/api/run";

function populateSignals() {
  const xPreference = xBrowser.selected?.id ?? "time";
  const yPreference = yBrowser.selected?.id ?? "signal:0";
  xBrowser.setSignals(choice.signals ?? [], xPreference,
    choice.time_label ?? "time (s)");
  yBrowser.setSignals(choice.signals ?? [], yPreference);
  updatePlotSignals();
}

function browserSignal(browser) {
  const index = browser.selected?.index;
  if (index === null || index === undefined) return null;
  const signal = choice?.signals?.[index];
  if (!signal || viewerDocument?.choice_name !== "Mode" || modeScale === 1) {
    return signal;
  }
  return { ...signal, values: signal.values.map((value) => value * modeScale) };
}

function formatModeScale(value) {
  if (value >= 100 || value < 0.01) return `${value.toPrecision(3)}×`;
  if (value >= 10 || value < 0.1) return `${value.toFixed(2)}×`;
  return `${value.toFixed(3).replace(/0+$/, "").replace(/\.$/, "")}×`;
}

function setModeScale(value) {
  modeScale = Number.isFinite(value) && value > 0 ? value : 1;
  elements.modeScaleInput.value = String(Math.log10(modeScale));
  elements.modeScaleReadout.value = formatModeScale(modeScale);
  mechanism.setModeScale(modeScale);
  updatePlotSignals();
  setSample(sampleIndex);
}

function updatePlotSignals() {
  plot.setSignals(choice?.times ?? [], browserSignal(xBrowser),
    browserSignal(yBrowser), choice?.time_label ?? "time (s)");
  plot.setSample(Math.round(sampleIndex));
}

function setSample(index) {
  if (!choice) return;
  sampleIndex = Math.max(0, Math.min(Number(index), choice.times.length - 1));
  elements.timeSlider.value = String(sampleIndex);
  mechanism.update(sampleIndex);
  const nearestIndex = Math.round(sampleIndex);
  plot.setSample(nearestIndex);
  const time = interpolatedSample(choice.times, sampleIndex);
  elements.timeReadout.value = choice.time_label === "static history sample"
    ? "" : choice.time_label === "mode phase"
      ? `phase ${(100 * time).toFixed(1)}%` : `${time.toFixed(3)} s`;
  const xSignal = browserSignal(xBrowser);
  const ySignal = browserSignal(yBrowser);
  const xValue = xSignal ? interpolatedSample(xSignal.values, sampleIndex) : time;
  const yValue = ySignal ? interpolatedSample(ySignal.values, sampleIndex) : NaN;
  elements.xValueReadout.value = Number.isFinite(xValue)
    ? `X ${Number(xValue).toPrecision(6)}` : "X —";
  elements.yValueReadout.value = Number.isFinite(yValue)
    ? `Y ${Number(yValue).toPrecision(6)}` : "Y —";
  elements.sampleReadout.value =
    `${nearestIndex + 1} / ${choice.times.length}`;
}

function selectChoice(index, preserveFollow = true, resetView = false) {
  const followedBody = preserveFollow ? mechanism.getFollowTarget() : null;
  choice = viewerDocument.choices[index];
  playing = false;
  elements.playButton.textContent = "Play";
  elements.timeSlider.max = String(choice.times.length - 1);
  elements.timeSlider.step = choice.time_label === "static history sample"
    ? "1" : "any";
  elements.playButton.disabled = choice.times.length < 2;
  elements.timeSlider.disabled = choice.times.length < 2;
  elements.modeScaleControl.hidden = viewerDocument.choice_name !== "Mode";
  mechanism.load(choice, viewerDocument.appearance, viewerDocument.dimension,
    resetView, viewerDocument.choice_name === "Mode");
  mechanism.setModeScale(modeScale);
  if (followedBody && mechanism.followTargetNames().includes(followedBody)) {
    mechanism.setFollowTarget(followedBody);
  }
  graphicsControl.setPaths(mechanism.graphicPaths(),
    mechanism.followTargetNames(), mechanism.deformationPaths());
  populateSignals();
  setSample(0);
}

function loadDocument(documentValue, fileName, options = {}) {
  const previousChoice = Number(elements.choiceSelect.value || 0);
  const previousSample = sampleIndex;
  const previousAtEnd = choice && sampleIndex >= choice.times.length - 1;
  viewerDocument = prepareDocument(documentValue);
  if (options.resetGraphics !== false) mechanism.resetGraphicSettings();
  elements.title.textContent = viewerDocument.title || "SimpView";
  elements.fileName.textContent = fileName;
  elements.choiceName.textContent = viewerDocument.choice_name || "Result";
  elements.choiceSelect.replaceChildren();
  for (const [index, item] of viewerDocument.choices.entries()) {
    const option = document.createElement("option");
    option.value = String(index);
    option.textContent = item.label;
    elements.choiceSelect.append(option);
  }
  elements.choiceControl.hidden = viewerDocument.choices.length < 2;
  elements.dropMessage.hidden = true;
  elements.fitButton.disabled = false;
  const choiceIndex = Math.min(options.choiceIndex ?? previousChoice,
    viewerDocument.choices.length - 1);
  selectChoice(choiceIndex, options.resetGraphics === false,
    options.resetView !== false);
  if (options.followLatest || (previousAtEnd && !options.resetSample)) {
    setSample(choice.times.length - 1);
  } else if (options.preserveSample) {
    setSample(Math.min(previousSample, choice.times.length - 1));
  }
}

function configureRunControls(documentValue, model) {
  openedModel = model;
  staticRun = null;
  dynamicRun = null;
  stateRun = null;
  elements.runToolbar.hidden = !model;
  elements.downloadButton.hidden = true;
  elements.runStatus.textContent = model ? "Ready" : "";
  if (!model) return;
  const settings = documentValue.run_settings ?? {};
  const requested = settings.analysis_mode === "static" ? "static"
    : settings.analysis_mode === "modal" ? "modal" : "dynamic";
  elements.analysisSelect.value = requested;
  elements.endTimeInput.value = String(settings.end_time ?? 1);
  elements.frameRateInput.value = String(
    Math.max(1, Math.min(240, Math.round(settings.frames_per_second ?? 60))));
  updateAnalysisControls();
}

function updateAnalysisControls() {
  const analysis = elements.analysisSelect.value;
  const timed = analysis === "dynamic";
  elements.endTimeInput.disabled = !timed || runBusy;
  elements.frameRateInput.disabled = !timed || runBusy;
  elements.runButton.textContent = analysis === "initial_conditions"
    ? "Calculate ICs" : analysis === "static" ? "Run static"
      : analysis === "modal" ? "Calculate modes" : "Run dynamic";
}

async function openFile(file) {
  if (!file) return;
  try {
    clearTimeout(runPollTimer);
    activeRun = null;
    setModeScale(1);
    configureRunControls({}, null);
    elements.fitButton.disabled = true;
    elements.fileName.textContent = `Opening ${file.name}…`;
    const documentValue = await readViewerFile(file);
    loadDocument(documentValue, file.name);
    const lowerName = file.name.toLowerCase();
    const isModel = lowerName.endsWith(".toml") || lowerName.endsWith(".lua");
    configureRunControls(documentValue, isModel ? {
      name: file.name,
      source: await file.text(),
    } : null);
  } catch (error) {
    console.error(error);
    elements.fileName.textContent = error.message;
  }
}

function setRunBusy(busy) {
  runBusy = busy;
  elements.openButton.disabled = busy;
  elements.runButton.disabled = busy;
  elements.analysisSelect.disabled = busy;
  updateAnalysisControls();
}

function scheduleRunPoll(delay = 350) {
  clearTimeout(runPollTimer);
  runPollTimer = setTimeout(pollRun, delay);
}

function runStatusLabel(snapshot) {
  const count = snapshot.document?.choices?.[0]?.times?.length ?? 0;
  const label = snapshot.analysis_mode === "initial_conditions"
    ? "Initial conditions" : snapshot.analysis_mode === "static"
      ? "Static equilibrium" : snapshot.analysis_mode === "modal"
        ? "Modal analysis" : "Dynamic simulation";
  if (snapshot.status === "complete") return `Complete — ${count} frames`;
  if (snapshot.status === "failed") return `Failed — ${snapshot.message}`;
  if (snapshot.status === "interrupted") return "Interrupted";
  return `${label} — ${count} frames available`;
}

async function pollRun() {
  if (!activeRun) return;
  try {
    const response = await fetch(
      `${RUN_SERVICE}/${activeRun.id}?after=${activeRun.revision}`,
      { cache: "no-store" },
    );
    if (response.status === 204) {
      scheduleRunPoll();
      return;
    }
    const snapshot = await response.json().catch(() => ({}));
    if (!response.ok) {
      throw new Error(snapshot.error ?? "The running analysis could not be read.");
    }
    activeRun.revision = snapshot.revision;
    activeRun.resultName = snapshot.result_name;
    const modalResult = snapshot.status === "complete" &&
      snapshot.analysis_mode === "modal";
    if (modalResult) setModeScale(1);
    loadDocument(snapshot.document, `${openedModel.name} — live`, {
      resetGraphics: false,
      resetView: false,
      choiceIndex: modalResult ? 0 : undefined,
      preserveSample: !modalResult,
      resetSample: modalResult,
    });
    elements.runStatus.textContent = runStatusLabel(snapshot);
    const finished = ["complete", "failed", "interrupted"]
      .includes(snapshot.status);
    if (finished) {
      elements.downloadButton.hidden = !snapshot.result_ready;
      setRunBusy(false);
      activeRun.finished = true;
      if (snapshot.status === "complete" &&
          activeRun.analysis === "static") {
        staticRun = { id: activeRun.id };
        dynamicRun = null;
        stateRun = staticRun;
      } else if (snapshot.status === "complete" &&
          activeRun.analysis === "dynamic") {
        dynamicRun = { id: activeRun.id };
        stateRun = dynamicRun;
      } else if (snapshot.status === "complete" &&
          activeRun.analysis === "initial_conditions") {
        staticRun = null;
        dynamicRun = null;
        stateRun = null;
      }
      updateAnalysisControls();
      return;
    }
    scheduleRunPoll();
  } catch (error) {
    console.error(error);
    elements.runStatus.textContent = error.message;
    setRunBusy(false);
  }
}

async function startRun() {
  if (!openedModel) return;
  const endTime = Number(elements.endTimeInput.value);
  const framesPerSecond = Number(elements.frameRateInput.value);
  const analysis = elements.analysisSelect.value;
  if (!Number.isFinite(endTime) || !Number.isFinite(framesPerSecond) ||
      framesPerSecond <= 0) {
    elements.runStatus.textContent =
      "End time and frames/s must be valid numbers.";
    return;
  }
  setRunBusy(true);
  elements.downloadButton.hidden = true;
  elements.runStatus.textContent = "Starting analysis…";
  try {
    const response = await fetch(RUN_SERVICE, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        name: openedModel.name,
        source: openedModel.source,
        end_time: endTime,
        frames_per_second: framesPerSecond,
        analysis,
        initial_run_id: analysis === "dynamic" || analysis === "modal"
          ? stateRun?.id : undefined,
      }),
    });
    const value = await response.json().catch(() => ({}));
    if (!response.ok) {
      throw new Error(value.error ?? "The analysis could not be started.");
    }
    activeRun = {
      id: value.id,
      revision: -1,
      resultName: value.result_name,
      analysis,
    };
    pollRun();
  } catch (error) {
    console.error(error);
    elements.runStatus.textContent = error.message;
    setRunBusy(false);
  }
}

async function downloadResult() {
  if (!activeRun?.finished) return;
  elements.downloadButton.disabled = true;
  try {
    const result = await saveResultFile(
      `${RUN_SERVICE}/${activeRun.id}/result`, activeRun.resultName);
    if (result === "saved") elements.runStatus.textContent = "Result saved";
  } catch (error) {
    console.error(error);
    elements.runStatus.textContent = error.message;
  } finally {
    elements.downloadButton.disabled = false;
  }
}

function findSampleAtTime(time) {
  let low = 0;
  let high = choice.times.length - 1;
  while (low < high) {
    const middle = Math.floor((low + high) / 2);
    if (choice.times[middle] < time) low = middle + 1;
    else high = middle;
  }
  if (low === 0 || choice.times[low] === time) return low;
  const lower = low - 1;
  const interval = choice.times[low] - choice.times[lower];
  return interval > 0 ? lower + (time - choice.times[lower]) / interval : low;
}

function animationFrame(now) {
  if (!playing || !choice) return;
  const firstTime = choice.times[0];
  const lastTime = choice.times.at(-1);
  const duration = lastTime - firstTime;
  const playbackRate = choice.time_label === "mode phase" ? 0.5 : 1;
  const elapsed = resultStart - firstTime +
    playbackRate * (now - animationStart) / 1000;
  const resultTime = duration > 0
    ? firstTime + (elapsed % duration)
    : firstTime;
  setSample(findSampleAtTime(resultTime));
  requestAnimationFrame(animationFrame);
}

elements.openButton.addEventListener("click", () => elements.fileInput.click());
elements.fileInput.addEventListener("change", () =>
  openFile(elements.fileInput.files[0]));
elements.choiceSelect.addEventListener("change", () =>
  selectChoice(Number(elements.choiceSelect.value)));
elements.fitButton.addEventListener("click", () => mechanism.fit());
elements.analysisSelect.addEventListener("change", updateAnalysisControls);
elements.runButton.addEventListener("click", startRun);
elements.downloadButton.addEventListener("click", downloadResult);
elements.timeSlider.addEventListener("input", () => {
  playing = false;
  elements.playButton.textContent = "Play";
  setSample(Number(elements.timeSlider.value));
});
elements.modeScaleInput.addEventListener("input", () =>
  setModeScale(10 ** Number(elements.modeScaleInput.value)));
elements.modeScaleInput.addEventListener("dblclick", () => setModeScale(1));
elements.modeScaleReset.addEventListener("click", () => setModeScale(1));
elements.playButton.addEventListener("click", () => {
  if (!choice) return;
  if (playing) {
    playing = false;
    elements.playButton.textContent = "Play";
    return;
  }
  if (sampleIndex >= choice.times.length - 1) setSample(0);
  playing = true;
  elements.playButton.textContent = "Pause";
  animationStart = performance.now();
  resultStart = interpolatedSample(choice.times, sampleIndex);
  requestAnimationFrame(animationFrame);
});

window.addEventListener("dragover", (event) => {
  event.preventDefault();
  document.body.classList.add("dragging");
});
window.addEventListener("dragleave", (event) => {
  if (event.relatedTarget === null) document.body.classList.remove("dragging");
});
window.addEventListener("drop", (event) => {
  event.preventDefault();
  document.body.classList.remove("dragging");
  openFile(event.dataTransfer.files[0]);
});

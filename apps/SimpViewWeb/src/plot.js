function finiteExtent(values) {
  let low = Infinity;
  let high = -Infinity;
  for (const value of values) {
    if (!Number.isFinite(value)) continue;
    low = Math.min(low, value);
    high = Math.max(high, value);
  }
  return [low, high];
}

function niceStep(value) {
  const exponent = Math.floor(Math.log10(value));
  const power = 10 ** exponent;
  const fraction = value / power;
  const niceFraction = fraction < 1.5 ? 1 : fraction < 3 ? 2 :
    fraction < 7 ? 5 : 10;
  return niceFraction * power;
}

export function niceAxis(values, targetDivisions = 5) {
  let [low, high] = finiteExtent(values);
  if (!Number.isFinite(low)) [low, high] = [-1, 1];
  if (low === high) {
    const halfRange = low === 0 ? 1 : Math.max(Math.abs(low) * 0.05, 1e-12);
    low -= halfRange;
    high += halfRange;
  }
  const step = niceStep((high - low) / targetDivisions);
  const lowerIndex = Math.floor(low / step + 1e-12);
  const upperIndex = Math.ceil(high / step - 1e-12);
  const ticks = Array.from({ length: upperIndex - lowerIndex + 1 },
    (_, offset) => {
      const value = (lowerIndex + offset) * step;
      if (Math.abs(value) < step * 1e-12) return 0;
      return Number(value.toPrecision(14));
    });
  return { low: ticks[0], high: ticks.at(-1), step, ticks };
}

export function numberLabel(value, step) {
  if (value === 0) return "0";
  if (Math.abs(value) >= 1e5 || Math.abs(value) < 1e-4) {
    return value.toExponential(2).replace(/\.?(?:0+)(?=e)/, "");
  }
  const decimalPlaces = Math.max(0, -Math.floor(Math.log10(step)));
  return value.toFixed(Math.min(decimalPlaces, 12));
}

export class SignalPlot {
  constructor(canvas) {
    this.canvas = canvas;
    this.context = canvas.getContext("2d");
    this.xValues = [];
    this.yValues = [];
    this.xLabel = "";
    this.yLabel = "";
    this.sample = 0;
    this.resizeObserver = new ResizeObserver(() => this.draw());
    this.resizeObserver.observe(canvas);
  }

  setSignals(times, xSignal, ySignal, timeLabel = "time (s)") {
    this.xValues = xSignal?.values ?? times;
    this.yValues = ySignal?.values ?? [];
    this.xLabel = xSignal?.name ?? timeLabel;
    this.yLabel = ySignal?.name ?? "";
    this.sample = 0;
    this.draw();
  }

  setSample(sample) {
    this.sample = sample;
    this.draw();
  }

  draw() {
    const bounds = this.canvas.getBoundingClientRect();
    if (bounds.width < 2 || bounds.height < 2) return;
    const ratio = window.devicePixelRatio || 1;
    const width = Math.round(bounds.width * ratio);
    const height = Math.round(bounds.height * ratio);
    if (this.canvas.width !== width || this.canvas.height !== height) {
      this.canvas.width = width;
      this.canvas.height = height;
    }
    const context = this.context;
    context.setTransform(ratio, 0, 0, ratio, 0, 0);
    context.clearRect(0, 0, bounds.width, bounds.height);
    context.fillStyle = "#fff";
    context.fillRect(0, 0, bounds.width, bounds.height);
    if (this.xValues.length < 2 ||
        this.yValues.length !== this.xValues.length) {
      context.fillStyle = "#7a848f";
      context.font = "13px system-ui";
      context.fillText("Select a signal to plot", 18, 28);
      return;
    }

    const margin = { left: 78, right: 22, top: 16, bottom: 42 };
    const plotWidth = Math.max(bounds.width - margin.left - margin.right, 1);
    const plotHeight = Math.max(bounds.height - margin.top - margin.bottom, 1);
    const xAxis = niceAxis(this.xValues);
    const yAxis = niceAxis(this.yValues);
    const { low: xLow, high: xHigh } = xAxis;
    const { low: yLow, high: yHigh } = yAxis;
    const x = (value) => margin.left +
      ((value - xLow) / Math.max(xHigh - xLow, Number.EPSILON)) * plotWidth;
    const y = (value) => margin.top +
      (1 - (value - yLow) / (yHigh - yLow)) * plotHeight;

    context.strokeStyle = "#d6dce2";
    context.lineWidth = 1;
    context.fillStyle = "#64707c";
    context.font = "11px ui-monospace, monospace";
    context.textBaseline = "middle";
    context.textAlign = "right";
    for (const value of yAxis.ticks) {
      const yy = y(value);
      context.beginPath();
      context.moveTo(margin.left, yy);
      context.lineTo(margin.left + plotWidth, yy);
      context.stroke();
      context.fillText(numberLabel(value, yAxis.step), margin.left - 8, yy);
    }

    context.textBaseline = "top";
    const xTicks = xAxis.ticks;
    for (const [index, value] of xTicks.entries()) {
      const xx = x(value);
      context.beginPath();
      context.moveTo(xx, margin.top);
      context.lineTo(xx, margin.top + plotHeight);
      context.stroke();
      context.textAlign = index === 0 ? "left" :
        index === xTicks.length - 1 ? "right" : "center";
      context.fillText(numberLabel(value, xAxis.step), xx,
        margin.top + plotHeight + 5);
    }

    context.strokeStyle = "#27678f";
    context.lineWidth = 1.8;
    context.beginPath();
    let drawing = false;
    for (let index = 0; index < this.xValues.length; index += 1) {
      const xValue = this.xValues[index];
      const yValue = this.yValues[index];
      if (!Number.isFinite(xValue) || !Number.isFinite(yValue)) {
        drawing = false;
        continue;
      }
      const xx = x(xValue);
      const yy = y(yValue);
      if (!drawing) context.moveTo(xx, yy);
      else context.lineTo(xx, yy);
      drawing = true;
    }
    context.stroke();

    const selected = Math.max(0,
      Math.min(this.sample, this.xValues.length - 1));
    if (Number.isFinite(this.xValues[selected]) &&
        Number.isFinite(this.yValues[selected])) {
      context.fillStyle = "#c24b37";
      context.beginPath();
      context.arc(x(this.xValues[selected]), y(this.yValues[selected]),
        3.8, 0, 2 * Math.PI);
      context.fill();
    }

    context.fillStyle = "#64707c";
    context.font = "11px system-ui";
    context.textAlign = "center";
    context.textBaseline = "bottom";
    context.fillText(this.xLabel, margin.left + plotWidth / 2,
      bounds.height - 3);
    context.save();
    context.translate(13, margin.top + plotHeight / 2);
    context.rotate(-Math.PI / 2);
    context.fillText(this.yLabel, 0, 0);
    context.restore();
  }
}

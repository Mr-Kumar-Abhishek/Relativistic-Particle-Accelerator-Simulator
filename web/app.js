(() => {
  "use strict";

  const PARTICLE_STRIDE = 8;
  const PHASE_STRIDE = 6;
  const DIAG_LEN = 20;

  let wasm = null;
  let running = true;
  let trailHistory = [];
  let telemetryHistory = [];

  const $ = (id) => document.getElementById(id);

  const ctrls = {
    mode: $("ctrl-mode"),
    species: $("ctrl-species"),
    energy: $("ctrl-energy"),
    b: $("ctrl-b"),
    quad: $("ctrl-quad"),
    rf: $("ctrl-rf"),
    count: $("ctrl-count"),
    substeps: $("ctrl-substeps"),
  };

  function syncLabels() {
    $("val-energy").textContent = ctrls.energy.value;
    $("val-b").textContent = Number(ctrls.b.value).toFixed(2);
    $("val-quad").textContent = Number(ctrls.quad.value).toFixed(2);
    $("val-rf").textContent = ctrls.rf.value;
    $("val-count").textContent = ctrls.count.value;
    $("val-substeps").textContent = ctrls.substeps.value;
  }

  function resetSimulation() {
    syncLabels();
    if (!wasm) return;
    wasm.sim_init(
      Number(ctrls.mode.value),
      Number(ctrls.species.value),
      Number(ctrls.count.value),
      Number(ctrls.energy.value),
      Number(ctrls.b.value),
      Number(ctrls.quad.value),
      Number(ctrls.rf.value)
    );
    trailHistory = [];
    telemetryHistory = [];
    renderFrame();
  }

  function resizeCanvas(canvas) {
    const dpr = window.devicePixelRatio || 1;
    const rect = canvas.getBoundingClientRect();
    const w = Math.max(200, Math.floor(rect.width * dpr));
    const h = Math.max(180, Math.floor(rect.height * dpr));
    if (canvas.width !== w || canvas.height !== h) {
      canvas.width = w;
      canvas.height = h;
    }
    return { ctx: canvas.getContext("2d"), w, h, dpr };
  }

  function formatSci(v, digits = 3) {
    if (!Number.isFinite(v)) return "--";
    if (v === 0) return "0.000";
    if (Math.abs(v) < 1e-2 || Math.abs(v) >= 1e4) return v.toExponential(digits - 1);
    return v.toFixed(digits);
  }

  function renderOrbit(particles, count, ringRadius, mode) {
    const { ctx, w, h, dpr } = resizeCanvas($("cv-orbit"));
    ctx.clearRect(0, 0, w, h);

    const cx = w * 0.5;
    const cy = h * 0.5;
    const scale = (Math.min(w, h) * 0.40) / Math.max(ringRadius, 0.1);

    // Subtle polar grid
    ctx.strokeStyle = "#162238";
    ctx.lineWidth = 1;
    for (let k = 1; k <= 4; k++) {
      ctx.beginPath();
      ctx.arc(cx, cy, (Math.min(w, h) * 0.40 * k) / 4, 0, Math.PI * 2);
      ctx.stroke();
    }

    if (mode === 0) {
      // Synchrotron FODO Ring & 8 FODO magnet sectors
      const cells = 8;
      for (let c = 0; c < cells; c++) {
        const a0 = (c / cells) * Math.PI * 2;
        const da = (Math.PI * 2) / cells;
        // Dipole bending arcs
        ctx.strokeStyle = "rgba(56, 189, 248, 0.25)";
        ctx.lineWidth = 10 * dpr;
        ctx.beginPath();
        ctx.arc(cx, cy, ringRadius * scale, a0 + da * 0.12, a0 + da * 0.40);
        ctx.stroke();
        ctx.beginPath();
        ctx.arc(cx, cy, ringRadius * scale, a0 + da * 0.62, a0 + da * 0.90);
        ctx.stroke();

        // QF & QD quadrupole markers
        ctx.strokeStyle = "rgba(245, 158, 11, 0.7)";
        ctx.lineWidth = 12 * dpr;
        ctx.beginPath();
        ctx.arc(cx, cy, ringRadius * scale, a0, a0 + da * 0.08);
        ctx.stroke();

        ctx.strokeStyle = "rgba(168, 85, 247, 0.7)";
        ctx.beginPath();
        ctx.arc(cx, cy, ringRadius * scale, a0 + da * 0.50, a0 + da * 0.58);
        ctx.stroke();
      }
      // RF Cavity indicator
      ctx.fillStyle = "#22c55e";
      const rfA = (0.45 / cells) * Math.PI * 2;
      ctx.fillRect(
        cx + Math.cos(rfA) * ringRadius * scale - 5 * dpr,
        cy + Math.sin(rfA) * ringRadius * scale - 5 * dpr,
        10 * dpr,
        10 * dpr
      );
    } else {
      // Cyclotron Dees + Extraction radius
      const R = ringRadius * scale;
      ctx.strokeStyle = "rgba(56, 189, 248, 0.35)";
      ctx.lineWidth = 2 * dpr;
      ctx.beginPath();
      ctx.arc(cx, cy, R, Math.PI * 0.52, Math.PI * 1.48);
      ctx.stroke();
      ctx.beginPath();
      ctx.arc(cx, cy, R, -Math.PI * 0.48, Math.PI * 0.48);
      ctx.stroke();

      ctx.setLineDash([4 * dpr, 4 * dpr]);
      ctx.strokeStyle = "rgba(245, 158, 11, 0.4)";
      ctx.beginPath();
      ctx.moveTo(cx, cy - R);
      ctx.lineTo(cx, cy + R);
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // Draw reference particle trail
    if (count > 0) {
      trailHistory.push([particles[0], particles[1]]);
      if (trailHistory.length > 350) trailHistory.shift();
      ctx.strokeStyle = "rgba(56, 189, 248, 0.45)";
      ctx.lineWidth = 1.6 * dpr;
      ctx.beginPath();
      for (let i = 0; i < trailHistory.length; i++) {
        const px = cx + trailHistory[i][0] * scale;
        const py = cy - trailHistory[i][1] * scale;
        if (i === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      }
      ctx.stroke();
    }

    // Draw particles colored by kinetic energy spread
    for (let i = 0; i < count; i++) {
      const base = i * PARTICLE_STRIDE;
      const px = cx + particles[base] * scale;
      const py = cy - particles[base + 1] * scale;
      ctx.fillStyle = i === 0 ? "#f8fafc" : "#38bdf8";
      ctx.beginPath();
      ctx.arc(px, py, (i === 0 ? 3.5 : 2.1) * dpr, 0, Math.PI * 2);
      ctx.fill();
    }
  }

  function renderPhaseSpace(phase, count, diag) {
    const { ctx, w, h, dpr } = resizeCanvas($("cv-phase"));
    ctx.clearRect(0, 0, w, h);

    const cx = w * 0.5;
    const cy = h * 0.5;

    ctx.strokeStyle = "#1b2942";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, cy);
    ctx.lineTo(w, cy);
    ctx.moveTo(cx, 0);
    ctx.lineTo(cx, h);
    ctx.stroke();

    let maxX = 1e-4;
    let maxXp = 1e-4;
    for (let i = 0; i < count; i++) {
      const b = i * PHASE_STRIDE;
      maxX = Math.max(maxX, Math.abs(phase[b]) * 1000.0);
      maxXp = Math.max(maxXp, Math.abs(phase[b + 1]) * 1000.0);
    }
    const scaleX = (w * 0.42) / (maxX * 1.15);
    const scaleXp = (h * 0.42) / (maxXp * 1.15);

    // Draw Twiss RMS Ellipse
    const alpha = diag[10];
    const beta = diag[11];
    const emit = diag[8] * 1e-3; // mm*mrad -> mm*rad? Actually diag[8] is mm*mrad -> m*rad * 1e6.
    // Convert x (mm) = 1000 * sqrt(eps_m_rad * beta), xp (mrad) = -1000 * sqrt(eps_m_rad / beta) * (alpha*cos + sin)
    const eps_m_rad = diag[8] * 1e-6;
    if (eps_m_rad > 0 && beta > 0) {
      ctx.strokeStyle = "rgba(168, 85, 247, 0.85)";
      ctx.lineWidth = 2 * dpr;
      ctx.beginPath();
      for (let k = 0; k <= 96; k++) {
        const th = (k / 96) * Math.PI * 2;
        const xMm = 1000.0 * Math.sqrt(eps_m_rad * beta) * Math.cos(th);
        const xpMrad =
          -1000.0 * Math.sqrt(eps_m_rad / beta) * (alpha * Math.cos(th) + Math.sin(th));
        const sx = cx + xMm * scaleX;
        const sy = cy - xpMrad * scaleXp;
        if (k === 0) ctx.moveTo(sx, sy);
        else ctx.lineTo(sx, sy);
      }
      ctx.stroke();
    }

    // Scatter particles
    ctx.fillStyle = "#38bdf8";
    for (let i = 0; i < count; i++) {
      const b = i * PHASE_STRIDE;
      const xMm = phase[b] * 1000.0;
      const xpMrad = phase[b + 1] * 1000.0;
      ctx.beginPath();
      ctx.arc(cx + xMm * scaleX, cy - xpMrad * scaleXp, 2.2 * dpr, 0, Math.PI * 2);
      ctx.fill();
    }
  }

  function renderRfSpace(phase, count) {
    const { ctx, w, h, dpr } = resizeCanvas($("cv-rf"));
    ctx.clearRect(0, 0, w, h);

    const cx = w * 0.5;
    const cy = h * 0.5;
    ctx.strokeStyle = "#1b2942";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, cy);
    ctx.lineTo(w, cy);
    ctx.moveTo(cx, 0);
    ctx.lineTo(cx, h);
    ctx.stroke();

    // Separatrix bucket sketch
    ctx.strokeStyle = "rgba(34, 197, 94, 0.35)";
    ctx.lineWidth = 1.5 * dpr;
    ctx.beginPath();
    for (let k = 0; k <= 120; k++) {
      const phi = -Math.PI + (k / 120) * Math.PI * 2;
      const amp = h * 0.35 * Math.cos(phi * 0.5);
      const sx = cx + (phi / Math.PI) * (w * 0.44);
      if (k === 0) ctx.moveTo(sx, cy - amp);
      else ctx.lineTo(sx, cy - amp);
    }
    for (let k = 120; k >= 0; k--) {
      const phi = -Math.PI + (k / 120) * Math.PI * 2;
      const amp = h * 0.35 * Math.cos(phi * 0.5);
      const sx = cx + (phi / Math.PI) * (w * 0.44);
      ctx.lineTo(sx, cy + amp);
    }
    ctx.stroke();

    let maxDe = 1e-3;
    for (let i = 0; i < count; i++) {
      maxDe = Math.max(maxDe, Math.abs(phase[i * PHASE_STRIDE + 5]));
    }
    const scaleDe = (h * 0.40) / (maxDe * 1.15);

    ctx.fillStyle = "#22c55e";
    for (let i = 0; i < count; i++) {
      const b = i * PHASE_STRIDE;
      const dphi = phase[b + 4];
      const de = phase[b + 5];
      const sx = cx + (dphi / Math.PI) * (w * 0.44);
      const sy = cy - de * scaleDe;
      ctx.beginPath();
      ctx.arc(sx, sy, 2.2 * dpr, 0, Math.PI * 2);
      ctx.fill();
    }
  }

  function renderHistory(diag) {
    telemetryHistory.push({
      sx: diag[6],
      sy: diag[7],
      ek: diag[2],
    });
    if (telemetryHistory.length > 240) telemetryHistory.shift();

    const { ctx, w, h, dpr } = resizeCanvas($("cv-history"));
    ctx.clearRect(0, 0, w, h);

    if (telemetryHistory.length < 2) return;
    let maxEnv = 1e-3;
    let maxEk = 1e-3;
    let minEk = Infinity;
    for (const pt of telemetryHistory) {
      maxEnv = Math.max(maxEnv, pt.sx, pt.sy);
      maxEk = Math.max(maxEk, pt.ek);
      minEk = Math.min(minEk, pt.ek);
    }
    const ekSpan = Math.max(maxEk - minEk, 1e-3);

    const drawSeries = (getter, color) => {
      ctx.strokeStyle = color;
      ctx.lineWidth = 2 * dpr;
      ctx.beginPath();
      for (let i = 0; i < telemetryHistory.length; i++) {
        const px = (i / (telemetryHistory.length - 1)) * (w - 20 * dpr) + 10 * dpr;
        const py = h - 14 * dpr - getter(telemetryHistory[i]) * (h - 28 * dpr);
        if (i === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      }
      ctx.stroke();
    };

    drawSeries((p) => p.sx / maxEnv, "#38bdf8");
    drawSeries((p) => p.sy / maxEnv, "#a855f7");
    drawSeries((p) => (p.ek - minEk) / ekSpan, "#22c55e");

    ctx.font = `${11 * dpr}px sans-serif`;
    ctx.fillStyle = "#38bdf8";
    ctx.fillText("RMS sigma_x", 12 * dpr, 16 * dpr);
    ctx.fillStyle = "#a855f7";
    ctx.fillText("RMS sigma_y", 96 * dpr, 16 * dpr);
    ctx.fillStyle = "#22c55e";
    ctx.fillText("Mean Ek (MeV)", 180 * dpr, 16 * dpr);
  }

  function renderFrame() {
    if (!wasm) return;
    const count = wasm.sim_get_particle_count();
    const buf = wasm.memory.buffer;
    const particles = new Float64Array(buf, wasm.sim_get_particles_ptr(), count * PARTICLE_STRIDE);
    const phase = new Float64Array(buf, wasm.sim_get_phase_space_ptr(), count * PHASE_STRIDE);
    const diag = new Float64Array(buf, wasm.sim_get_diagnostics_ptr(), DIAG_LEN);

    $("m-gamma").textContent = diag[0].toFixed(6);
    $("m-beta").textContent = diag[1].toFixed(6);
    $("m-ek").textContent = `${formatSci(diag[2], 4)} MeV`;
    $("m-brho").textContent = `${formatSci(diag[3], 4)} T*m`;
    $("m-rl").textContent = `${formatSci(diag[4], 4)} m`;
    $("m-fc").textContent = `${formatSci(diag[5], 3)} MHz`;
    $("m-rms").textContent = `${diag[6].toFixed(2)} / ${diag[7].toFixed(2)} mm`;
    $("m-emit").textContent = `${formatSci(diag[8], 3)} mm*mrad`;
    $("m-twiss").textContent = `${diag[11].toFixed(2)} m / ${diag[10].toFixed(2)}`;
    $("m-drift").textContent = diag[13].toExponential(2);
    $("m-power").textContent = `${formatSci(diag[14], 3)} W`;
    $("m-turns").textContent = `${diag[15].toFixed(1)} | ${(diag[17] * 100).toFixed(0)}%`;

    $("orbit-meta").textContent = `R = ${diag[19].toFixed(3)} m`;
    $("twiss-meta").textContent = `|bg - a^2 - 1| = ${diag[18].toExponential(2)}`;

    renderOrbit(particles, count, diag[19], Number(ctrls.mode.value));
    renderPhaseSpace(phase, count, diag);
    renderRfSpace(phase, count);
    renderHistory(diag);
  }

  function tick() {
    if (running && wasm) {
      wasm.sim_step(Number(ctrls.substeps.value));
      renderFrame();
    }
    requestAnimationFrame(tick);
  }

  function renderMarkdown(md) {
    const escaped = md
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
    return escaped
      .replace(/```([\s\S]*?)```/g, "<pre><code>$1</code></pre>")
      .replace(/^### (.*)$/gm, "<h3>$1</h3>")
      .replace(/^## (.*)$/gm, "<h2>$1</h2>")
      .replace(/^# (.*)$/gm, "<h1>$1</h1>")
      .replace(/\*\*(.*?)\*\*/g, "<strong>$1</strong>")
      .replace(/`([^`]+)`/g, "<code>$1</code>")
      .replace(/\n\n/g, "<br/><br/>")
      .replace(/\n/g, "<br/>");
  }

  async function openTab(tabName) {
    document.querySelectorAll(".tab-btn").forEach((b) => {
      b.classList.toggle("active", b.dataset.tab === tabName);
    });
    if (tabName === "sim") {
      $("tab-sim").classList.remove("hidden");
      $("tab-docs").classList.add("hidden");
      return;
    }
    $("tab-sim").classList.add("hidden");
    $("tab-docs").classList.remove("hidden");
    $("doc-filename").textContent = `docs/${tabName}`;
    $("doc-raw-link").href = `docs/${tabName}`;
    $("doc-content").innerHTML = "Loading document...";
    try {
      const resp = await fetch(`docs/${tabName}`);
      const text = await resp.text();
      $("doc-content").innerHTML = renderMarkdown(text);
    } catch (err) {
      $("doc-content").textContent = `Unable to load docs/${tabName}: ${err}`;
    }
  }

  async function bootWasm() {
    const pill = $("engine-status");
    try {
      let instance;
      if (WebAssembly.instantiateStreaming) {
        try {
          const res = await WebAssembly.instantiateStreaming(
            fetch("particle_accelarator.wasm"),
            {}
          );
          instance = res.instance;
        } catch (_) {
          const bytes = await (await fetch("particle_accelarator.wasm")).arrayBuffer();
          instance = (await WebAssembly.instantiate(bytes, {})).instance;
        }
      } else {
        const bytes = await (await fetch("particle_accelarator.wasm")).arrayBuffer();
        instance = (await WebAssembly.instantiate(bytes, {})).instance;
      }
      wasm = instance.exports;
      pill.textContent = "WASM f64 Engine Active";
      pill.classList.remove("loading");
      pill.classList.add("ready");
      resetSimulation();
      requestAnimationFrame(tick);
    } catch (err) {
      pill.textContent = "WASM Load Error";
      console.error(err);
    }
  }

  Object.values(ctrls).forEach((el) => {
    el.addEventListener("input", () => {
      if (el === ctrls.substeps) syncLabels();
      else resetSimulation();
    });
  });

  $("btn-toggle").addEventListener("click", () => {
    running = !running;
    $("btn-toggle").textContent = running ? "Pause" : "Resume";
  });

  $("btn-step").addEventListener("click", () => {
    if (wasm) {
      wasm.sim_step(Number(ctrls.substeps.value));
      renderFrame();
    }
  });

  $("btn-reset").addEventListener("click", resetSimulation);

  document.querySelectorAll(".tab-btn").forEach((btn) => {
    btn.addEventListener("click", () => openTab(btn.dataset.tab));
  });

  syncLabels();
  bootWasm();
})();

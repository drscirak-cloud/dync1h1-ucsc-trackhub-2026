(function () {
  "use strict";

  const AA3 = {
    A: "Ala", C: "Cys", D: "Asp", E: "Glu", F: "Phe",
    G: "Gly", H: "His", I: "Ile", K: "Lys", L: "Leu",
    M: "Met", N: "Asn", P: "Pro", Q: "Gln", R: "Arg",
    S: "Ser", T: "Thr", V: "Val", W: "Trp", Y: "Tyr"
  };

  let viewer3d = null;
  let currentRes = 3384;
  let currentAlt = "";
  let selectedVariant = null;
  let filteredVariants = [];

  function getBundle() {
    return window.DYNC1H1_EXPLORER_DATA || null;
  }

  function domainForAa(bundle, aaPos) {
    for (let i = 0; i < bundle.domains.length; i++) {
      const d = bundle.domains[i];
      if (aaPos >= d.start && aaPos <= d.end) return d;
    }
    return bundle.domains[0] || { name: "Tail", region: "Tail", color: "#58a6ff", start: 1, end: 4646 };
  }

  function exonForAa(bundle, aaPos) {
    const cdsNt = (aaPos - 1) * 3 + 1;
    let cum = 0;
    for (let i = 0; i < bundle.exons.length; i++) {
      const ex = bundle.exons[i];
      cum += ex.cds_length_nt;
      if (cdsNt <= cum) return ex;
    }
    return bundle.exons[bundle.exons.length - 1] || { exon_number: 78, cds_length_nt: 129, exon_frame_class: "Terminal_Exon78" };
  }

  function initTheme() {
    const saved = localStorage.getItem("dync1h1_theme");
    if (saved === "dark" || saved === "light") {
      document.documentElement.setAttribute("data-theme", saved);
    }
    document.getElementById("btn-theme-toggle").addEventListener("click", () => {
      const cur = document.documentElement.getAttribute("data-theme");
      const next = cur === "dark" ? "light" : "dark";
      document.documentElement.setAttribute("data-theme", next);
      localStorage.setItem("dync1h1_theme", next);
    });
  }

  function initDialog() {
    const dlg = document.getElementById("dialog-ucsc");
    document.getElementById("btn-ucsc-modal").addEventListener("click", () => dlg.showModal());
    document.getElementById("btn-close-modal").addEventListener("click", () => dlg.close());
    dlg.addEventListener("click", (ev) => {
      const rect = dlg.getBoundingClientRect();
      const inside = (
        ev.clientX >= rect.left && ev.clientX <= rect.right &&
        ev.clientY >= rect.top && ev.clientY <= rect.bottom
      );
      if (!inside) dlg.close();
    });
  }

  function buildHomodimerPdb(bundle, stateKey) {
    const coords = stateKey === "OFF" ? bundle.ca_off_phi_9292 : bundle.ca_on_dda_9292;
    const sens = bundle.pure3d_sens_4646;
    const lines = [];
    let serial = 1;
    for (let chainIdx = 0; chainIdx < 2; chainIdx++) {
      const chainId = chainIdx === 0 ? "A" : "B";
      const baseIdx = chainIdx * 4646;
      for (let i = 0; i < 4646; i++) {
        const aaPos = i + 1;
        const xyz = coords[baseIdx + i];
        if (!xyz) continue;
        const bfac = Math.min(99.99, Math.max(0, (sens[i] || 0) * 100));
        const line = (
          "ATOM  " +
          String(serial).padStart(5, " ") +
          "  CA  ALA " +
          chainId +
          String(aaPos).padStart(4, " ") +
          "    " +
          Number(xyz[0]).toFixed(3).padStart(8, " ") +
          Number(xyz[1]).toFixed(3).padStart(8, " ") +
          Number(xyz[2]).toFixed(3).padStart(8, " ") +
          "  1.00" +
          bfac.toFixed(2).padStart(6, " ") +
          "           C  "
        );
        lines.push(line);
        serial++;
      }
      lines.push("TER");
    }
    lines.push("END");
    return lines.join("\n");
  }

  function valueToHex(val, minV, maxV, lowHex, midHex, highHex) {
    const t = Math.max(0, Math.min(1, (val - minV) / Math.max(1e-6, maxV - minV)));
    function hex2rgb(h) {
      const n = parseInt(h.replace("#", ""), 16);
      return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
    }
    const c1 = t < 0.5 ? hex2rgb(lowHex) : hex2rgb(midHex);
    const c2 = t < 0.5 ? hex2rgb(midHex) : hex2rgb(highHex);
    const u = t < 0.5 ? t * 2 : (t - 0.5) * 2;
    const r = Math.round(c1[0] + (c2[0] - c1[0]) * u);
    const g = Math.round(c1[1] + (c2[1] - c1[1]) * u);
    const b = Math.round(c1[2] + (c2[2] - c1[2]) * u);
    return "#" + [r, g, b].map((v) => v.toString(16).padStart(2, "0")).join("");
  }

  function render3DHomodimer() {
    const bundle = getBundle();
    if (!bundle || typeof $3Dmol === "undefined") return;
    const container = document.getElementById("mol-viewer");
    if (!viewer3d) {
      viewer3d = $3Dmol.createViewer(container, { backgroundColor: "#0d1117" });
    }
    viewer3d.clear();
    const stateKey = document.getElementById("select-state").value;
    const colorMode = document.getElementById("select-color-mode").value;
    const pdbText = buildHomodimerPdb(bundle, stateKey);
    viewer3d.addModel(pdbText, "pdb");

    viewer3d.setStyle({}, {
      trace: {
        radius: 0.55,
        colorfunc: function (atom) {
          const aaPos = atom.resi;
          const idx = aaPos - 1;
          if (idx < 0 || idx >= 4646) return "#888888";
          if (colorMode === "domain") {
            return domainForAa(bundle, aaPos).color || "#58a6ff";
          } else if (colorMode === "d3d") {
            return valueToHex(bundle.pure3d_sens_4646[idx] || 0, 0.08, 0.55, "#2b5c8f", "#f4d35e", "#e63946");
          } else if (colorMode === "gnm") {
            return valueToHex(bundle.gnm_strain_4646[idx] || 0, 0.0, 0.45, "#1d3557", "#48cae4", "#ff6b6b");
          } else {
            return valueToHex(bundle.loco_oe_4646[idx] || 0, 0.20, 0.85, "#3a86ff", "#e9c46a", "#9b5de5");
          }
        }
      }
    });

    if (currentRes >= 1 && currentRes <= 4646) {
      const lbl = selectedVariant && (selectedVariant.hgvs_p || selectedVariant.hgvs_c)
        ? (selectedVariant.hgvs_p || selectedVariant.hgvs_c)
        : `Codon ${currentRes}${currentAlt ? " -> " + (AA3[currentAlt] || currentAlt) : ""}`;
      viewer3d.addStyle({ resi: currentRes }, { sphere: { radius: 4.8, color: "#ffd166" } });
      viewer3d.addLabel(lbl + " (Chain A)", {
        fontSize: 12,
        fontColor: "#0d1117",
        backgroundColor: "#ffd166",
        backgroundOpacity: 0.92
      }, { chain: "A", resi: currentRes });
    }

    viewer3d.zoomTo();
    viewer3d.render();
    updateViewerLegend(colorMode, bundle);
  }

  function updateViewerLegend(colorMode, bundle) {
    const leg = document.getElementById("viewer-legend");
    if (colorMode === "domain") {
      const seen = new Set();
      const items = [];
      bundle.domains.forEach((d) => {
        if (!seen.has(d.region)) {
          seen.add(d.region);
          items.push(`<span><i class="legend-dot" style="background:${d.color}"></i>${d.region}</span>`);
        }
      });
      items.push(`<span><i class="legend-dot" style="background:#ffd166"></i>Selected Codon (Chains A &amp; B)</span>`);
      leg.innerHTML = items.join("");
    } else {
      const labels = {
        d3d: "Low (Blue) &rarr; Intermediate (Gold) &rarr; High Track 1A Pure 3D Sensitivity (Red, 0 Cons)",
        gnm: "Low Strain (Navy) &rarr; Collective 9,292-Node GNM Hinge Peak (Cyan/Coral)",
        oe: "Tolerated Window (Blue) &rarr; Strong LOCO 31-Codon CpG Depletion (Violet)"
      };
      leg.innerHTML = `<span>${labels[colorMode] || ""}</span> &middot; <span><i class="legend-dot" style="background:#ffd166"></i>Selected Codon</span>`;
    }
  }

  function shortQuad(q) {
    if (!q) return "Q_IV";
    if (q.startsWith("Q_I_")) return "Q I (3D+ / T2+)";
    if (q.startsWith("Q_IIa")) return "Q IIa (Splice/PTV/NMD)";
    if (q.startsWith("Q_IIb")) return "Q IIb (Constrained Missense)";
    if (q.startsWith("Q_III")) return "Q III (3D Pocket/Seam)";
    return "Q IV (Concordant Low)";
  }

  function updateInspector() {
    const bundle = getBundle();
    if (!bundle) return;
    const idx = Math.max(0, Math.min(4645, currentRes - 1));
    const dom = domainForAa(bundle, currentRes);
    const ex = exonForAa(bundle, currentRes);

    const p3d = bundle.pure3d_sens_4646[idx] || 0;
    const gnm = bundle.gnm_strain_4646[idx] || 0;
    const loco = bundle.loco_oe_4646[idx] || 0;

    const v = selectedVariant;
    const titleStr = v
      ? `${v.hgvs_p || v.hgvs_c || v.id} (${v.tier})`
      : `Codon ${currentRes}${currentAlt ? " -> " + (AA3[currentAlt] || currentAlt) : " (Wild-Type Profile)"}`;
    document.getElementById("inspector-title").textContent = titleStr;
    document.getElementById("inspector-domain").textContent =
      `${dom.name} (${dom.region}) \u00B7 Exon ${ex.exon_number}/78 (${ex.exon_frame_class.replace(/_/g, " ")})`;

    const t1a = v ? v.t1_pure3d : p3d;
    const t1b = v ? v.t1_cons : p3d;
    const t2 = v ? v.t2_intr : loco;
    const joint = v ? v.m3_joint : (0.45 * p3d + 0.55 * loco);
    const quadRaw = v ? v.quad : (
      t1a >= 0.5 && t2 >= 0.5 ? "Q_I_Concordant" :
      t1a < 0.5 && t2 >= 0.5 ? "Q_IIb_Constrained" :
      t1a >= 0.5 && t2 < 0.5 ? "Q_III_3D" : "Q_IV_Low"
    );

    const badge = document.getElementById("inspector-quadrant-badge");
    badge.className = "badge " + (
      quadRaw.startsWith("Q_I_") ? "badge-q1" :
      quadRaw.startsWith("Q_II") ? "badge-q2" :
      quadRaw.startsWith("Q_III") ? "badge-q3" : "badge-q4"
    );
    badge.textContent = shortQuad(quadRaw);

    document.getElementById("dl-track1").innerHTML = `
      <dt>Structural Domain</dt><dd>${dom.name} (${dom.start}&ndash;${dom.end} AA)</dd>
      <dt>Codon Pure 3D Sensitivity</dt><dd>${p3d.toFixed(3)}</dd>
      <dt>9,292-AA GNM Hinge Strain</dt><dd>${gnm.toFixed(4)}</dd>
      <dt>Track 1A Pure 3D OOF</dt><dd>${t1a.toFixed(3)}</dd>
      <dt>Track 1B (3D + Cons) OOF</dt><dd>${t1b.toFixed(3)}</dd>
    `;

    document.getElementById("dl-track2").innerHTML = `
      <dt>LOCO 31-Codon CpG (1&minus;O/E)</dt><dd>${loco.toFixed(3)}</dd>
      <dt>Canonical Exon (${ex.exon_number}/78)</dt><dd>${ex.cds_length_nt} nt (${ex.is_internal_frame_preserving_3n ? "Internal 3n" : (ex.exon_number === 78 ? "Terminal 3n" : "Phase-Shifting")})</dd>
      <dt>Empirical GTEx WT/WT Swing</dt><dd>Muscle 7.70 pp &middot; Nerve 4.28 pp &middot; Cortex 0.00 pp</dd>
      <dt>Track 2-Intrinsic OOF (0 AF)</dt><dd>${t2.toFixed(3)}</dd>
      <dt>Joint Dual-Track M3B OOF</dt><dd><strong>${joint.toFixed(3)}</strong></dd>
    `;

    let assayRec = "Low dual-track pathogenicity probability (Q IV); deprioritize in clinical reporting.";
    if (quadRaw.startsWith("Q_IIa")) {
      assayRec = "Q IIa (SpliceAI >= 0.20 / PTV / NMD): Route to patient blood or fibroblast RT-PCR, RNA-seq, or minigene splicing confirmation.";
    } else if (quadRaw.startsWith("Q_IIb")) {
      assayRec = "Q IIb (Constrained coding missense / regulatory, SpliceAI < 0.20): Route to cargo-adaptor co-immunoprecipitation (BICD2 / HOOK3 / JIP3) or recombinant dynein assay.";
    } else if (quadRaw.startsWith("Q_III")) {
      assayRec = "Q III (3D pocket, stalk heptad, or Phi seam in lower-depletion exon): Route to single-molecule TIRF motility, LIS1 arrest, or DSF thermal unfolding assay.";
    } else if (quadRaw.startsWith("Q_I_") || joint >= 0.80) {
      assayRec = "Q I Concordant High-Confidence Pathogenic Profile: Prioritize clinical segregation and phenotype correlation across cortical MRI and neuromuscular electrophysiology.";
    }
    document.getElementById("synthesis-text").textContent = assayRec;
  }

  function populateSelectors() {
    const bundle = getBundle();
    if (!bundle) return;
    const altSel = document.getElementById("select-alt-aa");
    Object.keys(AA3).sort().forEach((k) => {
      const opt = document.createElement("option");
      opt.value = k;
      opt.textContent = `Substitute -> ${AA3[k]} (${k})`;
      altSel.appendChild(opt);
    });

    const presetSel = document.getElementById("select-preset");
    const archetypes = [
      { label: "p.Glu1518Lys (9E0Z LIS1 Pocket, +7.06 kcal/mol, 0 clashes)", aa: 1518 },
      { label: "p.Gly2598Val (9BLZ AAA3 Turn, +21.94 kcal/mol, 0 clashes)", aa: 2598 },
      { label: "p.Pro3173Arg (9BMD Stalk CC1 Repaired, +13.43 kcal/mol, 0 clashes)", aa: 3173 },
      { label: "p.Arg3384Gln (MTBD Microtubule Cleft Archetype)", aa: 3384 },
      { label: "p.Arg3695Trp (9BLZ AAA5 H-Bond Loss, +1.68 kcal/mol)", aa: 3695 },
      { label: "p.Leu196Ser (5OWO N-Terminal Dimerization Core)", aa: 196 }
    ];
    archetypes.forEach((p, idx) => {
      const opt = document.createElement("option");
      opt.value = String(idx);
      opt.textContent = p.label;
      presetSel.appendChild(opt);
    });
    presetSel.addEventListener("change", () => {
      const idx = parseInt(presetSel.value, 10);
      if (isNaN(idx) || !archetypes[idx]) return;
      const targetAa = archetypes[idx].aa;
      const match = bundle.variants.find((v) => v.aa === targetAa);
      if (match) {
        selectVariantRow(match);
      } else {
        selectedVariant = null;
        currentRes = targetAa;
        document.getElementById("input-residue").value = currentRes;
        updateInspector();
        render3DHomodimer();
        drawLinearSvgCursor();
      }
    });
  }

  function drawLinearSvg() {
    const bundle = getBundle();
    if (!bundle) return;
    const svg = document.getElementById("linear-profile-svg");
    const W = 1160, padL = 52, padR = 16;
    const plotW = W - padL - padR;
    const xOf = (aa) => padL + ((aa - 1) / 4645) * plotW;

    let html = "";
    bundle.domains.forEach((d) => {
      const x1 = xOf(d.start);
      const x2 = xOf(d.end);
      html += `<rect x="${x1.toFixed(1)}" y="8" width="${Math.max(1.5, x2 - x1).toFixed(1)}" height="20" rx="2" fill="${d.color}" opacity="0.88"><title>${d.name} (${d.start}-${d.end} AA)</title></rect>`;
      if (x2 - x1 > 52) {
        html += `<text x="${((x1 + x2) / 2).toFixed(1)}" y="22" font-size="9.2" fill="#ffffff" text-anchor="middle" font-weight="600">${d.name.replace(/_/g, " ")}</text>`;
      }
    });

    function buildPath(arr, yTop, yBot, maxVal) {
      const step = 12;
      const pts = [];
      for (let i = 0; i < arr.length; i += step) {
        let sum = 0, cnt = 0;
        for (let k = i; k < Math.min(arr.length, i + step); k++) {
          sum += arr[k] || 0;
          cnt++;
        }
        const v = Math.min(1, Math.max(0, (sum / cnt) / maxVal));
        const x = xOf(i + 1 + step / 2);
        const y = yBot - v * (yBot - yTop);
        pts.push(`${x.toFixed(1)},${y.toFixed(1)}`);
      }
      return pts.join(" ");
    }

    const pT1 = buildPath(bundle.pure3d_sens_4646, 40, 92, 0.65);
    const pGNM = buildPath(bundle.gnm_strain_4646, 102, 150, 0.55);
    const pT2 = buildPath(bundle.loco_oe_4646, 160, 208, 1.00);

    html += `<text x="6" y="70" font-size="9.5" fill="currentColor" font-weight="600">T1A 3D</text>`;
    html += `<polyline fill="none" stroke="#1b9e77" stroke-width="1.6" points="${pT1}"/>`;
    html += `<line x1="${padL}" y1="93" x2="${W - padR}" y2="93" stroke="#888" stroke-opacity="0.25"/>`;

    html += `<text x="6" y="130" font-size="9.5" fill="currentColor" font-weight="600">GNM v2</text>`;
    html += `<polyline fill="none" stroke="#e76f51" stroke-width="1.5" points="${pGNM}"/>`;
    html += `<line x1="${padL}" y1="151" x2="${W - padR}" y2="151" stroke="#888" stroke-opacity="0.25"/>`;

    html += `<text x="6" y="188" font-size="9.5" fill="currentColor" font-weight="600">LOCO O/E</text>`;
    html += `<polyline fill="none" stroke="#7570b3" stroke-width="1.6" points="${pT2}"/>`;

    [1, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4646].forEach((tick) => {
      const xt = xOf(tick);
      html += `<text x="${xt.toFixed(1)}" y="222" font-size="9.5" fill="currentColor" opacity="0.7" text-anchor="middle">${tick}</text>`;
    });

    html += `<line id="svg-cursor-line" x1="${xOf(currentRes).toFixed(1)}" y1="4" x2="${xOf(currentRes).toFixed(1)}" y2="210" stroke="#ffd166" stroke-width="2.2"/>`;
    svg.innerHTML = html;

    svg.addEventListener("click", (ev) => {
      const rect = svg.getBoundingClientRect();
      const relX = ((ev.clientX - rect.left) / rect.width) * W;
      const aa = Math.round(1 + ((relX - padL) / plotW) * 4645);
      if (aa >= 1 && aa <= 4646) {
        selectedVariant = null;
        currentRes = aa;
        document.getElementById("input-residue").value = currentRes;
        updateInspector();
        render3DHomodimer();
        drawLinearSvgCursor();
      }
    });
  }

  function drawLinearSvgCursor() {
    const ln = document.getElementById("svg-cursor-line");
    if (!ln) return;
    const W = 1160, padL = 52, padR = 16;
    const x = padL + ((currentRes - 1) / 4645) * (W - padL - padR);
    ln.setAttribute("x1", x.toFixed(1));
    ln.setAttribute("x2", x.toFixed(1));
  }

  function updateStoichiometry() {
    const nmdPct = parseFloat(document.getElementById("range-nmd").value);
    const eqtlPct = parseFloat(document.getElementById("range-eqtl").value);
    document.getElementById("lbl-nmd-val").textContent = `${nmdPct.toFixed(0)}%`;
    document.getElementById("lbl-eqtl-val").textContent = `${eqtlPct >= 0 ? "+" : ""}${eqtlPct.toFixed(0)}%`;

    const fNmd = nmdPct / 100.0;
    const beta = eqtlPct / 100.0;
    const eMut = Math.max(0, (1.0 + beta) * (1.0 - fNmd));
    const pMut = eMut / (1.0 + eMut);
    const fWTWT = Math.pow(1.0 - pMut, 2) * 100.0;
    const fWTMut = 2.0 * pMut * (1.0 - pMut) * 100.0;
    const fMutMut = Math.pow(pMut, 2) * 100.0;
    const totalDose = ((1.0 + eMut) / 2.0) * 100.0;

    document.getElementById("stoich-bars").innerHTML = `
      <div class="bar-item">
        <div class="bar-top"><span>Intact Wild-Type Homodimers (WT / WT)</span><strong>${fWTWT.toFixed(2)}%</strong></div>
        <div class="bar-track"><div class="bar-fill" style="width:${fWTWT.toFixed(1)}%;background:#1b9e77"></div></div>
      </div>
      <div class="bar-item">
        <div class="bar-top"><span>Heterodimeric Poisoned Complexes (WT / Mut)</span><strong>${fWTMut.toFixed(2)}%</strong></div>
        <div class="bar-track"><div class="bar-fill" style="width:${fWTMut.toFixed(1)}%;background:#f4a261"></div></div>
      </div>
      <div class="bar-item">
        <div class="bar-top"><span>Biallelic Mutant Homodimers (Mut / Mut)</span><strong>${fMutMut.toFixed(2)}%</strong></div>
        <div class="bar-track"><div class="bar-fill" style="width:${fMutMut.toFixed(1)}%;background:#e63946"></div></div>
      </div>
      <p class="muted">Total steady-state DYNC1H1 subunit pool: <strong>${totalDose.toFixed(1)}%</strong> of diploid baseline (poisoned dimer fraction: <strong>${(fWTMut + fMutMut).toFixed(2)}%</strong>).</p>
    `;
  }

  function selectVariantRow(v) {
    selectedVariant = v;
    if (v.aa && v.aa >= 1 && v.aa <= 4646) {
      currentRes = v.aa;
      document.getElementById("input-residue").value = currentRes;
    }
    currentAlt = "";
    document.getElementById("select-alt-aa").value = "";
    updateInspector();
    render3DHomodimer();
    drawLinearSvgCursor();
  }

  function renderVariantTable() {
    const bundle = getBundle();
    if (!bundle) return;
    const q = document.getElementById("input-search").value.trim().toLowerCase();
    const tierFilter = document.getElementById("filter-tier").value;

    filteredVariants = bundle.variants.filter((v) => {
      if (tierFilter === "ClinVar_PLP_1StarPlus_Published" && v.tier !== "ClinVar_PLP_1StarPlus_Published") return false;
      if (tierFilter === "ClinVar_PLP_1StarPlus" && !v.tier.startsWith("ClinVar_PLP_1StarPlus")) return false;
      if (tierFilter === "Published_DYNC1H1_LOVD_Verified" && v.tier !== "Published_DYNC1H1_LOVD_Verified") return false;
      if (tierFilter === "VUS_High_Priority" && !(v.tier === "ClinVar_VUS_or_Conflicting" && v.m3_joint >= 0.80)) return false;
      if (tierFilter === "Q_IIa" && !v.quad.startsWith("Q_IIa")) return false;
      if (tierFilter === "Q_IIb" && !v.quad.startsWith("Q_IIb")) return false;
      if (tierFilter === "Q_III" && !v.quad.startsWith("Q_III")) return false;
      if (tierFilter === "ClinVar_BLB_1StarPlus" && v.tier !== "ClinVar_BLB_1StarPlus") return false;

      if (!q) return true;
      const hay = `${v.hgvs_p} ${v.hgvs_c} ${v.pos} exon${v.exon} exon ${v.exon} ${v.dom} ${v.tier} ${v.quad}`.toLowerCase();
      return hay.includes(q);
    });

    document.getElementById("table-count-label").textContent =
      `Showing ${Math.min(250, filteredVariants.length)} of ${filteredVariants.length.toLocaleString()} matching variants (click any row to load into 3D Homodimer Inspector)`;

    const tbody = document.getElementById("variant-tbody");
    tbody.innerHTML = "";
    filteredVariants.slice(0, 250).forEach((v) => {
      const tr = document.createElement("tr");
      const label = v.hgvs_p || v.hgvs_c || v.id;
      tr.innerHTML = `
        <td><strong>${label}</strong>${v.hgvs_p && v.hgvs_c ? `<br><small class="muted">${v.hgvs_c}</small>` : ""}</td>
        <td>chr14:${v.pos}</td>
        <td>Exon ${v.exon} &middot; ${v.dom}</td>
        <td>${v.tier.replace(/_/g, " ")}</td>
        <td>${v.t1_pure3d.toFixed(3)}</td>
        <td>${v.t1_cons.toFixed(3)}</td>
        <td>${v.t2_intr.toFixed(3)}</td>
        <td><strong>${v.m3_joint.toFixed(3)}</strong></td>
        <td>${shortQuad(v.quad)}</td>
      `;
      tr.addEventListener("click", () => selectVariantRow(v));
      tbody.appendChild(tr);
    });
  }

  function exportFilteredCsv() {
    if (!filteredVariants.length) return;
    const cols = ["id", "hgvs_p", "hgvs_c", "pos", "aa", "exon", "dom", "tier", "t1_pure3d", "t1_cons", "t2_intr", "m3_joint", "quad"];
    const rows = [cols.join(",")];
    filteredVariants.forEach((v) => {
      rows.push(cols.map((c) => `"${String(v[c] ?? "").replace(/"/g, '""')}"`).join(","));
    });
    const blob = new Blob([rows.join("\n")], { type: "text/csv;charset=utf-8;" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = "DYNC1H1_DualTrack_v2_Filtered_Variants.csv";
    a.click();
    URL.revokeObjectURL(a.href);
  }

  window.addEventListener("DOMContentLoaded", () => {
    initTheme();
    initDialog();

    document.getElementById("btn-inspect-residue").addEventListener("click", () => {
      const pos = Math.max(1, Math.min(4646, parseInt(document.getElementById("input-residue").value, 10) || 1));
      selectedVariant = null;
      currentRes = pos;
      currentAlt = document.getElementById("select-alt-aa").value;
      document.getElementById("input-residue").value = pos;
      updateInspector();
      render3DHomodimer();
      drawLinearSvgCursor();
    });

    document.getElementById("select-state").addEventListener("change", render3DHomodimer);
    document.getElementById("select-color-mode").addEventListener("change", render3DHomodimer);
    document.getElementById("input-search").addEventListener("input", renderVariantTable);
    document.getElementById("filter-tier").addEventListener("change", renderVariantTable);
    document.getElementById("btn-export-csv").addEventListener("click", exportFilteredCsv);

    document.getElementById("range-nmd").addEventListener("input", updateStoichiometry);
    document.getElementById("range-eqtl").addEventListener("input", updateStoichiometry);
    document.querySelectorAll(".pill").forEach((btn) => {
      btn.addEventListener("click", () => {
        document.getElementById("range-nmd").value = btn.dataset.nmd;
        document.getElementById("range-eqtl").value = btn.dataset.eqtl;
        updateStoichiometry();
      });
    });

    updateStoichiometry();
    populateSelectors();
    drawLinearSvg();
    renderVariantTable();
    updateInspector();
    if (typeof $3Dmol !== "undefined") {
      render3DHomodimer();
    } else {
      setTimeout(render3DHomodimer, 650);
    }
  });
})();

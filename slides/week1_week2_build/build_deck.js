// Week 1 + Week 2 project update deck (pptxgenjs). Run: node build_deck.js <out.pptx> <applyThemePath>
const pptxgen = require("pptxgenjs");
const path = require("path");
const OUT = process.argv[2] || "week1_week2_update.pptx";
const { applyTheme } = require(process.argv[3]);
const DIR = __dirname;

const THEME = {
  name: "Hematoxylin and Eosin",
  headFontFace: "Cambria",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "1F1A2E", lt1: "FFFFFF", dk2: "3D2B6B", lt2: "F3EEF8",
    accent1: "3D2B6B", accent2: "C2456E", accent3: "1F8A7E", accent4: "C77C02", accent5: "8C7AA8", accent6: "5F6B7A",
    hlink: "1F8A7E", folHlink: "8C7AA8",
  },
};
const HEX = THEME.colors;

const pres = new pptxgen();
pres.layout = "LAYOUT_16x9"; // 10 x 5.625 in
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "Which Truth Does the Slide See? Week 1 + Week 2 update";
pres.author = "Aman Uniyal, Aman Srivastava";
const C = pres.SchemeColor;
// semantic colours
const H = C.accent1, EOS = C.accent2, RNA = C.accent3, IHC = C.accent4, MUTED = C.accent6;

// ---------------- layouts
pres.defineSlideMaster({
  title: "TITLE",
  background: { color: C.text2 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.6, y: 1.35, w: 8.8, h: 1.1, fontSize: 38, bold: true,
      color: C.background1, valign: "bottom", align: "left" }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.6, y: 2.5, w: 8.8, h: 0.8, fontSize: 17,
      color: C.background2, valign: "top", align: "left" }, text: "" } },
  ],
});
pres.defineSlideMaster({
  title: "CONTENT",
  background: { color: C.background1 },
  margin: [0.3, 0.4, 0.4, 0.4],
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.4, y: 0.22, w: 7.9, h: 0.62, fontSize: 24, bold: true,
      color: C.text2, valign: "middle", align: "left" }, text: "" } },
    { text: { text: "Which Truth Does the Slide See?  ·  TCGA-BRCA  ·  Uniyal & Srivastava",
      options: { x: 0.4, y: 5.28, w: 6.5, h: 0.25, fontSize: 9, color: C.accent6, align: "left" } } },
  ],
  slideNumber: { x: 9.15, y: 5.28, w: 0.45, h: 0.25, fontSize: 9, color: C.accent6, align: "right" },
});

// ---------------- helpers
let objN = 0;
const nm = (s) => `${s}_${++objN}`;
function pill(slide, text, fill) {
  slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 8.45, y: 0.33, w: 1.15, h: 0.36, rectRadius: 0.18,
    fill: { color: fill }, line: { color: fill }, objectName: nm("week_tag") });
  slide.addText(text, { x: 8.45, y: 0.33, w: 1.15, h: 0.36, fontSize: 10.5, bold: true, color: C.background1,
    align: "center", valign: "middle", margin: 0, isTextBox: true, objectName: nm("week_tag_text") });
}
function card(slide, x, y, w, h, title, items, o = {}) {
  slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.08, fill: { color: o.fill || C.background2 },
    line: { color: o.line || o.fill || C.background2, width: 0.75 }, objectName: nm("card") });
  if (o.dot) {
    slide.addShape(pres.shapes.OVAL, { x: x + 0.14, y: y + 0.13, w: 0.2, h: 0.2, fill: { color: o.dot }, line: { color: o.dot },
      objectName: nm("card_dot") });
  }
  const tx = x + (o.dot ? 0.42 : 0.15);
  slide.addText(title, { x: tx, y: y + 0.07, w: w - (tx - x) - 0.1, h: 0.32, fontSize: o.titleSize || 13, bold: true,
    color: o.titleColor || C.text2, valign: "middle", margin: 0, isTextBox: true, objectName: nm("card_title") });
  if (items && items.length) {
    const runs = items.map((it, i) => {
      const t = typeof it === "string" ? [{ text: it }] : it;
      return t.map((r, j) => ({ text: r.text, options: Object.assign({ bullet: j === 0 && !o.plain ? { indent: 12 } : false },
        r.options || {}, j === t.length - 1 && i < items.length - 1 ? { breakLine: true } : {}) }));
    }).flat();
    slide.addText(runs, { x: x + 0.12, y: y + 0.42, w: w - 0.22, h: h - 0.5, fontSize: o.size || 11.5, color: C.text1,
      valign: "top", paraSpaceAfter: 3, margin: [0, 2, 0, 0], isTextBox: true, objectName: nm("card_body") });
  }
}
function stat(slide, x, y, w, big, label, color) {
  slide.addText(big, { x, y, w, h: 0.5, fontSize: 26, bold: true, color, align: "left", valign: "bottom", margin: 0,
    isTextBox: true, objectName: nm("stat_value") });
  slide.addText(label, { x, y: y + 0.5, w, h: 0.42, fontSize: 10, color: C.accent6, align: "left", valign: "top", margin: 0,
    isTextBox: true, objectName: nm("stat_label") });
}
function arrow(slide, x1, y1, x2, y2) {
  // PowerPoint rejects negative extents: draw with positive w/h and flip instead
  const o = { x: Math.min(x1, x2), y: Math.min(y1, y2), w: Math.abs(x2 - x1), h: Math.abs(y2 - y1),
    line: { color: C.accent6, width: 1.5, endArrowType: "triangle" }, objectName: nm("arrow") };
  if (x2 < x1) o.flipH = true;
  if (y2 < y1) o.flipV = true;
  slide.addShape(pres.shapes.LINE, o);
}
const hdr = (t) => ({ text: t, options: { bold: true, color: C.background1, fill: { color: C.text2 }, align: "center", valign: "middle" } });
const cell = (t, o = {}) => ({ text: t, options: Object.assign({ align: "center", valign: "middle" }, o) });
const lcell = (t, o = {}) => cell(t, Object.assign({ align: "left" }, o));

// ================= 1. Title
pres.addSection({ title: "Title" });
let s = pres.addSlide({ masterName: "TITLE", sectionTitle: "Title" });
s.addText("Which Truth Does the Slide See?", { placeholder: "title" });
s.addText("Auditing what pathology foundation models learn when IHC, RNA and protein disagree", { placeholder: "body" });
s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.6, y: 0.55, w: 2.9, h: 0.42, rectRadius: 0.21, fill: { color: EOS }, line: { color: EOS },
  objectName: "title_tag" });
s.addText("Project update  ·  Week 1 + Week 2", { x: 0.6, y: 0.55, w: 2.9, h: 0.42, fontSize: 13, bold: true, color: C.background1,
  align: "center", valign: "middle", margin: 0, isTextBox: true, objectName: "title_tag_text" });
s.addText([
  { text: "Aman Uniyal  ·  Aman Srivastava", options: { bold: true, breakLine: true } },
  { text: "Topics in AI for Healthcare, IIIT Hyderabad", options: { breakLine: true } },
  { text: "Base project: Kamineni – Multimodal Breast Cancer Diagnostic & Outcome Assistant" },
], { x: 0.6, y: 3.9, w: 8.8, h: 1.0, fontSize: 13, color: C.background2, valign: "top", margin: 0, paraSpaceAfter: 3,
  isTextBox: true, objectName: "title_authors" });
s.addNotes("Week 1 slides were not presented, so this deck covers Week 1 (proposal, data, first baselines) and the Week 2 updates, ending with the Week 3 plan.");

// ================= 2. Proposal
pres.addSection({ title: "Week 1: proposal and baselines" });
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Week 1: proposal and baselines" });
s.addText("Proposal: problem, question and pipeline", { placeholder: "title" });
pill(s, "WEEK 1", H);
card(s, 0.4, 1.0, 4.55, 1.55, "Clinical problem", [
  [{ text: "ER/PR/HER2 status drives therapy and is measured by " }, { text: "IHC", options: { bold: true, color: IHC } }],
  "IHC often disagrees with mRNA; for HER2 0 vs 1+ pathologists agree in only ~26% of cases",
  [{ text: "H&E models learn whatever label they get: " }, { text: "which biological truth do they track?", options: { bold: true } }],
], { dot: IHC });
card(s, 0.4, 2.65, 4.55, 1.0, "Research question", [
  "When IHC, RNA and protein disagree, which one does an H&E foundation model follow, and does it depend on the slide type?",
], { fill: C.background2, dot: EOS, size: 11.5, plain: true });
card(s, 5.1, 1.0, 4.5, 2.65, "Five experiments", [
  [{ text: "1a  ", options: { bold: true, color: EOS } }, { text: "Is the IHC or the RNA label more predictable?" }],
  [{ text: "1b  ", options: { bold: true, color: EOS } }, { text: "On discordant cases, which side does the model take?" }],
  [{ text: "1c  ", options: { bold: true, color: EOS } }, { text: "Does it change from FFPE to frozen slides?" }],
  [{ text: "1d  ", options: { bold: true, color: EOS } }, { text: "Is HER2 0 vs 1+ visible in H&E at all?" }],
  [{ text: "1e  ", options: { bold: true, color: EOS } }, { text: "Do two frozen slides of one tumour agree?" }],
], { dot: H, size: 12 });
// pipeline flow
const steps = ["H&E slide\n(FFPE / frozen)", "Tissue\npatches", "Frozen encoder\n(UNI, CONCH, TITAN)", "Slide classifier\n(mean-pool, ABMIL)", "Compare labels\nsite-grouped CV"];
const bw = 1.62, gap = 0.2, by = 3.88;
steps.forEach((t, i) => {
  const x = 0.4 + i * (bw + gap);
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: by, w: bw, h: 0.68, rectRadius: 0.08, fill: { color: i === 4 ? H : C.background2 },
    line: { color: H, width: 1 }, objectName: nm("pipe_box") });
  s.addText(t, { x, y: by, w: bw, h: 0.68, fontSize: 10.5, bold: i === 4, color: i === 4 ? C.background1 : C.text2, align: "center",
    valign: "middle", margin: 2, isTextBox: true, objectName: nm("pipe_text") });
  if (i < 4) arrow(s, x + bw + 0.02, by + 0.34, x + bw + gap - 0.02, by + 0.34);
});
s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 5.86, y: 4.7, w: 2.8, h: 0.4, rectRadius: 0.08, fill: { color: C.background1 },
  line: { color: EOS, width: 1 }, objectName: "labels_box" });
s.addText([{ text: "Labels: " }, { text: "IHC", options: { color: IHC, bold: true } }, { text: "  |  " },
  { text: "RNA", options: { color: RNA, bold: true } }, { text: "  |  protein (RPPA)" }],
  { x: 5.86, y: 4.7, w: 2.8, h: 0.4, fontSize: 10.5, color: C.text1, align: "center", valign: "middle", margin: 0, isTextBox: true,
    objectName: "labels_text" });
arrow(s, 6.55, 4.7, 6.55, 4.58);
s.addNotes("Novelty is the evaluation, not a new architecture: same patients, slides and folds, evaluated against different reference standards.");

// ================= 3. Dataset
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Week 1: proposal and baselines" });
s.addText("Dataset: TCGA-BRCA", { placeholder: "title" });
pill(s, "WEEK 1", H);
stat(s, 0.45, 0.95, 1.6, "1,097", "patients (open access)", H);
stat(s, 2.15, 0.95, 1.6, "40", "tissue-source sites", H);
s.addText("Labels joined per patient: IHC, RNA (ESR1/PGR/ERBB2), PAM50, RPPA protein. Three slide types per patient:",
  { x: 0.45, y: 2.0, w: 3.9, h: 0.55, fontSize: 11.5, color: C.text1, valign: "top", margin: 0, isTextBox: true, objectName: "ds_intro" });
s.addTable([
  [hdr("Slide type"), hdr("Patients")],
  [lcell("FFPE diagnostic (DX)"), cell("1,060")],
  [lcell("Frozen top (TS)"), cell("1,043")],
  [lcell("Frozen bottom (BS)"), cell("427")],
], { x: 0.45, y: 2.62, w: 3.6, colW: [2.4, 1.2], rowH: 0.32, fontSize: 11, color: C.text1, border: { type: "solid", pt: 0.5, color: "D9D2E6" },
  fill: { color: C.background1 }, objectName: "ds_slides_table" });
s.addText("How often the “truths” disagree (no images)", { x: 4.55, y: 0.95, w: 5.0, h: 0.35, fontSize: 13, bold: true, color: C.text2,
  margin: 0, isTextBox: true, objectName: "ds_disc_title" });
s.addTable([
  [hdr("Marker"), hdr("Compared"), hdr("Discordant"), hdr("mRNA alone → IHC (AUC)")],
  [lcell("ER", { bold: true }), lcell("IHC vs RNA"), cell("5.5%"), cell("0.96")],
  [lcell("PR", { bold: true }), lcell("IHC vs RNA"), cell("16.7%"), cell("0.92")],
  [lcell("HER2", { bold: true }), lcell("IHC/FISH vs amplification"), cell("10.8%"), cell("0.87")],
], { x: 4.55, y: 1.35, w: 5.05, colW: [0.75, 1.95, 0.95, 1.4], rowH: 0.34, fontSize: 10.5, color: C.text1,
  border: { type: "solid", pt: 0.5, color: "D9D2E6" }, objectName: "ds_disc_table" });
card(s, 4.55, 2.95, 5.05, 1.75, "Notes", [
  "TCGA's HER2 field is noisy: we use a strict IHC/FISH rule",
  "ER discordance is mostly ER-low tumours (< 10% of cells stained)",
  [{ text: "All five experiments are feasible", options: { bold: true } }, { text: " on public data (smallest: 1d, 57 vs 261 cases)" }],
], { dot: RNA });
s.addNotes("Discordance is a minority of patients, which is why the paired, same-patient design matters.");

// ================= 4. Week 1 baselines
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Week 1: proposal and baselines" });
s.addText("Baselines and first results", { placeholder: "title" });
pill(s, "WEEK 1", H);
s.addText("3 frozen encoders on the DX slides of 1,055 patients → logistic regression on the slide embedding (TITAN) or the mean of all patches (UNI, CONCH); 5-fold site-grouped CV over seeds",
  { x: 0.4, y: 0.88, w: 9.2, h: 0.42, fontSize: 11, color: C.accent6, margin: 0, isTextBox: true, objectName: "w1_setup" });
s.addChart(pres.charts.BAR, [
  { name: "IHC label", labels: ["TITAN", "UNI v1", "CONCH v1"], values: [0.87, 0.87, 0.84] },
  { name: "RNA label", labels: ["TITAN", "UNI v1", "CONCH v1"], values: [0.92, 0.91, 0.87] },
], { x: 0.4, y: 1.35, w: 3.75, h: 2.55, barDir: "col", barGrouping: "clustered", chartColors: [HEX.accent4, HEX.accent3],
  showTitle: true, title: "1a  ER: AUC by reference label (n = 801)", titleFontSize: 11, titleColor: HEX.dk2, titleFontFace: "+mn-lt",
  showValue: true, dataLabelPosition: "outEnd", dataLabelFontSize: 9, dataLabelFormatCode: "0.00", dataLabelFontFace: "+mn-lt",
  valAxisMinVal: 0.7, valAxisMaxVal: 1.0, valAxisMajorUnit: 0.1, valAxisLabelFormatCode: "0.0", valAxisLabelFontSize: 9, catAxisLabelFontSize: 10,
  catAxisLabelColor: HEX.accent6, valAxisLabelColor: HEX.accent6, catAxisLabelFontFace: "+mn-lt", valAxisLabelFontFace: "+mn-lt",
  valGridLine: { color: "E6E1EE", size: 0.5 }, catGridLine: { style: "none" }, showLegend: true, legendPos: "b", legendFontSize: 9,
  legendFontFace: "+mn-lt", objectName: "w1_er_chart" });
s.addTable([
  [hdr(""), hdr("Task"), hdr("TITAN"), hdr("UNI v1"), hdr("CONCH v1")],
  [lcell("1a AUC", { bold: true }), lcell("PR: IHC → RNA (n=986)"), cell("0.80 → 0.79"), cell("0.80 → 0.73"), cell("0.78 → 0.75")],
  [lcell(""), lcell("HER2: clin. → amp. (n=478)"), cell("0.68 → 0.72"), cell("0.63 → 0.68"), cell("0.62 → 0.68")],
  [lcell("1b side", { bold: true }), lcell("ER discordant (n=54)"), cell("39:15**", { color: RNA, bold: true }), cell("34:20"), cell("35:19*", { color: RNA, bold: true })],
  [lcell(""), lcell("PR discordant (n=169)"), cell("69:100*", { color: IHC, bold: true }), cell("55:114***", { color: IHC, bold: true }), cell("61:108***", { color: IHC, bold: true })],
  [lcell(""), lcell("HER2 discordant (n=60)"), cell("35:25"), cell("44:16***", { color: RNA, bold: true }), cell("39:21*", { color: RNA, bold: true })],
  [lcell("Control", { bold: true }), lcell("Shuffled ER labels"), cell("0.46"), cell("0.51"), cell("0.46")],
  [lcell(""), lcell("HER2 0 vs low"), cell("0.40"), cell("0.50"), cell("0.43")],
], { x: 4.3, y: 1.35, w: 5.3, colW: [0.72, 1.88, 0.9, 0.9, 0.9], rowH: 0.29, fontSize: 9.5, color: C.text1,
  border: { type: "solid", pt: 0.5, color: "D9D2E6" }, objectName: "w1_table" });
s.addText([{ text: "1b: discordant cases siding molecular : IHC; " }, { text: "teal", options: { color: RNA, bold: true } },
  { text: " = molecular side wins, " }, { text: "amber", options: { color: IHC, bold: true } }, { text: " = IHC side wins (binomial *, **, *** = p < 0.05, 0.01, 0.001). HER2-low: mRNA alone 0.69. TITAN 5 seeds, UNI/CONCH 3." }],
  { x: 4.3, y: 3.72, w: 5.3, h: 0.4, fontSize: 8.5, color: C.accent6, margin: 0, isTextBox: true, objectName: "w1_footnote" });
card(s, 0.4, 4.18, 4.55, 0.98, "Consistent across all three encoders", [
  [{ text: "ER, HER2: ", options: { bold: true } }, { text: "H&E tracks the molecular label more than IHC" }],
  [{ text: "PR: ", options: { bold: true } }, { text: "the opposite, sides with IHC (to be checked)" }],
], { dot: RNA, size: 11 });
card(s, 5.1, 4.18, 4.5, 0.98, "Site leakage and HER2-low", [
  "Embeddings identify the hospital 82% of the time",
  "HER2 0 vs 1+: 0.63 (random folds) → chance (site-grouped)",
], { dot: EOS, size: 11 });
s.addNotes("These are the Week 1 numbers as first computed. The PR result is revisited in Week 2.");

// ================= 5. Week 2: ABMIL + PR re-check
pres.addSection({ title: "Week 2: updates" });
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Week 2: updates" });
s.addText("ABMIL agrees; PR was a cut-off artefact", { placeholder: "title" });
pill(s, "WEEK 2", EOS);
s.addText("Attention-MIL (ABMIL) on UNI and CONCH patches, 3 seeds, same site-grouped folds", { x: 0.4, y: 0.88, w: 4.5, h: 0.3,
  fontSize: 11, color: C.accent6, margin: 0, isTextBox: true, objectName: "w2a_setup" });
s.addTable([
  [hdr("ABMIL"), hdr("UNI v1"), hdr("CONCH v1")],
  [lcell("1a ER: IHC → RNA"), cell("0.86 → 0.91", { bold: true, color: RNA }), cell("0.85 → 0.90", { bold: true, color: RNA })],
  [lcell("1a HER2: clinical → amplif."), cell("0.67 → 0.69"), cell("0.67 → 0.73")],
  [lcell("1b ER discordant (n=54)"), cell("34:20"), cell("36:18*", { bold: true, color: RNA })],
  [lcell("1b HER2 discordant (n=60)"), cell("47:13***", { bold: true, color: RNA }), cell("34:26")],
  [lcell("Shuffled ER (control)"), cell("0.49"), cell("0.46")],
  [lcell("HER2 0 vs low"), cell("0.55"), cell("0.52")],
], { x: 0.4, y: 1.22, w: 4.5, colW: [2.1, 1.2, 1.2], rowH: 0.31, fontSize: 10.5, color: C.text1,
  border: { type: "solid", pt: 0.5, color: "D9D2E6" }, objectName: "w2a_abmil_table" });
card(s, 0.4, 3.55, 4.5, 1.4, "Takeaway", [
  "ABMIL agrees with mean-pooling: the RNA-defined ER label is more predictable, and HER2 calls side with amplification",
  "HER2-low stays at chance with attention too",
], { dot: RNA, size: 11 });
s.addChart(pres.charts.BAR, [
  { name: "GMM RNA cut-off", labels: ["TITAN", "UNI mean", "CONCH mean", "UNI ABMIL", "CONCH ABMIL"], values: [-0.006, -0.047, -0.027, -0.058, -0.029] },
  { name: "Prevalence-matched cut-off", labels: ["TITAN", "UNI mean", "CONCH mean", "UNI ABMIL", "CONCH ABMIL"], values: [0.033, -0.008, 0.010, -0.006, 0.011] },
], { x: 5.1, y: 0.95, w: 4.5, h: 2.75, barDir: "col", barGrouping: "clustered", chartColors: [HEX.accent5, HEX.accent3],
  showTitle: true, title: "PR: ΔAUC (RNA − IHC label) by RNA cut-off", titleFontSize: 11, titleColor: HEX.dk2, titleFontFace: "+mn-lt",
  showValue: true, dataLabelPosition: "outEnd", dataLabelFontSize: 8, dataLabelFormatCode: "+0.00;-0.00;0.00", dataLabelFontFace: "+mn-lt",
  valAxisMinVal: -0.08, valAxisMaxVal: 0.06, valAxisLabelFormatCode: "+0.00;-0.00;0", valAxisLabelFontSize: 9, catAxisLabelFontSize: 9,
  catAxisLabelColor: HEX.accent6, valAxisLabelColor: HEX.accent6, catAxisLabelFontFace: "+mn-lt", valAxisLabelFontFace: "+mn-lt",
  valGridLine: { color: "E6E1EE", size: 0.5 }, catGridLine: { style: "none" }, showLegend: true, legendPos: "b", legendFontSize: 9,
  legendFontFace: "+mn-lt", objectName: "w2a_pr_chart" });
card(s, 5.1, 3.82, 4.5, 1.28, "PR: no reliable preference", [
  "GMM cut-off labels more tumours PR-negative than IHC does; matching IHC prevalence removes the gap",
  "Discordant PR (n = 133) then splits evenly: all five models p > 0.3",
], { dot: IHC, size: 10.5 });
s.addNotes("Week 1 plan items: finish ABMIL; test the PR RNA cut-off before interpreting the reversal. The apparent PR 'IHC preference' disappears with a prevalence-matched cut-off in all five models.");

// ================= 6. Week 2: ER discordant + attention
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Week 2: updates" });
s.addText("ER-low tumours are read as ER-negative", { placeholder: "title" });
pill(s, "WEEK 2", EOS);
s.addText("Model trained on concordant patients; which label does it agree with on discordant ones? (RNA : IHC)",
  { x: 0.4, y: 0.88, w: 4.75, h: 0.3, fontSize: 10.5, color: C.accent6, margin: 0, isTextBox: true, objectName: "w2b_setup" });
s.addTable([
  [hdr("ER discordant subgroup"), hdr("n"), hdr("TITAN"), hdr("UNI"), hdr("CONCH")],
  [lcell("IHC+ / RNA−"), cell("34"), cell("24:10"), cell("21:13"), cell("22:12")],
  [lcell("IHC− / RNA+"), cell("20"), cell("15:5"), cell("13:7"), cell("13:7")],
  [lcell("ER-low (IHC+, < 10% stained)", { bold: true }), cell("11"), cell("9:2", { bold: true, color: RNA }), cell("9:2", { bold: true, color: RNA }), cell("9:2", { bold: true, color: RNA })],
  [lcell("ER 10–39% stained"), cell("9"), cell("7:2"), cell("6:3"), cell("7:2")],
], { x: 0.4, y: 1.22, w: 4.75, colW: [2.0, 0.45, 0.75, 0.75, 0.8], rowH: 0.31, fontSize: 10.5, color: C.text1,
  border: { type: "solid", pt: 0.5, color: "D9D2E6" }, objectName: "w2b_table" });
s.addText("ABMIL attention maps: 54 discordant + 10 reference patients (UNI)", { x: 0.4, y: 2.95, w: 4.75, h: 0.3, fontSize: 12, bold: true,
  color: C.text2, margin: 0, isTextBox: true, objectName: "w2b_att_title" });
stat(s, 0.4, 3.22, 1.5, "35 / 54", "discordant: side with RNA", RNA);
stat(s, 1.95, 3.22, 1.5, "8 / 10", "ER-low: side with RNA", RNA);
stat(s, 3.5, 3.22, 1.6, "r = 0.52", "P(ER+) vs ESR1 level", H);
s.addText([
  { text: "References correctly classified 8 / 10; attention is focused (5–7% of patches hold half of it) and similar across groups", options: { bullet: { indent: 12 }, breakLine: true } },
  { text: "Both ER-low exceptions are PAM50 HER2-enriched; all 6 Basal ER-low tumours are called ER−", options: { bullet: { indent: 12 } } },
], { x: 0.4, y: 4.18, w: 4.75, h: 0.98, fontSize: 10.5, color: C.text1, valign: "top", paraSpaceAfter: 3, margin: 0, isTextBox: true,
  objectName: "w2b_att_notes" });
// heatmaps
s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 5.35, y: 0.95, w: 4.25, h: 4.2, rectRadius: 0.08, fill: { color: C.background2 },
  line: { color: C.background2 }, objectName: "w2b_fig_bg" });
s.addImage({ path: path.join(DIR, "att_TCGA-A2-A4S1.png"), x: 5.55, y: 1.12, w: 3.85, h: 3.85 * 271 / 524, objectName: "w2b_heat1" });
s.addText("ER-low (IHC+ < 10%), RNA−, Basal → called ER−, P(ER+) = 0.00", { x: 5.5, y: 3.2, w: 3.95, h: 0.3, fontSize: 9.5,
  color: C.text1, align: "left", margin: 0, isTextBox: true, objectName: "w2b_cap1" });
s.addImage({ path: path.join(DIR, "att_TCGA-AR-A0TX.png"), x: 5.55, y: 3.62, w: 1.6, h: 1.6 * 438 / 505, objectName: "w2b_heat2" });
s.addText([{ text: "Exception: ER-low, RNA−, HER2-enriched → called ER+ (with IHC)", options: { breakLine: true } },
  { text: "Colour = attention percentile within the slide (bright = high); patch positions only", options: { color: C.accent6, fontSize: 8.5 } }],
  { x: 7.3, y: 3.66, w: 2.2, h: 1.4, fontSize: 9.5, color: C.text1, valign: "top", margin: 0, isTextBox: true, objectName: "w2b_cap2" });
s.addNotes("Week 1 plan items: ER-low vs ER-high discordant cases; attention maps on discordant cases. The ER RNA preference is driven by ER-low tumours, which the model reads as ER-negative; HER2-enriched tumours are the exception.");

// ================= 7. Week 2: fixed splits + frozen pipeline
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Week 2: updates" });
s.addText("Fixed splits; frozen slides encoded", { placeholder: "title" });
pill(s, "WEEK 2", EOS);
card(s, 0.4, 0.95, 3.1, 2.15, "Fixed site-grouped CV splits", [
  "Two scikit-learn versions gave different folds for the same seed",
  "One shared split file (data/splits/folds.csv) for every script; library pinned",
  "All linear baselines rerun on Ada: conclusions unchanged",
], { dot: H, size: 10.5 });
stat(s, 0.55, 3.2, 2.9, "+0.058", "TITAN ER ΔAUC on fixed splits (95% CI +0.031, +0.084)", RNA);
// frozen pipeline
s.addText("Frozen-slide feature pipeline (Ada, RTX 2080 Ti)", { x: 3.75, y: 0.95, w: 5.85, h: 0.3, fontSize: 13, bold: true, color: C.text2,
  margin: 0, isTextBox: true, objectName: "w2c_pipe_title" });
const fs = ["Manifest\n1,516 TS/BS slides\n373 TS+BS pairs", "Download\nIDC copy\n(Google Cloud)", "TRIDENT\ntissue seg. +\n256 px @ 20×", "UNI v1\n(frozen)", "Mean + 4,096-\npatch bag;\ndelete slide"];
const fw = 1.02, fg = 0.19, fy = 1.33;
fs.forEach((t, i) => {
  const x = 3.75 + i * (fw + fg);
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: fy, w: fw, h: 0.95, rectRadius: 0.08, fill: { color: i === 1 ? EOS : C.background2 },
    line: { color: i === 1 ? EOS : H, width: 1 }, objectName: nm("fpipe_box") });
  s.addText(t, { x, y: fy, w: fw, h: 0.95, fontSize: 9, color: i === 1 ? C.background1 : C.text2, bold: i === 1, align: "center",
    valign: "middle", margin: 1, isTextBox: true, objectName: nm("fpipe_text") });
  if (i < 4) arrow(s, x + fw + 0.01, fy + 0.475, x + fw + fg - 0.01, fy + 0.475);
});
stat(s, 3.75, 2.45, 1.85, "1,515", "frozen slides encoded (1,044 TS, 471 BS)", H);
stat(s, 5.75, 2.45, 1.85, "~75 s", "per slide; 3 GPUs in parallel", H);
stat(s, 7.75, 2.45, 1.85, "1.0000", "cosine to TRACE features (46 DX slides re-extracted)", RNA);
card(s, 3.75, 3.6, 2.85, 1.55, "Why IDC, not GDC", [
  "GDC from Ada: 0.02–0.3 MB/s per connection (430 GB → weeks)",
  "Same slides in IDC (DICOM) via Google Cloud: ~10 MB/s",
], { dot: EOS, size: 10 });
card(s, 6.75, 3.6, 2.85, 1.55, "Ready for 1c / 1e", [
  "Our pipeline reproduces the DX features exactly, so DX vs frozen differ only in tissue",
  "UNI weights access granted; all slides done",
], { dot: RNA, size: 10 });
s.addNotes("Week 1 logistics: fixed CV splits; UNI weights access; download + encode the frozen slides. Week 1 plan: start frozen slides for 1c and 1e. Features are ready; the 1c/1e analyses are Week 3.");

// ================= 8. Plan for Week 3
pres.addSection({ title: "Week 3 plan" });
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Week 3 plan" });
s.addText("Plan for the next week (Week 3)", { placeholder: "title" });
pill(s, "WEEK 3", RNA);
card(s, 0.4, 0.95, 2.95, 2.55, "Done in Week 2", [
  "✓  ABMIL (UNI, CONCH) for 1a, 1b, 1d",
  "✓  PR cut-off test: no reliable preference",
  "✓  ER-low vs ER-high discordant cases",
  "✓  Attention maps on discordant cases",
  "✓  Fixed CV splits; baselines rerun",
  "✓  UNI access; 1,515 frozen slides encoded",
], { dot: H, size: 10.5, fill: C.background2, plain: true });
card(s, 3.52, 0.95, 2.95, 2.55, "Next experiments", [
  [{ text: "1c: ", options: { bold: true } }, { text: "FFPE vs frozen × IHC vs RNA label; interaction test; % tumour cells as covariate" }],
  [{ text: "1e: ", options: { bold: true } }, { text: "TS vs BS agreement on 373 pairs (ICC, κ, % flipped calls)" }],
  [{ text: "Statistics: ", options: { bold: true } }, { text: "DeLong tests for 1a; combined 1b evidence across encoders" }],
], { dot: EOS, size: 10.5 });
card(s, 6.65, 0.95, 2.95, 2.55, "Analysis and logistics", [
  "Overlay attention on the real slide images (does it sit on tumour?)",
  "Describe discordant patients: PAM50, % of cells ER-stained, histology",
  "Rerun ABMIL on the fixed splits for one set of final tables",
], { dot: RNA, size: 10.5 });
s.addTable([
  [hdr("Compute on Ada (11 GB GPUs)"), hdr("Estimate")],
  [lcell("Frozen-slide features (1,515 slides)"), cell("done, ~15 h on 2–3 GPUs")],
  [lcell("1c / 1e analyses (linear probes, bootstrap)"), cell("< 1 h CPU")],
  [lcell("ABMIL rerun on fixed splits (UNI + CONCH)"), cell("~10 GPU-hours")],
], { x: 1.7, y: 3.75, w: 6.6, colW: [4.2, 2.4], rowH: 0.29, fontSize: 10.5, color: C.text1,
  border: { type: "solid", pt: 0.5, color: "D9D2E6" }, objectName: "w3_compute" });
s.addNotes("Everything listed under Next is Week 3 work; Week 2 delivered all items planned in Week 1.");

(async () => {
  await pres.writeFile({ fileName: OUT });
  await applyTheme(OUT, THEME);
  console.log("wrote", OUT);
})();

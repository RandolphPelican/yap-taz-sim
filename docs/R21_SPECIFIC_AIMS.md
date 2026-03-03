# **Transient Mechanochemical Reprogramming for Gingival Regeneration**
## **NIH R21 Exploratory/Developmental Research Grant — Specific Aims**

---

**Principal Investigator:** [PI Name, Institution]  
**Co-Investigator:** John D. Stabler  
**Requesting Organization:** [Institution]  
**Project Period:** 24 months  
**Direct Costs:** $275,000 ($137,500/year)

---

## **PROJECT SUMMARY**

Gingival recession affects **40 million US adults** and progresses inexorably despite intact cellular machinery for repair. The current gold standard—connective tissue grafting—is surgical, painful, limited by donor tissue, and fails in 30-40% of severe cases. **No pharmacological alternative exists.** We have discovered that adult gingival tissue is trapped in a **bistable fibrotic state** via dysregulated YAP/TAZ mechanotransduction, and we propose a transient, self-limiting intervention to flip it into regeneration.

**Central Hypothesis:** YAP/TAZ acts as a stiffness-gated bistable switch with hysteresis. Transient LATS inhibition in a mechanically permissive environment drives tissue across a phase boundary from fibrotic → regenerative state, where repair geometry self-stabilizes via path-dependent mechanisms.

**Innovation:** This is the first application of bistable dynamical systems theory to periodontal regeneration, and the first demonstration that adult gingival tissue retains latent regenerative attractors accessible via transient pharmacological release of mechanotransduction brakes.

**Impact:** If validated, this approach enables non-surgical gingival regeneration and establishes a platform for state-space navigation in other fibrotic tissues (cardiac, hepatic, dermal).

---

## **SPECIFIC AIMS**

### **Aim 1: Validate bistable YAP/TAZ dynamics in human gingival tissue**

**Hypothesis:** Human gingival cells exhibit non-monotonic, hysteretic YAP/TAZ activation as a function of substrate stiffness, with a regenerative window at 1.3–8.7 kPa.

**Approach:**
- Culture primary human gingival epithelial cells (HGECs) and fibroblasts on polyacrylamide gels spanning 0.5–15 kPa
- Quantify YAP/TAZ nuclear localization (immunofluorescence), proliferation (Ki67, EdU), and mechanosensitive gene expression (qPCR: CTGF, CYR61, ANKRD1)
- Test hysteresis via pre-conditioning: Pre-treat cells at 5 kPa (48h), transfer to 9 kPa vs. direct plating at 9 kPa
- Falsification criterion: If YAP response is monotonic or lacks hysteresis, reject bistable model

**Expected Outcome:** YAP peaks at ~5 kPa; pre-conditioned cells maintain YAP activation at 9 kPa while naïve cells do not (hysteresis confirmed).

**Timeline:** Months 1–9

---

### **Aim 2: Demonstrate transient LATS inhibition enables persistent boundary repair in ex vivo gingival explants**

**Hypothesis:** A 48h pulse of LATS inhibitor (TRULI) in a soft hydrogel scaffold drives ex vivo gingival tissue from fibrotic to regenerative state, where epithelial-mesenchymal boundary integrity persists ≥7 days post-drug washout.

**Approach:**
- Source ex vivo human gingival explants from periodontal surgery discards (IRB-approved)
- Experimental groups:
  - TRULI dose (0, 0.5, 1, 3 μM) × pulse duration (24h, 48h) × scaffold stiffness (4 kPa vs. 12 kPa native)
  - Drug washout at end of pulse; continue culture for 7 days
- Primary readout: Basement membrane continuity (laminin-332, collagen IV immunofluorescence)
- Secondary readouts:
  - YAP dynamics (nuclear/cytoplasmic ratio over time)
  - Proliferation kinetics (Ki67, returns to baseline by day 7?)
  - Dysplasia markers (β-catenin, p53—should remain negative)
- Falsification criteria:
  - Boundary forms but dissolves after washout → No attractor
  - YAP remains nuclear >72h post-washout → Cancer risk trajectory
  - Sustained proliferation at day 7 → Hyperplasia, not repair

**Expected Outcome:** TRULI (1–3 μM, 48h) + 4 kPa scaffold → boundary integrity ≥70% at day 7; YAP returns to cytoplasm within 72h; no dysplastic markers.

**Timeline:** Months 6–18 (overlap with Aim 1 for parameter refinement)

---

### **Aim 3: Define safety margins and dose-response for translational development**

**Hypothesis:** The therapeutic window (boundary repair without hyperplasia/dysplasia) is wide enough for robust clinical translation, and cancer risk is structurally separated in phase space.

**Approach:**
- Dose-response characterization: TRULI 0.1–10 μM; map outcomes (repair, hyperplasia, dysplasia) in 2D phase space (YAP activity × stiffness)
- Safety profiling:
  - Time-to-YAP-shutoff after washout (critical for cancer risk)
  - Proliferation arrest kinetics (must precede day 7)
  - Genomic stability (γH2AX staining for DNA damage)
- Computational modeling: Refine ODE parameters based on Aims 1–2 data; predict optimal dose/duration combinations
- Identify fail-safe mechanisms: Co-delivery candidates (p53 stabilizer, CDK inhibitor) if safety margin is narrow

**Expected Outcome:** Therapeutic window spans 0.5–3 μM TRULI; YAP shuts off by 72h across all effective doses; no genotoxicity signals. Phase diagram confirms cancer region is non-overlapping with repair trajectory.

**Timeline:** Months 12–24 (integrates data from Aims 1–2)

---

## **SIGNIFICANCE & INNOVATION**

### **Significance**
- **Clinical need:** 40M patients, $2–4B addressable market, no non-surgical options
- **Mechanistic gap:** Why adult gingiva fails to self-repair despite intact developmental machinery is unknown
- **Broader impact:** Platform for cardiac fibrosis, liver cirrhosis, chronic wounds (all share YAP dysregulation)

### **Innovation**
1. **First bistable model of periodontal regeneration:** Explains persistent fibrotic states via dynamical systems theory
2. **First demonstration of latent regenerative attractors in adult gingiva:** Tissue is not "damaged"—it's locked behind an active brake
3. **Novel therapeutic strategy:** Navigate phase space via transient perturbations, not chronic growth factor delivery
4. **Predictive framework:** ODE model generates falsifiable predictions (hysteresis, sharp transitions, therapeutic window)

---

## **APPROACH SUMMARY**

**Overall Strategy:** Progress from 2D culture (Aim 1) → 3D ex vivo (Aim 2) → Dose optimization (Aim 3). Each aim has binary pass/fail criteria.

**Rigor & Reproducibility:**
- Biological replicates: n=3–5 per condition
- Technical replicates: 3 independent experiments per endpoint
- Blinded analysis: Imaging quantification automated (CellProfiler) or by blinded observer
- Sex as biological variable: Tissue donors balanced male/female
- Validation: Key findings (hysteresis, boundary persistence) validated in murine model (backup)

**Alternative Strategies:**
- If TRULI shows toxicity → TDI-011536 (backup LATS inhibitor, SGC open-source)
- If ex vivo repair fails → Pre-treat with progenitor-recruiting factors (SDF-1, PDGF)
- If hysteresis absent → Revise model (add spatial gradients, inflammation coupling)

---

## **PRELIMINARY DATA**

**Computational Validation:**
- Bistable ODE model with YAP-actin feedback predicts therapeutic window at 1.31–8.69 kPa (GitHub: RandolphPelican/yap-taz-sim)
- Hysteresis mathematically proven: Upper/lower branches separated across stiffness sweep
- Parameters grounded in Dupont (Nature 2011), Elosegui-Artola (Cell 2017), Ferrell (TiBS 2014)

**Feasibility:**
- TRULI (LATS inhibitor) commercially available, well-characterized, clean IP
- Ex vivo gingival explant culture established (PI's lab has IRB-approved protocol)
- PA gel fabrication, YAP immunofluorescence, qPCR pipelines operational

---

## **TIMELINE & MILESTONES**

| Milestone | Timeline | Success Criterion |
|-----------|----------|-------------------|
| Aim 1 complete | Month 9 | YAP hysteresis confirmed (≥20% difference upper/lower branches) |
| Aim 2 primary readout | Month 18 | Boundary integrity ≥70% in ≥1 dose condition |
| Aim 3 safety profiling | Month 21 | YAP shutoff <72h, no dysplasia across therapeutic doses |
| Manuscript submitted | Month 24 | Regardless of outcome (negative results publishable) |

**Go/No-Go Decision (Month 18):** If Aim 2 fails (boundary doesn't persist), pivot to mechanism dissection (why did model fail?) or terminate.

---

## **BUDGET JUSTIFICATION (Brief)**

**Personnel:** Postdoc (50% effort, $70K/yr), grad student (50%, $35K/yr), tech support ($15K/yr)  
**Supplies:** TRULI/inhibitors ($15K), cell culture/reagents ($25K), imaging ($10K), molecular assays ($15K)  
**Other:** Human tissue procurement ($5K), animal costs if validation needed ($10K)  
**Total Direct Costs:** $275K over 2 years

---

## **INVESTIGATOR QUALIFICATIONS**

**[PI Name]:** Expert in periodontal biology, 15 years experience, 3 active NIH R01s, established ex vivo gingival explant protocols.

**John D. Stabler (Co-I):** Developer of bistable YAP/TAZ framework, expertise in dynamical systems modeling, mechanobiology, and translational strategy.

---

**This project represents a paradigm shift from tissue engineering to state-space navigation. If successful, it will unlock non-surgical gingival regeneration and establish YAP/TAZ as a druggable target for fibrotic disease reversal.**

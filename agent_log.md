# Photon-hunter agent log

Remote session(s) by the photon-hunter agent. Newest entries at the bottom. Data:
`Desktop/HOI_freshstaret/data/raw/troubleshoot/`, figures: `.../data/processed/troubleshoot_figures/`.

## 2026-10-20 session 1 (remote; Tristan: 45 kHz triggers idler open, 5 kHz idler blocked; spot (501, 454))

### E0 (14:55, no camera) Re-analysis of hunt_20261020_134544 at the known spot
- Hypothesis: the hunt scored its coarse windows at an auto-found "spot" (486, 430), not at the beam
  position (501, 454), so it could have missed a heralded excess at the real spot.
- Method: re-scored every saved wide-gate and coarse run at both positions (aperture r=10 px,
  ring 20-35 px; hits = pixels >= 200 ADU above the readout's median).
- Result: at (501, 454) the coarse 1 us windows (27 ns .. 10 us, 160 readouts x 5000 acc) are flat:
  ring-subtracted spot ADU 1540-1760 per readout (+-45), spot hits 0.12-0.30 per readout, with no
  window standing out. Wide gate at (501, 454): 98148 ADU, contrast 44.6 at 12:21 (high pump);
  9446 ADU, contrast 4.2 at 13:24; 4480 ADU, contrast 1.4 at 13:46 (the hunt reference).
  The ~1700 ADU at 1 us is the same at every delay and in the 27 ns window, so it is fixed pattern,
  not light.
- Conclusion: the wrong spot did not hide anything; the old coarse scan is flat at the real spot too.
  The wide-gate beam at (501, 454) was already fading between 13:24 and 13:46.

### E1 (14:58) Trigger rate vs threshold and edge; wide-gate beam check (agent_scripts/exp1_trigger_and_beam.py)
- Hypothesis: if the herald pulse is only ~0.5 V, the rate should fall steeply above 0.5 V and climb
  below it (noise). A flat rate means clean pulses.
- Settings: gain 1, 3 ns gate, 100,000 acc x 4 readouts per point, rising edge 0.2-1.0 V, falling 0.5 V.
- Result: 38.3, 37.7, 37.7, 38.2, 39.6, 41.3, 41.7, 43.0 kHz at 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0 V
  (rising); 43.4 kHz falling 0.5 V. Flat within the slow drift of the source (the readings go up in
  time order, not with threshold).
- Wide gate 309 ns + 100 us, 5000 acc x 40 readouts, gain 9000: ring-subtracted spot at (501, 454)
  694 +- 408 ADU per readout (1.7 sigma), compactness contrast -0.1: no compact beam spot.
  Figure: troubleshoot_figures/agent_exp1_agent_wide100us.png.
- Incident: the script did not set the Agg backend, so plt.show() blocked after the first wide
  run; I stopped my own process and re-opened the camera in E2 (which closes it properly).

### E2 (15:15) Thresholds above 1 V; long wide gate at the light limit (agent_scripts/exp2_threshold_high_and_long_wide.py)
- Result, trigger rate: 41.8, 40.9, 39.9, 39.0 kHz at 1.5, 2.0, 2.5, 3.0 V, then 38.4 kHz at 0.5 V.
  The rate is flat from 0.2 V to 3.0 V: the camera sees clean pulses taller than 3 V (the "0.5 V"
  scope reading cannot be the pulse that reaches TRIGGER IN, or the pulse is much taller than
  measured). Marginal triggering and noise triggers are ruled out; threshold 0.5 V is fine.
- Wide gate 309 ns + 20 us (shorter than the ~22 us mean trigger spacing), 50,000 acc
  (1.0 s gate-open per readout) x 100 readouts, gain 9000: spot (501, 454) ring-subtracted
  -6124 +- 513 ADU per readout; aperture hits 30.9 per readout. The image shows a ~50 px wide
  ring centred near (500, 450) with a dark middle, and a bright blob near (750, 440).
  Figure: troubleshoot_figures/agent_exp2_wide20us_1s.png.
- Analysis (radial profiles, all wide-gate runs, per 1 s of gate-open time): the ring at r = 20-40 px
  and the blob at (750, 440) are equally present in the 12:25 signal-arm-BLOCKED run, so they are
  not signal-arm light. The beam itself is the compact core (r < 6 px minus r = 10-16 px):
  12:21 high pump 123,000 +- 2,800; blocked -2,700 +- 500; 13:24 13,500 +- 1,100;
  13:46 6,900 +- 1,100; now 2,250 +- 540 (E1) and 790 +- 290 (E2) ADU per s gate-open.
  Relative to the blocked run the beam is still there now at ~3,500-5,000 ADU/s (6-7 sigma),
  i.e. ~25-35x weaker than at 12:21 and ~2x weaker than during the last hunt.
- Conclusion: triggers are clean (~40 kHz) but very little signal-arm light reaches the camera:
  of order a few detected photons per second. Heralded photons per trigger are at most
  (signal detections/s) / (triggers/s) ~ 1e-4 even if every signal photon were heralded, so
  1 us windows (background ~3e-5 hit pixels per gate per us in the spot) need ~1000s of readouts.
  Narrow (100 ns) windows with 50,000 accumulations per readout give ~10x better S/B at the
  same trigger time; next step is a 100 ns-window scan over the fibre-delay estimate region.
- Background budget used for planning (old runs): spot hits per readout 0.05 with 10 ns x 5000 acc,
  ~0.2 with 1 us x 5000 acc; ~31 per 1 s gate-open (E2). ROI-wide hits are mostly per readout.

### E3 (15:22-15:46) 100 ns windows, 27-2027 ns, interleaved (agent_scripts/window_scan.py --tag exp3)
- Hypothesis: with the beam this weak, 1 us windows are background-limited; 100 ns windows with
  50,000 acc per readout should cut the in-gate background ~10x at the same trigger time.
- Settings: gain 9000, 20 windows of 100 ns (27, 127, ... 1927 ns), 50,000 acc, 3 passes x 15
  readouts in shuffled order -> 45 readouts (2.25 M gates, ~60 s of triggers) per window; spot (501, 454).
- Raw result (spot hits r<=10, >=200 ADU above the readout median): 3.4-5.8 hits per readout in every
  window; best 1827 ns at +3.9 sigma, earliest windows low (27 ns: -3.1, r<=20: -6.8 sigma).
- New finding, background: hits per readout track the readout interval (corr 0.86 over exp3 runs;
  ~200 ROI hits/s whatever the gate width: 100 ns x 50k and 10 ns x 5k give the same hits per
  second). So the dominant background is NOT in the gate: it accumulates with wall-clock time
  (MCP dark events or light leaking past the gated photocathode, plus CCD). The trigger rate
  drifts (readout 1.16-1.47 s per 50k triggers), which moved the raw scores; the early-window
  deficit is drift, the same on the whole sensor (also seen at r > 50 px and in the pedestal).
  Consequence: narrow windows do not lower the background much; wider windows cost little.
- Drift-corrected score (spot hits minus area-scaled hits at r > 50 px in the same readout):
  flat. All windows between -0.66 and +0.77 net hits per readout (errors ~0.3); best 1827 ns at
  +2.9 sigma (r<=10) / +1.9 sigma (r<=6), driven by one pass (+0.41, +0.37, +1.52). Not significant
  for 20 trials.
- Sensitivity: beam now ~12-15 hit pixels/s (gate open) in r<=10, so heralded <= 15 x eta hits/s,
  eta = probability that a detected signal photon's twin made a trigger. ~60 s per window at
  ~3 hits/s background excludes eta >~ 0.1 for a delay in 27-2027 ns.
- Event heights: photon events at gain 9000 peak at only ~200-300 ADU per pixel (median ~210-270),
  so the 200 ADU threshold cuts ~40% of real events; but background events have similar heights
  (median ~190), so raising/lowering the threshold does not separate them.
- Figures: troubleshoot_figures/agent_exp3_20261020_152224/ (window_scan.png, delay_scan_*.png).
- Code: window_scan.py now scores with the drift-corrected net hits (tested offline on exp3 runs:
  reproduces 1827 ns +0.77, z 2.9).

### E4 (15:47-16:14) 500 ns windows, 27-5027 ns, interleaved (window_scan.py --tag exp4)
- Hypothesis: since the background is mostly ungated (E3), wider windows cost little; 500 ns
  windows cover 5 us with the same per-window sensitivity as E3.
- Settings: gain 9000, 10 windows of 500 ns (27 ... 4527 ns), 50,000 acc (25 ms gate-open per
  readout), 6 passes x 12 readouts shuffled -> 72 readouts per window.
- Result (drift-free net spot hits per readout, r<=10): -1.40 .. -0.21 for all windows; best 527 ns
  at +1.3 sigma vs the others (527 ns was best in passes 1, 5, 6 at +1.3..+1.6; never > 2 sigma).
  Raw hits still rise with delay (7.5 at 27 ns -> 9.8 at 4527 ns) and r<=20 shows 4527 at +3.8:
  that is the dead-time effect below, not light.
- Camera-side check (free from the timing): readout interval = 1.70 s + 1.17 x (delay x 50k).
  Slope ~1 means every gate blocks triggers for about its delay, i.e. the camera does open one
  gate per trigger at the programmed delay. Longer delay -> fewer gates/s -> longer readout ->
  more ungated background per readout. This is why raw scores rise with delay.
- Source drift: the trigger rate fell through the session: ~38-43 kHz (14:58-15:15), ~36-38 kHz
  (15:23), ~29 kHz (15:50) and ~22-26 kHz (16:10) (from 50k acc per readout timing, dead time
  removed). Together with the beam fading 30x since 12:21 this suggests the pump power or the
  SPDC/idler alignment is drifting down, not only Tristan's pump change.
- Conclusion: no heralded window between 27 ns and 5027 ns at this sensitivity.
- Figures: troubleshoot_figures/agent_exp4_20261020_154719/.

### E5 (16:14-16:39) On/off-delay A/B at the fibre-delay estimate (agent_scripts/ab_test.py --tag exp5)
- Hypothesis: if heralded photons arrive 450-600 ns after the trigger, a 300-900 ns gate shows more
  net spot hits than the same gate 10 us later.
- Settings: gain 9000, gate 600 ns, 50,000 acc (30 ms gate-open per readout), on = 300 ns,
  off = 10,300 ns, alternated in blocks of 12 readouts, 23 cycles -> 276 readouts each
  (13.8 M gates each; 2.28 s per on readout -> ~22 kHz accepted triggers).
- Result: drift-free net spot hits -0.77 (on) vs -0.91 (off) per readout, z = +0.5. Raw hits
  10.8 vs 14.3 (z = -12): the off runs read out slower (2.78 vs 2.28 s, dead time of the 10 us
  delay), so they collect more ungated background; the net metric removes that.
- Upper limit: on - off < ~0.7 hit pixels per readout (2 sigma) -> < ~1.4e-5 heralded hit pixels
  per trigger in r <= 10 px at 300-900 ns.
- Figure: troubleshoot_figures/agent_exp5_20261020_161415/ab_test.png.

### Session summary (camera ~100 min, closed with ht.shutdown after every script)
- No heralded photons found. Nothing reached 5 sigma: best candidates 1827 ns (+2.9, E3),
  527 ns (+1.3, E4), 300-900 ns on/off +0.5 (E5).
- Ruled out remotely: marginal/noisy triggering (rate flat 0.2-3.0 V, both edges); gates not
  opening per trigger or ignoring the delay (dead-time slope ~1 vs delay); wrong scoring spot
  for the old hunt (re-scored at (501, 454): flat); threshold-limited triggering.
- Found: (1) very little signal-arm light reaches (501, 454): beam core ~30x weaker than at
  12:21 and still fading; the ring and the (750, 440) blob are in the blocked-arm run too, so
  they are not the beam. (2) The dominant background is ungated, wall-clock-proportional
  (~200 hits/s over the 150x150 ROI) so raw hit scores follow the trigger-rate drift and the
  dead time; scans must use the drift-free (spot minus far-field) score. (3) The trigger rate
  fell from ~40 kHz to ~22-25 kHz during the session.
- Needs a person in the lab (cannot be tested remotely):
  1. Measure the idler-signal coincidence rate and delay on the time tagger, with the signal
     fibre's camera end on an APD. If coincidences/idler singles are << 1%, the 40 kHz triggers
     are mostly not SPDC idlers or the signal arm is badly coupled; fix that first.
  2. Check what drifts: pump power at the crystal, SPDC/idler alignment (trigger rate fell
     ~40 -> ~25 kHz in 1.5 h, beam at the camera fell ~30x since 12:21).
  3. Re-couple / re-align the signal fibre output onto the camera; at the current level
     (~10 detected photons/s at the wide gate) even a perfect herald gives < 1 heralded photon/s.
  4. Check whether the idler path to TRIGGER IN has extra delay (long cable, delay box, time
     tagger output): if the trigger arrives after the signal photon, no gate delay can catch it.
  5. Scope TRIGGER IN: the rate is flat up to 3 V, so the pulse is > 3 V there, not ~0.5 V.

## Data reorganized (20 Oct, ~16:55)
The 20 Oct runs and figures were moved from `data/raw/troubleshoot` and
`data/processed/troubleshoot_figures` into one folder per experiment under
`Desktop/HOI_freshstaret/data/troubleshoot/2026-10-20/`: `1221_notebook_session`,
`1346_photon_hunt`, `1459_agent_exp1_wide_gate`, `1515_agent_exp2_wide_long`,
`1522_agent_exp3_100ns_scan`, `1547_agent_exp4_500ns_scan`, `1614_agent_exp5_on_off_test`
(each with `runs/` and `figures/`). Figure paths in the entries above refer to the old layout.

## 2026-10-20 session 2 (remote; second agent; camera free, Tristan's kernel 13928 not holding it)

### E6 (17:00-17:39) Wide 5 us gate on/off test + short-gate reference (agent_scripts/multi_ab.py --tag wide_ab)
- Hypothesis: most background is ungated (E3), so a 5 us gate from 27 ns catches the heralded photon
  at any plausible delay (110-130 m fibre + cables) for little extra background.
- Trigger rate at the start (3 ns x 100k acc): 25.4 kHz.
- Settings: gain 9000, 50,000 acc; arms on = 27 ns + 5 us (0.25 s gate-open/readout), off = 15027 ns + 5 us
  (shift 15 us > width; gate 20 us after the trigger is still well inside the ~40 us mean trigger spacing),
  short = 15027 ns + 50 ns (beam/stray-light zero point); blocks 12/12/6 readouts, order rotated,
  36 cycles -> 432/432/216 readouts (21.6 M gates per arm). s/readout: on 1.53, off 2.31, short 2.04.
- Expected before: heralded 1e-5..3e-5 spot hits/trigger -> 0.5-1.5 hits/readout; per-readout noise ~4-5
  hits -> sigma(on-off) ~0.3 after 430 readouts each, i.e. ~2-5 sigma.
- Result (offline, agent_scripts/analyze_ab.py and a per-run per-pixel-baseline version; the off arm
  has 0.78 s more wall time per readout, fitted out with a hits-per-second term):
  r<=10 disk on-off = -0.14 +- 0.41 hits/readout; beam-profile matched filter -0.28 +- 0.17 (z -1.6);
  r<=6 -0.51 +- 0.23 (z -2.2, wrong sign for signal). Online drift-free net score: z +0.03. FLAT.
- Upper limit (2 sigma, 27-5027 ns): < 1.4e-5 heralded spot hits per trigger (matched filter;
  1.8e-5 with the r<=10 disk). Same per-trigger limit as E5 but for the whole plausible delay range.
- NEW finding, in-gate light: with the gate open the camera images a lit scene, not just the beam
  (figures/full_frame_in_gate.png): a square mount with a bright ring around (501, 454) and a dark
  centre, bright columns below, a blob at (750, 440), brighter bottom-left. ~0.04 hits/px per 0.25 s
  gate-open (= 0.16 hits/px/s, ~50 hits/s in r<=10), i.e. continuous stray (room/ambient or scatter)
  light that leaks into the imaging path. For a 5 us gate at 25 kHz this is ~3x the ungated background,
  so the 'wide gates are cheap' idea only holds for gates <~1 us. The ring + dark centre also explain
  the 'hole' at the spot seen in E3/E5 hit maps (the scene is dark exactly at the fibre image), so it is
  not intensifier gain loss as I first suspected (commit 65c5eb1 message).
- Beam today: no compact excess at (501, 454) on top of the scene; core (r<=6) minus annulus (10-16 px)
  -0.31 +- 0.10 hits per readout (0.25 s gate-open) -> signal-arm light < ~1 hit/s in the core, i.e.
  below even the E2 level. With no detectable signal-arm light, the expected heralded rate is
  ~0 and no gate setting can find it.
- Figures: data/troubleshoot/2026-10-20/1700_agent_wide_ab/figures/ (multi_ab.png, in_gate_light.png,
  full_frame_in_gate.png, records.json, analysis.json).

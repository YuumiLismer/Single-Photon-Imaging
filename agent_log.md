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

### E7 (17:57-18:39) Early window 27-527 ns on/off (multi_ab.py --tag early_ab; Tristan: "start at ~100 ns")
- Hypothesis: the heralded delay could be ~100 ns (shorter than the fibre estimate, e.g. if the trigger
  path is long too). A 500 ns gate has ~10x less in-gate stray light than 5 us (E6), so it is the most
  sensitive way to cover 27-527 ns.
- Trigger rate at the start: 21.4 kHz (down from 25.4 kHz at 17:00).
- Settings: gain 9000, 50,000 acc (25 ms gate-open/readout), on = 27 ns + 500 ns, off = 5027 ns + 500 ns
  (shift 10x the width; smaller dead-time mismatch than a 15 us shift), blocks of 12, order alternated,
  39 cycles -> 468 readouts per arm (23.4 M gates each). s/readout on 2.26, off 2.53.
- Result: online drift-free net r<=10: on -0.32 vs off -0.63, z +1.5. Offline (analyze_ab, per-run
  baseline + interval regression): r<=6 -0.04 +- 0.11, r<=10 -0.39 +- 0.17, matched filter -0.15 +- 0.08
  (z -1.9). FLAT; the two analyses disagree in sign at the 1.5-2 sigma level, which is the size of the
  dead-time systematic, so nothing here is a candidate.
- Upper limit 27-527 ns (2 sigma, matched filter): < 6.5e-6 heralded spot hits per trigger (r<=6: 7.3e-6),
  about 2x tighter than E5 (300-900 ns) and E6 (27-5027 ns), and below the 1e-5..3e-5 expected from the
  afternoon beam brightness. But see E6: the beam itself is not detectable today, so the expectation
  today is ~0.
- Figures: data/troubleshoot/2026-10-20/1757_agent_early_ab/figures/ (multi_ab.png, records.json, analysis.json).

### Session 2 summary (camera ~85 min: 17:00-17:39 and 17:57-18:39, closed with ht.shutdown each time)
- No heralded photons. Limits (2 sigma, heralded spot hits per trigger): 27-527 ns < 6.5e-6;
  27-5027 ns < 1.4e-5.
- New: with the gate open the camera sees a lit scene (stray light, ~0.16 hits/px/s) whose dark centre is
  exactly at (501, 454); the signal-arm beam is not detectable on top of it today (core minus annulus
  -0.31 +- 0.10 hits per 0.25 s gate-open). Wide gates are not cheap (5 us: stray light ~3x the ungated
  background); keep gates <= ~500 ns once the delay range is known.
- Stopped early: `git push` was refused by the session's permission system, so per the git rules I
  stopped and did not run agent_scripts/bg_probe.py (ungated-background origin + MCP bracket gating test;
  written, untested, not committed).
- Needs a person in the lab: (1) make the signal beam visible again at a wide gate (it was contrast ~40 at
  12:21, gone now): check pump power, SPDC/fibre coupling, the signal fibre end at the camera;
  (2) find and block the stray light (room lights/LEDs/pump scatter lighting the mount around the fibre
  output; cover the camera-to-fibre path); (3) the coincidence/delay check on the time tagger with an APD
  (session 1 list); (4) trigger rate fell again 25.4 -> 21.4 kHz in 1 h.

## 2026-10-21 session 3 (remote; Tristan: "try it with an internal trigger" / "just see if the camera works")
Data: `Desktop/HOI_freshstaret/data/troubleshoot/2026-10-21/` (0855_agent_health_timing, 0856_agent_health_switch,
0857_agent_health_gain_gate, 0911_agent_health_closed_bg). Scripts: agent_scripts/health_check.py (acquire),
health_analyze.py (figures + health.json). Camera free at start; sensor -20.0 C Locked throughout; no saturated pixels.

### H1 (08:55) Internal trigger takes
- Read-back Trigger Source = Internal / Frequency as set on every run; switching External <-> Internal twice works
  (read-back follows, timing follows). Readout interval = accumulations / rate to 0.1 %: 0.501, 1.001, 2.001 s at 1 kHz
  (500/1000/2000 acc), 1.001 s at 10 kHz x 10k and 100 kHz x 100k; floor 0.226 s = readout time (frame transfer).
- External (APD) in between: 0.546 s per 20k acc -> >= 36.6 kHz triggers (up from 21.4 kHz at 18:00 yesterday).
- Pedestal 640 ADU, read noise ~10 ADU (gain 1). Figure: 0856_agent_health_switch/figures/timing.png.

### H2 (08:57-09:10) Intensifier response, gate closing, single photons (internal, gate delay 27 ns)
- Gain (10 kHz x 10k acc, 10 us gate = 0.1 s open, minus 3 ns gate): in-gate light 6.6, 16.2, 28.6, 41.4 ADU/px/readout at
  gain 1000/3000/6000/9000 (~proportional; 0.3 +- 1.0 at 100). Read-back sub-gains: ADC Quality "Electron Multiplied",
  emICCD gain mode Optimal (EM register is in use).
- Gate (gain 9000, 1 s wall): hit px/px 0.0063 (3 ns), 0.0068 (100 ns), 0.0088 (1 us), 0.030 (10 us), 0.116 (50 us):
  linear in gate-open time, ~0.22-0.24 hits/px per s open. Same 0.1 s open at 1 / 10 / 100 kHz: 0.029 / 0.030 / 0.030,
  and 3 ns gates at 1 vs 100 kHz equal: the gate opens and closes properly; the number of gates does not matter.
- Events: open-gate events are multi-pixel blobs (median 3 px, peak 253 ADU, summed 532 ADU, long charge tail) = intensifier
  photon events; closed-gate events are mostly SINGLE pixels (median 1 px, summed 280 ADU) -> not photocathode events.
- Beam: no compact spot at (501, 454). Core(r<6) minus annulus(10-16) = -6,700 +- 1,000 ADU per s gate-open (0.5 s runs)
  vs +123,000 at 12:21 (high pump) and -2,700 with the arm blocked. The scene is the session-2 stray-light scene (mount,
  bright ring, dark centre at the spot, stripes, blob ~(756, 446)); no compact spot elsewhere beyond scene structure.
- Figures: 0857_agent_health_gain_gate/figures/{gain,gate,full_frame,single_photons,events}.png, health.json.

### H3 (09:11-, in progress at handback) Closed-gate background vs wall time
- 3 ns gate, gain 9000, 1 kHz: hits/readout (ROI) 13, 40, 126, 476, 1005 at 0.25, 0.5, 1, 2, 3 s (4 s: ~1650 in H2);
  ROI mean rises 635 -> 697 ADU. SUPERLINEAR (~t^1.8), not proportional to time, independent of gate count.
  Intensifier-off and gain-1 controls were still running at handback (script closes the camera with ht.shutdown).
- Hypothesis: EM-register-amplified CCD signal (single-pixel events, gain dependent) growing with exposure, e.g. a
  light leak onto the CCD/phosphor or CCD dark/heating; to be decided by the intensifier-off run.

## 2026-10-21 session 4 (remote; Tristan: "image these single photons with an EXTERNAL trigger", gates <= 10 us, readouts <= 1 s)
Data: `Desktop/HOI_freshstaret/data/troubleshoot/2026-10-21/` (0949_agent_rate_check, 0950_agent_dscan, 1047_agent_dscan2,
1053_agent_dscan3). Camera free at start (Tristan's kernel 13928 not holding it); -20 C Locked; no saturation; closed with
ht.shutdown after every script. Camera time ~75 min (09:49-11:05).

### S1 (09:49) External trigger rate (agent_scripts/rate_check.py; 0.5 V, rising, 50 ohm, DC; 3 ns gate, 20k/30k acc)
- 54.7 / 55.0 kHz (0.37 s per 20k, 0.55 s per 30k triggers). Up from 36.6 kHz at 08:55 and 21-25 kHz yesterday.

### S2 (09:50-11:05) Delay search 527-3027 ns, five 500 ns windows, interleaved ON/OFF (multi_ab.py)
- Hypothesis: with 110-130 m fibre (~550-640 ns) the heralded photon sits somewhere in 527-3027 ns (27-527 ns done in E7).
- Settings: gain 9000, External, ON = d .. d+500 ns, OFF = d+5000 ns (same width), d = 527/1027/1527/2027/2527; blocks of
  20 readouts, 10 arms rotated per cycle. Pass 1 40k acc (20 ms open/readout), 12 cycles; pass 2 30k acc, 1 cycle; pass 3
  20k acc, 3 cycles. ~11 min per window in total.
- Trigger rate was NOT stable: 52 kHz (cycles 1-3) -> 18 kHz (cycles 6-7, ~10:12-10:22) -> 42 kHz -> ~25-30 kHz at 11:00.
  So in pass 1 readouts went to 2.1-2.4 s (cycles 5-9), breaking the 1 s rule. I could not stop the run (process kill was
  refused by the permission system), so I added a readout-interval guard to multi_ab.py (dc1a2db), which stopped passes 2
  and 3 cleanly at 1.13 / 1.05 s. OFF arms always take ~25% longer per readout (later gate -> more missed triggers), fitted
  out with the readout-interval regression.
- Result (window_figures.py: per-run baseline, beam-profile matched filter, regression on arm + readout interval;
  combine_windows.py: per-trigger inverse-variance over the 3 passes). Net heralded spot hits per trigger, all passes:
  527 ns +5.1e-6 +- 2.8e-6 (z +1.8); 1027 ns +5.7e-6 +- 2.7e-6 (z +2.1); 1527 ns +1.7e-6 +- 2.7e-6 (z +0.6);
  2027 ns +5.3e-6 +- 2.7e-6 (z +1.9); 2527 ns -4.3e-6 +- 2.8e-6 (z -1.5). Pass 1 alone: z +1.7/+1.8/+0.5/+1.9/-2.2;
  pass 1 with only readouts <= 1.05 s: z +1.2/+2.0/-1.1/+1.0/-0.9. No window >= 4 sigma, so no halving and no target image.
- Conclusion: FLAT from 527 to 3027 ns. 2-sigma upper limits (spot hits per trigger): 527 ns 1.1e-5, 1027 ns 1.1e-5,
  1527 ns 7.0e-6, 2027 ns 1.1e-5, 2527 ns 5.5e-6. The similar +5e-6 in three windows is common to non-adjacent
  windows, so it is a shared systematic (dead-time/background correction), not a localized heralded photon.
- Figures: 0950_agent_dscan/figures/window_<d>ns.png (pass 1, ON / OFF / ON-OFF), window_<d>ns_le1s.png (readouts <= 1.05 s),
  summary_vs_delay.png (pass 1), summary_vs_delay_combined.png + windows_combined.json (all passes);
  1047_agent_dscan2/figures and 1053_agent_dscan3/figures: the same per pass.
- Lab checks: (1) the herald/idler rate drifts 18-55 kHz within an hour (is the idler counter or threshold behaving?);
  (2) the signal beam is still not visible at a wide gate (session 3 H2), so ~0 heralded photons are expected at the camera
  whatever the delay; (3) time-tagger coincidence check of the signal-idler delay.

## 2026-10-21 session 5 (remote; Tristan in the lab: signal arm re-coupled, "increased the counts in the signal arm")
Data: `Desktop/HOI_freshstaret/data/troubleshoot/2026-10-21/` (1400_agent_beam5 [stopped], 1400_agent_beam5_2, ...).
Camera free at start; -20 C Locked; no saturation; ht.shutdown after every script.

### B1 (14:00-14:04) Beam check, External trigger (agent_scripts/beam_check.py, new)
- Hypothesis: after re-coupling, the compact signal spot is back (12:21 on 20 Oct: +123,000 ADU per s gate-open).
- Trigger rate: 68.0 kHz (14:00), 66.4 kHz (14:01), 88.4 kHz (14:04) (20k acc x 3 ns gate).
- First try stopped itself: 54k acc x 10 us took 1.36 s/readout, because a 10 us gate blocks triggers while
  open (40 k gates/s at 68 kHz = R/(1+R*10 us)). Accumulations now come from the dead-time-corrected rate.
- Settings: gain 9000, 31,000 acc; W = 27 ns + 10 us (0.31 s open), Woff = 10027 ns + 10 us, C = 3 ns gate;
  3 cycles x 10 readouts per arm (0.66 / 0.97 / 0.36 s per readout).
- Result: BEAM IS BACK, compact, still at (501, 454). Core(r<6) minus annulus(10-16): W-C +123,000 +- 15,000 ADU/s
  gate-open (~230 detected photons/s at 532 ADU/photon; DoG z 99), Woff-C +117,000 +- 26,000 (same light: mostly
  unheralded, as expected), W-Woff +6,900 +- 4,000 at the spot (1.7 sigma; if real, ~1.3e-4 heralded hits per gate).
  Same in-gate stray scene as sessions 2-3 (mount square, blob (750, 450), column below), now with the spot on top.
- Expectation: heralded hits/trigger = 230/s x (idler-side heralding efficiency) / trigger rate; for 2-20 % at 66 kHz
  that is 7e-5..7e-4 per trigger, 5-50x above the session 2-4 limits -> a 500 ns window should give 5 sigma in
  ~20-60 readouts of 30k acc.
- Figure: 1400_agent_beam5_2/figures/beam.png, beam.json.

### D1 (14:03-14:09) Delay search 27-3027 ns, six 500 ns windows, ON/OFF interleaved (multi_ab.py --tag dscan5)
- Settings: gain 9000, External, 30,000 acc (15 ms open/readout), ON = d..d+500 ns, OFF = d+5000 ns, d = 27..2527;
  blocks of 20, 12 arms rotated, 3 cycles -> 60 readouts per arm; 0.24-0.34 s (ON) / 0.41-0.48 s (OFF) per readout
  (trigger rate ~90-120 k gates/s; no readout near the 1 s limit).
- Analysis: window_figures.py, new default estimator (matched filter minus 4 control apertures, ON-OFF, Welch errors;
  the old readout-interval regression is collinear with short readouts and gave -2..-4 sigma everywhere).
- Result (net spot hits/readout, z): 527 ns +1.58 +- 0.43 (z +3.7); 1027 -0.48 +- 0.42; 1527 -0.48 +- 0.49;
  2027 -0.43 +- 0.48; 2527 +0.85 +- 0.44 (+1.9). Net image at 527 ns: compact red spot exactly on (501, 454).
- NEW, important: the 27-527 ns ON gate collects NO light at all, not even the stray scene (far field 0.59e-3 hits/px
  = the closed-gate floor, vs 2.5e-3 in every other gate; core 0.07 vs 1.5-2 hits). Read-back is correct
  (27 ns, 500 ns, External). With the External trigger, gates that end before ~527 ns after the trigger do not open.
- Figures: 1403_agent_dscan5/figures/window_<d>ns.png, summary_vs_delay.png, windows.json.

### D2 (14:10-14:17) 250 ns ladder (multi_ab.py --tag ladder250)
- Settings: 50,000 acc, ON 277/527/777 ns + 250 ns, OFF +5000 ns, 5 cycles -> 100 readouts per arm (0.39-0.44 s ON).
- Result: 527-777 ns +2.34 +- 0.43 spot hits/readout, z +5.5 (4.7e-5 per trigger); 777-1027 -0.20 +- 0.34 (z -0.6);
  277-527 again dark (far field 0.94e-3 vs 3.0-3.6e-3: the gate does not open).
- Figures: 1410_agent_ladder250/figures/window_527ns.png, summary_vs_delay.png.

### D3 (14:17-14:26) 25 ns gates 527-802 ns (multi_ab.py --tag fine25; 50k acc, 40 readouts per arm)
- No window above +0.32 +- 0.21 spot hits/readout (all |z| < 1.6); sum of the 10 windows +0.75 +- 0.63 vs +2.34 +- 0.43
  expected from D2. Taken during the herald-rate spike (150-190 kHz, see H below), when most triggers were probably
  not idler photons, so the heralded fraction per trigger was diluted. Figures: 1417_agent_fine25/figures/.

### P1 (14:27) Gate-opening probe, 1 us gates at 27..1027 ns and 5027 ns (--tag gateprobe, 10 readouts each)
- All 1 us gates open (far field 3.9-4.4e-3 hits/px). So the dark gates of D1/D2 (27+500, 277+250) share only one
  feature: they END at 527 ns. Working hypothesis: with this herald input, gates that close before ~530 ns after the
  trigger are not opened (cause unknown; to be checked with 100 ns probes in D5).

### D4 (14:29-14:38) Confirmation (--tag confirm, 40k acc; stopped by the 1.05 s guard in cycle 4 as the rate fell)
- 60 readouts per arm: 527-777 ns +2.06 +- 0.34 spot hits/readout, z +6.1 (independent repeat of D2);
  527-1027 +1.80 +- 0.47 (z 3.9); 27-1027 +2.71 +- 0.67 (z 4.0); 777-1027 -0.15 +- 0.31.
- Figures: 1429_agent_confirm/figures/window_528ns.png (the 527-777 ns window), summary_vs_delay.png.

### D5a (14:39-14:50) 50 ns gates 527-777 ns (--tag fine50, 20k acc, 160 readouts per arm, 0.25 s/readout)
- 527-577 ns +0.456 +- 0.097 spot hits/readout, z +4.7; 577-627 +0.08 +- 0.06; 627-677 -0.04; 677-727 -0.03;
  727-777 +0.04. The heralded photon arrives in 527-577 ns (fibre 110-130 m predicts ~550-640 ns incl. cables).
- Figures: 1439_agent_fine50/figures/window_527ns.png, summary_vs_delay.png.

### H (offline) Herald rate vs time (agent_scripts/herald_rate_history.py, all 735 External runs since 20 Oct)
- 21 Oct: 55-60 kHz (09:50) -> 18 kHz (10:15) -> 45 kHz (10:40) -> 20 kHz (11:05); 70 kHz (14:00) -> 175 kHz
  (14:20) -> 50 kHz (14:40) -> 75 kHz (14:50). 20 Oct: 44 kHz (15:30) decaying to 22 kHz by 16:30.
- Figure: 1439_agent_herald_rate/figures/herald_rate_vs_time.png (+ .json).

### D5b (14:50-14:57) 10 ns gates 497-597 ns + 100 ns probes (--tag fine10, 40k acc, 40 readouts per arm)
- All 10 ns windows flat (|z| < 2.5, no positive excess). 100 ns probe gates ending at 427-552 ns: far field 0.36-0.46e-3
  = floor -> NOT open. Explanation found in P2/D6: 10 ns gates ending before ~577 ns do not open.

### P2 (14:58) Gate-end probe (--tag gateend): delay 27 ns, widths 500..1000 ns, 10 readouts, 30k acc
- Far field (hits/px x1e3): end 527 -> 0.55 (closed), end 552 -> 1.41 (partly), end 577 -> 2.91 (open), 602-1027 ->
  2.9-4.0. A gate opens fully from its programmed start only if it ENDS later than ~560-577 ns after the trigger.
  Cause unknown (camera/trigger-input behaviour with this herald signal); it hides any gate that closes earlier.

### D6 (15:00-15:11) Gate-start scans with fixed gate end (--tag endfix [end 592], endfix582 [end 582]; 40k acc, 60 readouts)
- End 592: start 447 +2.27 +- 0.47 (z 4.9), 487 +1.29 +- 0.33, 507 +2.46 +- 0.40 (z 6.2), 527 +1.64 +- 0.37,
  547 +2.45 +- 0.31 (z 7.8), 567 -0.03 +- 0.15.
- End 582: start 537 +2.16 +- 0.35 (z 6.3), 547 +1.97 +- 0.29 (z 6.8), 552 -0.05 +- 0.14, 557 +0.02, 562 -0.06,
  567 0.00, 572 +0.24 +- 0.12.
- => All heralded photons arrive between 547 and 552 ns (programmed delay) after the herald; nothing earlier (starting
  at 447 adds nothing) and nothing later. Edge sharper than the 5 ns step: FWHM < ~5 ns.
- Arrival-profile figure: 1512_agent_target/figures/arrival_profile.png (+ .json; agent_scripts/arrival_profile.py).

### T1 (15:12-15:17) TARGET heralded single-photon image (agent_scripts/target_image.py)
- Settings: External herald (0.5 V rising, 50 ohm), gain 9000, gate 542 ns + 40 ns, 80,000 acc/readout
  (3.2 ms open), ON / off-delay (+5 us) blocks of 50 alternated, 4 cycles -> 200 + 200 readouts;
  ON 0.37 s/readout (~215 k gates/s), OFF 0.77 s/readout. No saturation, -20 C Locked.
- Result: HERALDED SINGLE PHOTONS DETECTED: +3.17 +- 0.26 net spot hits per readout, z = 12.1 over the off-delay
  background (ON 4.18, OFF 1.01 spot hits/readout) = 4.0e-5 heralded spot hits per trigger, ~8.5 heralded photons/s
  of wall time. Compact spot at (501, 454) in ON - OFF; nothing elsewhere.
- threshold_sweep: significance flat ~18.7 sigma for 60-110 ADU above pedestal (best 76 ADU), falling above ~130 ADU;
  200 ADU (used above) is conservative. plot_wiener SNR is negative because plain ADU subtraction is dominated by the
  ungated background, which is larger in the OFF arm (0.77 vs 0.37 s per readout); use hit counting, not ADU sums.
- Figures: 1512_agent_target/figures/heralded_single_photon_image.png (main), target_threshold_sweep.png,
  target_plot_run.png, target_raw_vs_sum.png, target_wiener.png, target.json, arrival_profile.png.

### Session 5 summary (camera 14:00-15:17, ~77 min; ht.shutdown after every script)
- Beam back (+123,000 ADU/s, ~230 detected photons/s at the spot). Heralded photon delay 547-552 ns (programmed gate
  delay), FWHM < ~5 ns. Target image z = 12.1 at 542 + 40 ns. Herald rate swung 50-215 k gates/s during the session.
- Camera quirk: with the External herald, gates that end earlier than ~560-577 ns after the trigger never open.
  The 10 ns / 25 ns fine scans were blind for that reason. Any short gate must end >= ~580 ns.

## 2026-10-21 session 6 (remote; goal: narrowest gate that still images the heralded photons)
Data: `Desktop/HOI_freshstaret/data/troubleshoot/2026-10-21/` 1717_agent_attrs ... 1757_agent_narrow6. Camera 17:17-18:00
(~45 min), -20 C Locked, no saturation, ht.shutdown after every script; all changed settings read back and restored.

### R1 (17:19) Re-check (1719_agent_recheck6): herald 208 kHz (3 ns gate); 542+40 ns ON/OFF, 60k acc, 40+40 readouts: z 8.4. Spot still there.
### R2 (17:31) Trigger rate vs threshold (agent_scripts/threshold_rate.py, 1731_agent_thr_rate): flat 233 kHz from 0.2 to 1.3 V.
  The herald pulses are now well above 1.3 V (no longer marginal at 0.5 V); the rate is real pulses, not threshold noise.

### Q: the External gate-end quirk (runs gatefire, gateend6, intrate6, noclean6, fullframe6, gain6, cutofftime6)
- Internal trigger (50, 100, 200, 400 kHz): 27-527 ns and every other gate open normally. External: 27-527 never opens.
  => specific to the External trigger path, not to the gate pulser rate.
- End scan, start 27 ns (spot hits r<=10 per readout, open ~20): end 517-547 -> 0; 557 -> 2; 567 -> 12; 577 -> 19; 587 -> 22.
  Start 400 ns gives the same cutoff (ends 527/547 dark, 567 partial): it is the gate END that matters.
- No setting changes it: falling edge, threshold 0.3 / 0.8 V, AC coupling, high-impedance termination, Clean Until
  Trigger off, Full Frame readout, gain 3000/1000 all keep 27-527 dark and 27-577+ open. Gating Speed (Fast) and
  Intensifier Options (None) are read-only; SyncMASTER off. Full attribute dump: 1717_agent_attrs/figures/attrs.json.
- Per-readout it is all-or-nothing and switches in time: near the cutoff a run is either open or closed for whole
  readouts, e.g. 27-557 ns at 20k acc: one block open for 7 readouts (~1.6 s) then closed, other blocks closed throughout.
  Short gates are not partly open; the whole gate fires or not.
- Best explanation: a camera-internal latency of ~550-575 ns in handling an External trigger (jittering / switching
  between states over seconds, as if the trigger is re-synchronised to an internal clock); a gate that has to be
  finished before that never fires. Not fixable from software. It is not caused by the herald signal level.

### A: analysis fix (offline). window_figures/target_image used a per-pixel median baseline over each run's readouts;
  with the stronger, compact spot it removed about half the heralded hits. Default is now the per-readout pedestal.
  Session-5 target (542+40 ns) re-analysed: z 12.1 -> 17.0, 4.0e-5 -> 7.4e-5 heralded hits per trigger.

### S: gate-start scan, end fixed at 587 ns (1738_agent_startscan6; 100k acc, 50 readouts/arm, shared OFF 1546+41 ns)
- Net spot hits/readout: start 540-548 +8.8..+16.2 (z 10-14 each), 550 +2.9 +- 0.6, 552 +0.2, 554 -0.5, 556 -0.3.
  => photons arrive at 548-550 ns (~80 % of them between 548 and 550, the rest 550-552): FWHM < ~2 ns
  (step-limited). Heralded rate now ~1.2-1.8e-4 spot hits per trigger (2-4x session 5).
- Same run: 546+31 (end 577) +13.8 +- 1.2 (z 11); 546+26 (end 572) +8.9 (~55 %, one block fully closed).
- Arrival profile: 1751_agent_target31/figures/arrival_profile.png.

### P: readout post-selection by trigger rate (window_figures/target_image --rate-band). Within runs the herald rate
  is steady (target ON readouts 0.368-0.370 s); on the session-5 confirm run a 10 % band drops 2/3 of readouts and
  z falls 7.1 -> 4.1 (signal per readout not higher in the kept ones). Does not help; off by default.

### T2 (17:51-17:56) TARGET at 546 + 31 ns (agent_scripts/target_image.py, 1751_agent_target31)
- External herald ~199 k gates/s, gain 9000, 120,000 acc/readout (3.7 ms open), OFF = +1 us (1546-1577 ns),
  blocks of 50, 200 ON + 200 OFF readouts (0.60 / 0.73 s per readout), 4.5 min of readouts.
- Net +21.4 +- 0.65 spot hits/readout, z = 33.0 (ON 22.8, OFF 1.4) = 1.8e-4 per trigger, ~36 heralded photons/s.
  vs 40 ns session 5 (z 12.1 old baseline / 17.0 new in 4.6 min): z per sqrt(min) 15.7 vs 8.7.
- Figures: heralded_single_photon_image.png, target_plot_run.png, target_threshold_sweep.png (best 58 ADU),
  target_wiener.png, target_raw_vs_sum.png, arrival_profile.png.

### N (17:57) Narrower gates (1757_agent_narrow6; 120k acc, 40 readouts each, shared OFF 1546+31)
- 546+31 +25.6 +- 1.5 (z 16.9); 548+29 (end 577) +20.0 +- 1.2 (z 16.4); 548+26 (end 574) +10.3 +- 1.6 (z 6.3);
  549+24 (end 573) +8.9 +- 1.3 (z 6.8). Below an end of ~577 ns the gate drops out for whole readouts.
- Narrowest gate still imaging the photons: 24 ns (z 6.8 in 40 readouts, ~60 % of photons lost);
  narrowest efficient gate: 29-31 ns ending at 577 ns.
- A 3-10 ns gate needs the photons to arrive >= ~580 ns after the trigger: delay the signal by >= ~35 ns relative to
  the herald (shorter herald coax by ~7 m, or ~7 m more signal fibre).

## 2026-10-23 herald-rate-keeper (remote; pump laser back on; goal: S.Pos2 controls herald rate, hold 15-20 kHz)
Data: `Desktop/HOI_freshstaret/data/troubleshoot/2026-10-23/` 1149_herald_rate_lock ... 1156_herald_rate_lock. Camera 11:49-12:07
(rate readings only: gain 1, 3 ns gate, threshold 0.5 V rising; ht.shutdown after every script). Camera free at start.
XPS: S.Pos2 28.176 deg, group S status 13 (ready) at start; S.Pos1 0 deg, untouched; group P untouched.
Camera/XPS link via switch at 100 Mbps: readout floor still 0.226 s; only the first 20k-acc readout of each point was
readout-limited ('>='; excluded automatically). Good intervals 0.4-2.5 s, no change of accumulations/target needed.

### 11:49 Rate at the start angle, before any move (1149_herald_rate_lock, --measure): 624.9 +- 1.0 kHz at 28.176 deg.
  2.5x the highest rate seen before (~233 kHz). Laser clearly on. Default +45 deg calibration would cross the maximum,
  so a --max-rate ceiling (default 250 kHz) and negative --span / --start-angle were added before sweeping.
### 11:51 +2 deg test move (1151_waveplate_rate_check, ceiling 660 kHz): 28.176 -> 624.7, 30.153 -> 681.4 kHz (ceiling
  hit, returned). S.Pos2 moves fine; rate goes down for decreasing angle.
### 11:52 Response check (1151_waveplate_rate_check_2): offsets 0/-5/-10/-15/-22.5/0 -> 623.3 / 467.0 / 297.5 / 150.0 /
  59.0 / 631.0 kHz. Clear cos^2: **S.Pos2 is the pump waveplate**. Returned to 28.176.
### 11:54 Calibration (1154_herald_rate_lock, --start-angle 15 --span -30 --step 2.5, ceiling 250 kHz, --park):
  15.0 172 | 12.6 115 | 10.1 81.5 | 7.6 72.1 | 5.1 84.6 | 2.6 119 | 0.1 173 | -2.4 242 | -4.9 319 kHz (stopped).
  Fit B + A cos^2(2(th - th0)): min B = 74.6 kHz at 7.55 deg (mod 90); max at 52.55 deg (mod 90); fitted max 1.48 MHz is
  an extrapolation (only the bottom of the curve was swept), measured max so far 681 kHz at 30.15 deg.
  **Extinction floor ~70-75 kHz, ~10 % of 680+ kHz: 15-20 kHz is out of reach with S.Pos2.** The floor also moved:
  59 kHz at 5.7 deg at 11:53 vs 72-85 kHz at 5-7.6 deg at 11:55 (drifting pump polarization after warm-up and/or
  a pump-independent background). Figure: 1154_herald_rate_lock/figures/waveplate_calibration.png.
### 11:56-12:07 Hold 10 min (1156_herald_rate_lock, --hold, band 15-20 kHz): band out of reach, so the lock parked at the
  fitted minimum 7.546 deg (no moves needed). 87 points (5 readings each): 0 % in band; 71.4 kHz mean, 67.9-76.2 kHz,
  sd 2.2 kHz, slow +-4 kHz wander over minutes.
### Rate-vs-time plot (every reading of the session, band shaded, angle panel):
  1156_herald_rate_lock/figures/herald_rate_vs_time.png (also herald_rate_vs_time_lock.png via
  `agent_scripts/herald_rate_history.py --lock-logs --since 2026-10-23T11:00`).
### Left: S.Pos2 at 7.546 deg (~71 kHz, the minimum), group S ready (status 12, motors on). Total motion < 2 min.
### Next: the waveplate+PBS cannot go below ~70 kHz. To reach 15-20 kHz: check the pump polarization before the HWP
  (elliptical light or a poor PBS limits extinction), add an ND filter / lower laser current, or check for a
  pump-independent herald background (block the pump: rate should drop to dark counts).
### Code: trigger_rate_monitor.monitor(on_reading=, stop_after_good=); hold_herald_rate.py logs every reading/move to
  herald_rate_readings.csv and plots rate+angle vs time, --measure, --max-rate, negative --span, --start-angle, --park,
  50 deg move cap, status check, accumulations carried between points, no sub-0.01 deg lock moves;
  waveplate_rate_check.py uses the same logging, --max-rate; herald_rate_history.py --lock-logs.

## 2026-10-23 session 7 (remote; goal: re-image heralded photons at the lowest herald rate, then narrow the gate)
Data: `Desktop/HOI_freshstaret/data/troubleshoot/2026-10-23/` 1241_herald_rate_lock ... 1302_herald_rate_lock. Camera 12:41-13:03
(~22 min), free at start, -20 C Locked, no saturation, camera closed after every script. Only S.Pos2 moved (status 12).

### W1 (12:41-12:43) Re-find the pump minimum (hold_herald_rate.py --measure / --calibrate --start-angle 2.5 --span 10 --step 1.25 --park)
- 7.546 deg: 59.4 kHz. Sweep 2.5..12.5 deg: 81.1 / 67.8 / 60.2 / 59.3 / 64.2 / 73.8 / 89.4 / 110.3 / 136.2 kHz.
  Fit: min 59.4 kHz at 5.95 deg (moved 1.6 deg since 11:54). Parked at 5.95 deg. 1242_herald_rate_lock/figures/waveplate_calibration.png.

### T3 (12:44-12:49) Target at 546+31 ns at the minimum rate (target_image.py --acc 45000 --cycles 4 --off-shift 1000, 1243_agent_target31_lowrate)
- Hypothesis: the heralded spot is still there; per-trigger contrast improves at 1/3 of the session-6 rate.
- External ~61 k gates/s, gain 9000, 45,000 acc (1.4 ms open), 200 ON + 200 OFF readouts, 0.73 s/readout.
- Result: NO heralded signal: net +0.09 +- 0.09 spot hits/readout, z 1.0 (ON -0.04, OFF -0.12; session 6 expectation
  at 1.8e-4/trigger: ~8 per readout). ON and OFF hit images identical; no spot at (501, 454).
- Figures: 1243_agent_target31_lowrate/figures/heralded_single_photon_image.png, target_plot_run.png.

### B2 (12:50-12:52) Wide-gate beam check at 59 kHz (beam_check.py, 1250_agent_beam_lowrate)
- 29k acc x 10 us (0.29 s open), 30 readouts per arm W (27 ns + 10 us) / Woff (10027 ns + 10 us) / C (3 ns).
- Result: SIGNAL BEAM GONE. Old spot (501, 454) core-annulus W-C -1,474 +- 1,120 ADU/s gate-open (session 5: +123,000,
  ~230 photons/s at a similar herald rate of ~68 kHz). Only the usual stray scene (mount square, post at (350, 462), column),
  same in W and Woff. Woff readouts 1.08 s (slightly above 1 s, see code fix below). Figure: figures/beam.png.

### B3 (12:53-12:55) Beam check at higher pump (S.Pos2 14.5 deg, 175-188 kHz; 1253_agent_beam_188k)
- Hypothesis: the ~60 kHz floor is pump-independent counts, so at the minimum there are hardly any pairs; more pump
  should bring the beam back.
- Result: no beam either: old spot W-C -425 +- 203 ADU/s (< ~1 photon/s). Full-frame pump-dependence image
  (175 kHz minus 59 kHz W-C, figures/pump_dependence_W-C.png): no positive feature anywhere (max z 4.7, a single
  point); the stray scene is ~7 % DIMMER at the higher pump, so it is not SPDC light. No pump-scaling light reaches the sensor.
- (beam_check's own rate read 88 kHz here: readout-limited, 20k acc / 0.226 s. Fixed, see code.)

### T4 (12:56-12:59) Session-6 target repeated at the session-6 rate (546+31 ns, 120k acc, 100+100 readouts, 1256_agent_target31_175k)
- External ~167 k gates/s (session 6: 199 k), 0.72 s/readout. Net -0.11 +- 0.15 spot hits/readout, z -0.7
  (session 6: +21.4 +- 0.65, z 33). Upper limit ~2.5e-6 heralded hits per trigger (2 sigma), ~70x below session 6.
- => The heralded photons are gone together with the signal beam. Not a gate/delay or herald-rate effect.

### Code: beam_check.py re-measures the rate with more accumulations when readout-limited (88 -> 173 kHz) and sizes
  accumulations by the slower Woff arm (Woff 1.12 -> 0.77 s/readout at 168 kHz). Tested on the camera (1259/1300_agent_beamcheck_ratefix_test*).

### Herald rate this session: 59-60 kHz at the minimum (12:42), 175-188 kHz at 14.5 deg (12:53-12:59), 70.6 kHz back at
  5.95 deg (13:03: the floor drifts +-10 kHz). Plots: 1302_herald_rate_lock/figures/herald_rate_vs_time_lock_session7.png
  (every reading + angle), herald_rate_vs_time_runs_session7.png (run meta).
### Left: S.Pos2 at 5.95 deg (minimum, ~70 kHz), group S status 12. Camera closed.
### Not done: steps 3-4 of the plan (narrower gates, contrast vs herald rate) need the signal beam; with no signal
  photons at the camera there is nothing to compare.
### Lab checks (physical, in this order): (1) is the signal fibre still connected and coupled? Measure signal-arm singles
  on an APD / power meter, and signal-idler coincidences; since 21 Oct the rate at S.Pos2 28 deg rose from ~200 to 625 kHz,
  so something on the pump/idler side changed after the laser restart, and the crystal/phase-matching or the signal-fibre
  coupling may need re-alignment. (2) Any shutter, cap, filter or beam block in the signal path to the camera (object, lenses)?
  (3) Did the imaging optics move (spot would land elsewhere; but no pump-dependent light was found anywhere on the full frame)?
  (4) Is the ~60-70 kHz herald floor pump-independent (block the pump: rate should drop to dark counts)?

## 2026-10-26 session 8 (remote; Tristan: "run everything at the lowest possible kHz ... find the herald photons"; goal: image them)
Data: `Desktop/HOI_freshstaret/data/troubleshoot/2026-10-26/` 0921_herald_rate_lock ... 0947_herald_rate_lock. Camera 09:21-09:48
(~27 min), free at start, -20 C Locked, no saturation, camera closed after every script. Only S.Pos2 moved (status 12 throughout).

### R1 (09:21) Rate at the start angle, before moving: 303.5 +- 0.7 kHz at 16.787 deg (23 Oct: ~100 kHz at this angle).
### R2 (09:22-09:24) Fine sweep 13 -> -3 deg, 1 deg steps (0922_herald_rate_lock): 168.7 / 137.6 / 108.7 / 82.8 / 60.6 / 42.2 / 28.0 /
  19.3 / **15.3 (5.06 deg)** / 15.6 (4.04) / 20.8 / 30.9 / 46.5 / 65.4 / 87.2 / 112 / 141 kHz. Fit: min 16.0 kHz at 4.65 deg, max ~1.87 MHz.
  NEW: extinction is now ~0.8 % (floor 15-30 kHz) vs ~5-10 % on 23 Oct (floor 59-75 kHz): 15-20 kHz IS reachable now. The pump
  side changed since 23 Oct (rate at a given angle ~3x higher, much better extinction). Figure: 0922_herald_rate_lock/figures/waveplate_calibration.png.
### R3 Floor vs time at the minimum (angle stable at 4.55-4.65 deg): 15.3 kHz (09:23), 24.8 / 27.4 (09:25), 22.6 (sweep 09:26-09:27,
  0926_herald_rate_lock), 29-30 (09:38, threshold scan), 21-25 k gates/s (target, 09:39-09:45), 17.3 (sweep 09:46), 18.8 kHz (09:47).
  The floor wanders 15-30 kHz on minute scales; the minimum angle does not move. Re-parking at the minimum was enough.
  Plots: 0947_herald_rate_lock/figures/herald_rate_vs_time_lock_session8.png/herald_rate_vs_time_lock.png (waveplate/lock readings + angle),
  herald_rate_vs_time_runs_session8.png/herald_rate_vs_time.png (every External run of the session incl. 76/178 kHz beam checks).

### B1 (09:27) Wide-gate beam check at the minimum, ~25 kHz (beam_check.py, 0927_agent_agent_beam_min)
- Hypothesis: after 3 days the signal beam may be back (Tristan may have worked on the setup).
- Settings: gain 9000, External 0.5 V rising, 13,000 acc; W = 27 ns + 10 us, Woff = 10027 ns + 10 us, C = 3 ns; 3 cycles x 10 readouts (0.69-0.90 s/readout).
- First result looked like a NEW compact spot at (907, 781), W-C DoG z 99, W-Woff z 82 (+61,600 +- 61,400 ADU/s). The error bar = the
  value: it was a single transient (232k ADU core, max pixel 54,000 ADU) in ONE of the three W runs (W_c3); the other 8 runs show
  0.2-8k. Not light. Fixed beam_check: per-pixel median of the per-run frames (commit 80ae50b), re-analysed: beam_robust.png.
- Robust result: no spot. Old spot (501, 454): W-C -427 +- 2,190, W-Woff -129 +- 1,349 ADU/s gate-open (session 5: +123,000 at 66 kHz).
  W-C shows only the usual stray scene (mount square, post at (355, 460-507), static point (232, 159)), identical in Woff.
- Heralded reading of the same data (W-Woff covers heralded arrivals 27 ns - 10 us, full frame, gates that certainly open):
  old spot -0.03 +- 0.33 photons/readout vs ~2.3 expected at session-6 efficiency (1.8e-4/trigger x 13k gates): ~7 sigma below.
- Figures: 0927_agent_agent_beam_min/figures/beam.png (plain mean, with the transient), beam_robust.png.

### B2 (09:31-09:35) Pump dependence: same beam check at S.Pos2 9.46 deg (75.7 kHz, 24k acc) and 13.47 deg (177.9 kHz, 31k acc)
- 0931_agent_beam_70k, 0933_agent_beam_190k; analysis agent_scripts/pump_dependence.py (new, commit 40a05cb).
- Result: NO pump-dependent light anywhere on the full frame. Old spot (high-low) W-C -52 / -222 ADU/s, W-Woff -581 / +2,201 ADU/s;
  top DoG candidates z 5-7.5 all on the stray-lit post edge (negative/positive residuals of the scene). The stray scene gets
  DIMMER with pump: bright-pixel ratio 0.957 (76 vs 25 kHz) and 0.892 (178 vs 25 kHz), where SPDC light would be 3.0x / 7.2x (same
  as session 7: -7 %). At 178 kHz heralded W-Woff at the old spot: +1.3 +- 0.7 photons/readout vs ~5.6 expected.
- Expectation if the source were as on 21 Oct: ~123,000 ADU/s at 66 kHz -> ~300,000 ADU/s (~600 photons/s) at 178 kHz; measured
  < ~2,500 (2 sigma): > 100x below. Figure: 0933_agent_beam_190k/figures/pump_dependence.png (+ .json).

### C1 (09:37) Herald trigger health (threshold_rate.py, 0937_agent_thr_rate; at the minimum, gain 9000 x 3 ns gate)
- Rising 0.2/0.3/0.5/0.8/1.0/1.3 V: 30.3/29.8/29.4/29.7/29.6/29.8 kHz; falling: 29.8/29.7/29.9/29.4/29.3/28.7 kHz.
  Flat plateau on both edges: clean pulses > 1.3 V, no noise counting or double counting at 0.5 V.

### T1 (09:39-09:45) Heralded target at the minimum rate, 546 + 31 ns (target_image.py --acc 20000 --off-shift 1000, 0939_agent_target31_min)
- Hypothesis (Tristan's question): at ~20 kHz the per-trigger contrast should be at least as good as session 6 (z 33 at ~199 kHz).
- Settings: gain 9000, External ~21-25 k gates/s, 20,000 acc (0.62 ms open), OFF 1546+31 ns, blocks of 50, 200 ON + 200 OFF
  readouts, 0.80 (ON up to 0.94) s/readout.
- Result: net +0.06 +- 0.10 spot hits/readout, z 0.6 (ON -0.14, OFF -0.20; expected ~3.6 at 1.8e-4/trigger). ON and OFF hit
  images identical (same hot pixels). 2-sigma upper limit 1.4e-5 heralded spot hits per trigger, 13x below session 6.
- Figures: 0939_agent_target31_min/figures/heralded_single_photon_image.png, target_wiener.png, target_plot_run.png, target_threshold_sweep.png.

### C2 (offline) Camera sensitivity unchanged: stray-scene median W-C 150-250 ADU/px/s and the static point (232, 159) 12-22k ADU/s
  in session 5 (beam present), session 7 and today alike. The intensifier, External gating and readout work as when the photons
  were imaged; the missing light is upstream of the camera.

### Conclusion: the signal light (heralded or not) still does not reach the camera, at 25, 76 or 178 kHz, at any position on the
  sensor or any arrival time 27 ns - 10 us. Heralded imaging cannot proceed remotely; no narrower-gate / contrast-vs-rate
  comparison was possible. Stopped after ~27 min of camera time, as planned for this outcome.
### Left: S.Pos2 at 4.545 deg (minimum; 18.8 kHz at 09:47), group S status 12. Camera closed (checked free at 09:48).
### Lab checks (physical): (1) signal arm: is the signal fibre connected at both ends and coupled? Measure signal singles on an APD
  and signal-idler coincidences on the time tagger (expected delay ~550 ns herald-to-gate). (2) The herald side clearly changed
  since 21-23 Oct (rate at a fixed angle 200 -> 625 -> ~1000 kHz at 28 deg; extinction 10 % -> 0.8 %): what was changed (laser
  power, idler coupling/detector, PBS/HWP)? Check the idler counts are SPDC (coincidences with signal), not pump leakage or
  fluorescence. (3) Any cap/shutter/filter/block between the signal fibre output and the camera; did the imaging optics move.
  (4) Once the signal beam is visible again (beam_check.py: compact spot in W-C), re-run target_image.py at 546+31 ns at the minimum.

## 2026-10-26 session 9 (remote; after Tristan's APD check + signal-coupler touch-up + camera power cycle; goal: heralded image at the narrowest gate)
Data: `Desktop/HOI_freshstaret/data/troubleshoot/2026-10-26/` 1440_herald_rate_lock ... 1448_herald_rate_lock. Camera 14:42-14:49
(~7 min), free at start (now 192.168.0.81), -20 C Locked, no saturation, camera closed after every script. Only S.Pos2 moved (status 12).
Context (Tristan, APDs at the crystal-side couplers, max pump): 2 M idlers/s, 4 M signals/s, 200 k coincidences/s, accidentals ~8 k/s.

### W1 (14:42-14:43) Waveplate: found at 144.0 deg (moved by someone). Rate there 799.9 kHz (camera, gain 1, 3 ns gate; near the max).
- Moved to 102 deg (42 deg) and swept 1-deg steps down (1440_herald_rate_lock, ceiling 250 kHz): 102 103.4 | 101 93.2 | 100 86.6 |
  99 83.0 | 98 83.0 | 97 86.4 | 96 93.9 | 95 105.7 | 94 121.8 | 93 145 | 92 174 | 91 204 | 90 237 | 89 270 kHz (stopped).
- Fit: min 80.4 kHz at 8.60 deg (mod 90), max (extrapolated) 1.87 MHz at 53.6 deg. Parked at 98.596 deg.
- NEW: the minimum moved by +4 deg (4.6 -> 8.6 mod 90) and the floor rose from 15-30 kHz (session 8, this morning) to ~80 kHz.
  **The 15-20 kHz band is out of reach again**; the minimum itself is 4x above the 20 kHz trigger-error limit. Extinction ~4 %
  (80 / ~2000 kHz) vs 0.8 % this morning: the pump polarization into the HWP/PBS (or a pump-independent idler background) changed
  with today's lab work. Needs an ND filter / lower laser current, or a check of the polarization before the HWP.
- Figure: 1440_herald_rate_lock/figures/waveplate_calibration.png.

### B1 (14:43-14:45) Wide-gate beam check at the minimum (beam_check.py, 1443_agent_beam_min83k)
- Hypothesis: after the coupler realignment the signal beam reaches the camera again (APDs show 4 M signals/s at max pump).
- Settings: gain 9000, External 0.5 V rising, 24,000 acc; W = 27 ns + 10 us, Woff = 10027 ns + 10 us, C = 3 ns; 3 cycles x 10 readouts
  (0.57 / 0.83 / 0.30 s per readout). Herald 80.4 kHz at start, 55.9 kHz at end.
- Result: NO signal beam. Old spot (501, 454) core-annulus: W-C -139 +- 1,643, W-Woff +745 +- 1,838 ADU/s gate-open (~1.4 +- 3.5
  photons/s; 21 Oct at 66 kHz: +123,000). Full frame W-C: only the known stray scene (mount square, post (355, 454-504), static
  point (232, 159), present in Woff too and in sessions 7-8); W-Woff: noise, top DoG candidates z ~5 on the post edge.
  Stray-scene median W-C 176 ADU/px/s (this morning 182, 21 Oct 250): camera and imaging side unchanged.
- Heralded reading (W-Woff, 27 ns - 10 us): expected at session-6 efficiency 1.8e-4/trigger -> ~4.3 photons/readout = ~9,600 ADU/s;
  measured +745 +- 1,838 (~5 sigma below).
- Figure: 1443_agent_beam_min83k/figures/beam.png (+ beam.json).

### T1 (14:46-14:48) Heralded target at 546 + 31 ns (target_image.py --acc 45000 --cycles 2 --off-shift 1000, 1446_agent_target31_min83k)
- Rate right before / after: 76.1 / 80.3 kHz (98.596 deg). Run: ~78 k gates/s, 0.57 s/readout, 100 ON + 100 OFF readouts.
- Result: net +0.03 +- 0.13 spot hits/readout, z 0.3 (ON -0.04, OFF -0.07; expected ~8 at session-6 efficiency). ON and OFF
  hit images identical (same hot pixels). 2-sigma upper limit 6.4e-6 heralded spot hits per trigger (28x below session 6).
- Figures: 1446_agent_target31_min83k/figures/heralded_single_photon_image.png, target_wiener.png, target_plot_run.png, target_threshold_sweep.png.

### Herald rate this session: 800 kHz at 144 deg (found), 76-83 kHz at the minimum 98.6 deg, dipping to ~56 kHz for ~1 min (14:45).
  Imaging time above 20 kHz: 100 % (all ~7 min). Plots: 1440_herald_rate_lock/figures/herald_rate_vs_time_lock_session9.png/
  herald_rate_vs_time_lock.png (waveplate readings + angle), herald_rate_vs_time_runs_session9.png/herald_rate_vs_time.png (camera runs).

### Conclusion: after the realignment the signal light (heralded or not) still does not reach the camera: no spot anywhere on the full
  frame at 27 ns - 10 us, no heralded excess at 546+31 ns. Since the APDs see 4 M signals/s at the crystal-side coupler, the loss is
  between that coupler and the camera. No delay scan / narrow-gate series possible (nothing to narrow on). Stopped early as planned.
### Left: S.Pos2 at 98.596 deg (minimum, ~80 kHz), group S status 12. Camera closed (checked free 14:50).
### Lab checks (physical, in order):
  (1) Signal path after the APD test: is the signal fibre reconnected from the crystal-side coupler to the ~100 m delay fibre and from
      the delay fibre to the imaging mount (the APD test needed it unplugged)? Check each connector (FC/PC vs FC/APC mismatch, dust),
      and measure power/counts at the fibre end that goes into the imaging mount (e.g. APD on the delay-fibre output: expect ~MHz
      signals at max pump). That locates the loss to one segment.
  (2) At the imaging mount: fibre seated in the mount, caps off, filters/lens in place; any shutter/block between it and the camera.
  (3) Herald rate: the minimum is ~80 kHz now (15-30 kHz this morning), so the pump extinction got worse: check the pump
      polarization before the HWP, or add an ND filter / lower the laser current to reach 15-20 kHz at the minimum.
  (4) Then re-run beam_check.py (compact spot in W-C), then target_image.py at 546+31 ns and the narrow-gate series.

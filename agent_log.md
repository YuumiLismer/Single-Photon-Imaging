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

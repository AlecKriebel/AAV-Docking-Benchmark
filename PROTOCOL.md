# Glycan docking benchmark protocol v0.1

Prepared 24 October 2026. This is an operational pilot protocol and a proposed confirmatory design. Freeze the confirmatory version before scoring the test cohort. No performance threshold is a requirement that an algorithm must pass for publication: interpretable failures are results.

## Question and hypotheses

How do search information, receptor conformation, and glycan chemistry affect site recovery, pose accuracy, and experimentally measured specificity?

H1: Providing a true site improves pose recovery relative to predicted-site search; differences are larger at exposed interfaces. H2: Some failures reflect ranking of available near-native candidates, others insufficient returned sampling. H3: Size/charge-matched specificity comparisons reduce apparent performance driven by glycan size or sulfate count. H4: Removing neighboring subunits degrades interfacial-site predictions. These are hypotheses, not established outcomes.

## Minimum viable study

Development: AAV9 plus six unrelated, nonviral complexes, selected before docking to include monosaccharide, short linear, branched and sulfated chemistry where verified structures exist. AAV9 is a pipeline diagnostic, not a carbohydrate-torsion benchmark: a monosaccharide has no inter-sugar glycosidic linkage. AAV2 is an external site/evidence case, not an atomic RMSD reference based on low-resolution maps.

Confirmatory target: 24–36 curated complexes across at least 12 protein families, drawn from existing benchmarks; all eligible apo pairs within that cohort. This is a feasible starting target, not a power calculation. Estimate paired precision from development data and enlarge or narrow claims accordingly. Keep a later independent cohort if parameters or metrics change after seeing outcomes. Never handpick complexes because a method succeeded or failed.

Core methods: fpocket, P2Rank, AutoDock Vina, Vina-Carb, HADDOCK3. GNINA is initially a rescoring experiment on a common pose pool; audit training similarity before making generalization claims. GlycoTorch Vina is a justified Phase 4 addition because its GAG linkage/pucker support differs from Vina-Carb. No automatic expansion to many programs.

## Information regimes and equal comparisons

1. Full-surface/site-agnostic: ligand-free biological assembly or reference-independent symmetry sectors with complete neighbors. Candidate sites cannot use the native ligand or mutagenesis.
2. Predicted-site: top five unique predicted pockets, equal search allocation per pocket, same ligand inputs. Report site discovery separately from conditional docking.
3. True-site local: ligand-derived search box or true-interface restraints. This is an oracle diagnostic, not blind prediction.
4. Experimental-restraint: only independent prior evidence, with provenance and date. Avoid evaluating recovery solely against those same residues.

Cross holo versus apo with these regimes. Holo means ligand-stripped experimental bound protein. Apo means independently solved unbound protein, not a ligand-deleted holo structure. For HADDOCK, compare equally informed runs; do not compare its true-interface protocol directly with blind Vina. Report individual-pose ranks separately from cluster ranks.

For local docking, use ligand extent plus an 8 Å margin on each side, minimum side 20 Å, maximum 40 Å before review. AAV9 uses a fixed 20 Å cube. Search boxes can use the ligand only in the oracle arm. Predicted-site boxes use a predeclared ligand-size rule and predicted center. Preserve the same box for paired methods that support boxes. Guided restraints and boxes contain different information; label both rather than claiming exact equivalence.

## Structure and chemistry curation

Use mmCIF, with author and label chain/residue identifiers, model number, assembly transforms and revision history. Materialize the biological assembly once; do not expand an already complete capsid twice. Match apo/holo sequences explicitly, verify residue identity, align receptor backbone, retain transformations. Distinguish missing loops, alternate conformers, ions, covalent glycans, crystallographic contacts, and assembly-dependent recognition.

Primary pose references must have a chemically defined, noncovalent ligand with complete scored atoms and interpretable local density. Nominal resolution alone is insufficient. Review wwPDB validation, local density/occupancy, anomer, stereochemistry, ring puckering and glycosidic linkages. Reject covalent enzyme intermediates from noncovalent docking. Preserve structural waters/metals known to support recognition; specify a dry baseline and paired retained-water condition, not arbitrary deletion. Never treat unresolved glycan tails as exact coordinates.

Build ligands from CCD and explicit connectivity, verified against deposited atoms and glycan notation (GlycoCT/WURCS/IUPAC where available). Independently generate input conformers. Record anomer, reducing-end modification, linkage positions, sulfate positions, formal charge, protonation, total heavy atoms, glycosidic torsions and rotatable bonds. Record pH/ionic strength. No unsupported ligand parameters: unsupported chemistry is a reported failure category, not a silent exclusion.

Receptor preparation must retain all relevant subunits and report every atom addition, deletion or repaired bond. Freeze protonation rules; inspect histidine and acidic/basic residues locally. Run monomer versus complete-neighborhood ablations only as explicitly incomplete-model controls. Check the retained neighborhood covers the entire box plus interaction cutoff; evaluate final poses against the full assembly for clashes/accessibility. The pilot's reference-based neighborhood is not eligible for site-discovery claims.

## Search and compute

Pilot Vina: 1.2.7, standard Vina score, rigid receptor, Meeko 0.8.0, three seeds 104729/130363/155921, exhaustiveness 16, two CPUs, 20 returned poses, 0.5 Å minimum pose separation, 10 kcal/mol output window. Ligand embedding seed 20261024. No explicit waters. These settings are recorded in code/results, not tuned to the observed outcome.

Confirmatory initial proposal: same seeds, exhaustiveness 32, retain up to 50 poses; run a prespecified 16/32/64 convergence study only on development systems. For Vina-Carb reproduce its supported input conventions and pin its actual source commit; do not assume a modern Vina wrapper can invoke it. For HADDOCK reuse the authors' configurations and report all generated/refined counts, restraint definitions and native evaluation data access. First reproduce one released example. Software-specific search units are not equivalent: report both a standardized configuration comparison and a matched CPU-time sensitivity analysis. Log actual wall/CPU time, memory, hardware, threads, failures and retry costs.

## Endpoints

Primary site endpoint: Top-1/Top-5 recovery of distinct site classes, with predicted center within 4 Å of any native ligand heavy atom. Report 3/5 Å sensitivity, center-to-ligand distance, residue overlap and candidate coverage. A point near one tip of a long GAG is not complete-interface recovery; also report fraction of native interface residues covered. Merge symmetry-equivalent predictions before ranking. Baselines: random solvent-accessible surface patches and matched-size/charge patches; random points in a bounding box are unsuitable.

Primary pose endpoint: heavy-atom ligand RMSD ≤2 Å among Top-1 and Top-5; report ≤1/3 Å sensitivity and continuous RMSD. Evaluate both global-backbone and local-interface-backbone alignment; the pilot currently uses the fixed holo coordinate frame, not a separate local-interface fit. Respect chemical identity and graph/stereochemistry-valid atom permutations. Never independently align the ligand, which would remove pose error. Minimize over valid symmetry-equivalent sites only after consistent receptor transforms; do not count 60 capsid sites as 60 observations. Long glycans additionally require residue-wise/core and whole-ligand metrics, with the core defined from density before docking.

Contacts: heavy-atom pairs <4 Å; report recall, precision and residue-level overlap. Use the holo native contact identities for apo evaluation. Contact recall alone is inadequate for small glycans. Track steric clashes, ring geometry and linkage outliers independently of RMSD.

Sampling/ranking: record best returned RMSD, Top-1/5, rank of the closest pose, and success conditional on a near-native pose being present. Returned poses are a score-filtered subset; absent near-native returned poses do not prove that none were generated internally. A definitive sampling/scoring decomposition requires retained pre-ranking candidates or instrumentation. Native-pose scoring/local optimization is a diagnostic, never a docking success.

Specificity: define positives, tested negatives and untested alternatives separately. An unlisted array glycan is not a negative. Restrict each ranking to one assay/protein context, including glycan density/linker, concentration, valency, salt, pH and detection limit. Use AUROC and AUPRC only with tested labels and both classes; report prevalence. Use per-target Spearman/Kendall correlation for comparable quantitative series, and concordance with censoring where needed. Do not pool array RFU, KD and transduction. Do not infer KD from Vina score. Report raw and size/charge-matched rankings; do not normalize away a genuine mechanism without showing both.

## Statistics and leakage

Experimental unit is an independent protein family; runs, poses, symmetry copies and closely related structures are nested measurements. Report per-complex and family-macro results. Calculate paired differences and 95% family-cluster bootstrap intervals (10,000 resamples), averaging seeds within systems first. Resample whole families retaining algorithm pairing. Where data permit, use family-aware logistic/mixed models with algorithm, information regime, receptor state, heavy-atom count, charge, torsion count and exposure/geometry. Avoid a high-dimensional interaction model on 24 cases: limit prespecified covariates or treat associations descriptively. Structural-quality covariates include local uncertainty, missing residues and apo–holo site displacement. Associations do not establish mechanistic causation; paired ablations provide stronger evidence.

Cluster sequences and structural folds before splitting; report thresholds and cross-family exceptions. Audit P2Rank/GNINA training sets, related homologs and ligand overlaps; mark unknown provenance. Reference-based neighborhoods, water selection or restraints are leakage in a blind arm. Use reference structures only in scoring after test inputs are frozen. Deduplicate overlaps among reused datasets.

## Decision gates

Gate 1: reference chemistry, assemblies and atom maps pass independent checks; preparation errors resolved or explicitly excluded. Gate 2: positive metric tests and native-coordinate controls pass; six development systems execute without silent atom loss. Gate 3: sufficient distinct families and eligible negatives exist to support the intended endpoints. If no reliable specificity data, publish a site/pose study and defer specificity. Gate 4: changes survive uncertainty intervals and sensitivity to preparation, information regime and family weighting; otherwise emphasize scope/limitations rather than an algorithm leaderboard. Gate 5: external reviewer can reproduce one complex from archived inputs, versions and one command.

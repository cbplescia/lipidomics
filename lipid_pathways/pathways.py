"""Pathway definitions. Node: id -> (col,row,text,kind[s=this lipid,m=metabolite,x=other pathway,r=recycled]).
Edge: (src,dst,enzyme label,curve,kind[n|r])."""
K = {"s": "self", "m": "met", "x": "ext", "r": "rec"}


def build(nodes, edges):
    N = {k: (c, r, t, K[kd]) for k, (c, r, t, kd) in nodes.items()}
    E = []
    for e in edges:
        a, b, lab = e[0], e[1], e[2]
        rad = e[3] if len(e) > 3 else 0
        kind = e[4] if len(e) > 4 else ("r" if N[b][3] == "rec" else "n")
        E.append((a, b, lab, rad, kind))
    return N, E


P = {}

P["PC"] = dict(
    title="PC – Phosphatidylcholine", sub="Most abundant membrane phospholipid · ER / Golgi · Kennedy (CDP-choline) pathway",
    nodes={
        "cho": (0, 0, "Choline", "m"), "pcho": (0, 1, "Phosphocholine", "m"), "cdpc": (0, 2, "CDP-choline", "m"),
        "dg": (0, 3, "Diacylglycerol (DG)", "m"), "pe": (0, 4, "PE", "x"),
        "pc": (1, 1.5, "PC", "s"), "lpc": (1, 3.5, "LPC", "m"), "sm": (1, 0, "Sphingomyelin", "m"),
        "pa": (2, 0, "PA + choline", "m"), "dgp": (2, 1.5, "DG + phosphocholine", "m"),
        "gpc": (2, 3.5, "GPC + free FA", "m"), "fa": (2, 2.5, "Free FA / lyso-PC pool", "x"),
        "r1": (3, 0, "↺ Choline → CHKA\n(re-enters Kennedy)", "r"), "r2": (3, 1.5, "↺ DG → DGAT (TG)\nor CEPT/CHPT (PC)", "r"),
        "r3": (3, 3, "↺ FA → acyl-CoA\nLPCAT / β-oxidation", "r"), "r4": (3, 4, "↺ Glycerol-3-P → GPAT\n/ glycolysis", "r"),
    },
    edges=[("cho", "pcho", "CHKA/B"), ("pcho", "cdpc", "PCYT1 (CCT)\nrate-limiting"), ("cdpc", "pc", "CHPT1 / CEPT1", 0.0),
           ("dg", "pc", "(+ CDP-choline)", 0.0), ("pe", "pc", "PEMT (liver)\n3× SAM methylation", 0.0),
           ("pc", "lpc", "PLA2 (−FA)", 0.25), ("lpc", "pc", "LPCAT1–4\n+acyl-CoA", 0.25), ("pc", "sm", "SMS1/2\n(Cer → DG)"),
           ("pc", "pa", "PLD1/2", 0.0), ("pc", "dgp", "PLC"), ("lpc", "gpc", "lysophospholipase\n(PNPLA6)"),
           ("pa", "r1", "PLD product"), ("dgp", "r2", ""), ("gpc", "r3", ""), ("gpc", "r4", "GDPD5\n(+ choline)"),
           ("dgp", "r1", "", -0.15)],
)

P["PE"] = dict(
    title="PE – Phosphatidylethanolamine", sub="Inner-leaflet / mitochondrial & ER lipid · Kennedy (CDP-ethanolamine) and PS-decarboxylase routes",
    nodes={
        "eth": (0, 0, "Ethanolamine", "m"), "pet": (0, 1, "Phosphoethanolamine", "m"), "cdpe": (0, 2, "CDP-ethanolamine", "m"),
        "dg": (0, 3, "Diacylglycerol (DG)", "m"), "ps": (0, 4, "PS (mitochondria)", "x"),
        "pe": (1, 1.5, "PE", "s"), "pc": (1, 0, "PC", "m"), "lpe": (1, 3, "LPE", "m"),
        "nape": (1, 4.2, "NAPE → anandamide\n(NAPE-PLD)", "m"), "lc3": (1, 5.2, "PE-LC3/ATG8\n(autophagy)", "m"),
        "pa": (2, 0.3, "PA + ethanolamine", "m"), "dgp": (2, 1.5, "DG + phospho-\nethanolamine", "m"),
        "gpe": (2, 3, "GPE + free FA", "m"),
        "r1": (3, 0.3, "↺ Ethanolamine → ETNK\n(Kennedy)", "r"), "r2": (3, 1.5, "↺ P-Etn → PCYT2\n→ CDP-Etn → PE", "r"),
        "r3": (3, 3, "↺ FA → acyl-CoA\nLPEAT / β-oxidation", "r"), "r4": (3, 4.2, "↺ Glycerol-3-P → GPAT", "r"),
    },
    edges=[("eth", "pet", "ETNK1/2"), ("pet", "cdpe", "PCYT2\nrate-limiting"), ("cdpe", "pe", "CEPT1 / EPT1\n(+ DG)"),
           ("dg", "pe", ""), ("ps", "pe", "PISD\n(−CO₂)"),
           ("pe", "pc", "PEMT (liver)\n3× SAM"), ("pe", "lpe", "PLA2 (−FA)", 0.2), ("lpe", "pe", "LPEAT\n(MBOAT1/2)", 0.2),
           ("pe", "nape", "N-acylation\n(PLA2G4E)"), ("pe", "lc3", "ATG3 conjugation"),
           ("pe", "pa", "PLD"), ("pe", "dgp", "PLC"), ("lpe", "gpe", "lysophospholipase"),
           ("pa", "r1", ""), ("dgp", "r2", ""), ("gpe", "r3", ""), ("gpe", "r4", "GDPD\n(+ ethanolamine)"),
           ],
)

P["TG"] = dict(
    title="TG – Triacylglycerol", sub="Neutral storage lipid · ER → lipid droplets (LD) · VLDL / chylomicron export",
    nodes={
        "g3p": (0, 0, "Glycerol-3-P + FA-CoA", "m"), "lpa": (0, 1, "Lyso-PA (LPA)", "m"), "pa": (0, 2, "Phosphatidic acid", "m"),
        "dg": (0, 3, "Diacylglycerol (DG)", "m"), "mg": (0, 4, "2-MG (intestine)", "x"),
        "tg": (1, 1.7, "TG", "s"), "ld": (1, 0.2, "Lipid droplet storage\n(PLIN, CIDE)", "m"),
        "vldl": (1, 3.5, "VLDL / chylomicron\nsecretion (MTTP, apoB)", "m"),
        "dg2": (2, 0.2, "DG + FA", "m"), "mg2": (2, 1.7, "MG + FA", "m"), "gly": (2, 3.2, "Glycerol + FA", "m"),
        "lpl": (2, 4.5, "Plasma TG → FA\n(LPL, in tissues)", "x"), "lal": (2, 5.5, "Lipophagy: lysosomal\nLIPA → FA + glycerol", "x"),
        "r1": (3, 0.2, "↺ DG → TG\nor PC / PE", "r"), "r2": (3, 1.7, "↺ FA → ACSL → acyl-CoA\n→ β-oxidation or re-esterify", "r"),
        "r3": (3, 3.2, "↺ Glycerol → GK\n→ G3P (liver)", "r"),
    },
    edges=[("g3p", "lpa", "GPAT1-4"), ("lpa", "pa", "AGPAT1-5"), ("pa", "dg", "LIPIN1-3\n(PAP)"), ("dg", "tg", "DGAT1/2\n(+ acyl-CoA)"),
           ("mg", "dg", "MGAT2"), ("tg", "ld", "LD budding\n(SEIPIN)"), ("tg", "vldl", "MTTP lipidation", 0.0),
           ("tg", "dg2", "ATGL (PNPLA2)\n+ CGI-58"), ("dg2", "mg2", "HSL (LIPE)"), ("mg2", "gly", "MGL (MGLL)"),
           ("vldl", "lpl", ""), ("vldl", "lal", "", -0.1),
           ("dg2", "r1", "", 0.0), ("mg2", "r2", ""), ("gly", "r3", "glycerol kinase"), ("gly", "r2", "", 0.0)],
)

P["CL"] = dict(
    title="CL – Cardiolipin", sub="Mitochondrial inner-membrane signature lipid · synthesis, tafazzin remodeling, oxidation / mitophagy signalling",
    nodes={
        "pa": (0, 0, "PA (imported from ER)\nPRELID1–TRIAP1", "m"), "cdp": (0, 1, "CDP-DAG", "m"), "pgp": (0, 2, "PGP", "m"),
        "pg": (0, 3, "PG", "m"), "cdp2": (0, 4, "CDP-DAG (2nd)", "m"),
        "cl": (1, 1.5, "Nascent CL", "s"), "mlcl": (1, 3.2, "MLCL\n(monolyso-CL)", "m"), "mat": (1, 0, "Mature CL (tetra-18:2)", "s"),
        "ox": (1, 4.8, "Oxidised CL / externalised CL\n(mitophagy, apoptosis)", "m"),
        "dlcl": (2, 1.5, "DLCL (dilyso-CL) + FA", "m"), "pa2": (2, 3.2, "PA (mito) + DG", "m"),
        "gpg": (2, 0, "GPG + FA", "m"),
        "r1": (3, 0, "↺ Glycerol / G3P pool", "r"), "r2": (3, 1.5, "↺ FA → CoA → β-oxidation\nor re-acylation (ALCAT1)", "r"),
        "r3": (3, 3.2, "↺ PA → CDP-DAG\n→ PG / CL (mito)", "r"), "r4": (3, 4.8, "↺ Mitophagy → lysosome\n→ BMP / FA recycle", "r"),
    },
    edges=[("pa", "cdp", "TAMM41"), ("cdp", "pgp", "PGS1"), ("pgp", "pg", "PTPMT1"), ("pg", "cl", "CRLS1\n(+ CDP-DAG)"),
           ("cdp2", "cl", "", -0.0),
           ("cl", "mlcl", "iPLA2γ / PLA2G6\n(−FA)", 0.2), ("mlcl", "cl", "", 0.2), ("cl", "mat", "Tafazzin (TAZ)\nMLCLAT1 / ALCAT1"),
           ("mat", "ox", "Cyt c peroxidase\n(ROS)", 0.0), ("mlcl", "dlcl", "PLA2"),
           ("cl", "pa2", "PLD6 (MitoPLD)"), ("dlcl", "gpg", "PLA2"),
           ("gpg", "r1", ""), ("dlcl", "r2", ""), ("pa2", "r3", ""), ("ox", "r4", "")],
)

P["SM"] = dict(
    title="SM – Sphingomyelin", sub="ER (ceramide) → Golgi / plasma membrane (SM) · raft lipid · lysosomal degradation · salvage",
    nodes={
        "ser": (0, 0, "Serine + Palmitoyl-CoA", "m"), "kds": (0, 1, "3-Ketosphinganine →\nSphinganine", "m"), "dhc": (0, 2, "Dihydroceramide", "m"),
        "cer": (0, 3, "Ceramide (Cer)", "m"), "pc": (0, 4.2, "PC (phosphocholine donor)", "x"),
        "sm": (1, 2, "SM", "s"), "raft": (1, 0.3, "Lipid rafts / PM outer\nleaflet; cholesterol", "m"),
        "cer2": (2, 0.5, "Ceramide + phosphocholine", "m"), "sph": (2, 2.5, "Sphingosine + FA", "m"),
        "s1p": (2, 3.7, "S1P", "m"), "lyase": (2, 4.9, "Hexadecenal + P-Etn", "m"),
        "r1": (3, 0.5, "↺ Cer → SM / GlcCer", "r"), "r2": (3, 2.5, "↺ Sph → CerS → Cer\n(salvage)", "r"),
        "r3": (3, 3.7, "↺ Phosphocholine → CDP-choline", "r"), "r4": (3, 4.9, "↺ P-Etn → PE\nhexadecenal → palmitoyl-CoA", "r"),
    },
    edges=[("ser", "kds", "SPT, KDSR"), ("kds", "dhc", "CerS1-6"), ("dhc", "cer", "DEGS1"), ("cer", "sm", "SMS1 (Golgi)\nSMS2 (PM)"),
           ("pc", "sm", "+ PC → DG", 0.0), ("sm", "raft", ""),
           ("sm", "cer2", "Acid SMase (SMPD1)\nNeutral SMase (SMPD2/3)\nAlk SMase (ENPP7)"), ("cer2", "sph", "Ceramidases\nASAH1 / ASAH2 / ACER"),
           ("sph", "s1p", "SphK1/2"), ("s1p", "lyase", "S1P lyase (SGPL1)"),
           ("cer2", "r1", ""), ("sph", "r2", ""), ("cer2", "r3", "", -0.2), ("lyase", "r4", ""), ("s1p", "r2", "SGPP1/2", 0.0)],
)

P["LBPA"] = dict(
    title="LBPA / BMP – Bis(monoacylglycero)phosphate", sub="Late endosome / lysosome lipid (intralumenal vesicles) · unusual S,S stereochemistry · acid-hydrolase cofactor",
    nodes={
        "pg": (0, 0, "PG (from PA → CDP-DAG\n→ PGP → PG)", "m"), "lpg": (0, 1.8, "LPG\n(lyso-phosphatidylglycerol)", "m"),
        "pg2": (0, 3.6, "PG (2nd acyl donor)", "m"),
        "bmp": (1, 1.8, "BMP / LBPA", "s"), "ilv": (1, 0.1, "Intralumenal vesicles (ILV)\nacidic, PI3P → PI(3,5)P2", "m"),
        "fn": (1, 3.6, "Functions: activates lysosomal\nlipases, GBA1, NPC2 cholesterol export", "m"),
        "lbpa2": (2, 1.8, "Lyso-BMP + FA", "m"), "gpg": (2, 3.6, "GPG (glycerophospho-\nglycerol) + FA", "m"),
        "exo": (2, 0.1, "Exosome release\n(ILV → extracellular)", "x"),
        "r1": (3, 1.8, "↺ FA → acyl-CoA\n(lysosomal efflux)", "r"), "r2": (3, 3.6, "↺ Glycerol / G3P pool\n(→ glycerolipids)", "r"),
        "r3": (3, 0.1, "↺ Re-synthesised from PG", "r"),
    },
    edges=[("pg", "lpg", "PLA2G15 / PLA2\n(−FA)"), ("lpg", "bmp", "CLN5 (BMP synthase)\ntransacylation", 0.0),
           ("pg2", "bmp", "PLD3/4 (BMP synthase)", 0.0), ("bmp", "ilv", ""), ("bmp", "fn", ""),
           ("bmp", "lbpa2", "PLA2G15 (LPLA2)\n(slow – BMP is resistant)"), ("lbpa2", "gpg", "PLA2 / LYPLA"),
           ("ilv", "exo", ""), ("lbpa2", "r1", ""), ("gpg", "r2", "GDE / lysosomal\nglycerophosphodiesterase"),
           ("exo", "r3", "", 0.0)],
)

P["Cer"] = dict(
    title="Cer – Ceramide", sub="Central sphingolipid hub · de novo (ER) / sphingomyelinase / salvage routes",
    nodes={
        "ser": (0, 0, "Serine + Palmitoyl-CoA", "m"), "kds": (0, 1, "3-Ketosphinganine →\nSphinganine", "m"), "dhc": (0, 2, "Dihydroceramide", "m"),
        "sm": (0, 3.2, "SM  (SMase route)", "x"), "gc": (0, 4.2, "GlcCer (GBA1 route)", "x"), "sph": (0, 5.2, "Sphingosine (salvage)", "x"),
        "cer": (1, 2, "Ceramide", "s"), "sm2": (1, 0.2, "SM", "m"), "glc": (1, 1.2, "GlcCer / GalCer", "m"),
        "c1p": (1, 3.2, "Cer-1-phosphate", "m"), "acyl": (1, 4.4, "Acylceramide (skin)", "m"),
        "sph2": (2, 1.5, "Sphingosine + FA", "m"), "s1p": (2, 3, "S1P", "m"), "ly": (2, 4.5, "Hexadecenal + P-Etn", "m"),
        "r1": (3, 1.5, "↺ Sph → CerS → Cer\n(salvage, ~50–90 %)", "r"), "r2": (3, 3, "↺ SGPP → Sph", "r"),
        "r3": (3, 4.5, "↺ P-Etn → PE;\nhexadecenal → palmitoyl-CoA", "r"),
    },
    edges=[("ser", "kds", "SPT, KDSR"), ("kds", "dhc", "CerS1-6"), ("dhc", "cer", "DEGS1"), ("sm", "cer", "SMases", 0.0),
           ("gc", "cer", "GBA1 / GBA2"), ("sph", "cer", "CerS", 0.0),
           ("cer", "sm2", "SMS1/2"), ("cer", "glc", "UGCG / UGT8"), ("cer", "c1p", "CERK"), ("cer", "acyl", "CERS3 / ω-OH"),
           ("cer", "sph2", "Ceramidases\nASAH1 / ASAH2 / ACER1-3"), ("sph2", "s1p", "SphK1/2"), ("s1p", "ly", "SGPL1"),
           ("sph2", "r1", ""), ("s1p", "r2", ""), ("ly", "r3", "")],
)

P["PS"] = dict(
    title="PS – Phosphatidylserine", sub="Anionic inner-leaflet lipid · ER-MAM synthesis · PS decarboxylase in mitochondria · apoptotic “eat-me” signal",
    nodes={
        "pc": (0, 0, "PC + Serine", "m"), "pe": (0, 1.3, "PE + Serine", "m"),
        "ps": (1, 0.65, "PS", "s"),
        "flip": (1, 2.4, "Inner-leaflet asymmetry\n(flippases ATP11A/C)", "m"), "scr": (1, 3.6, "Apoptosis: scramblases\nTMEM16F / XKR8 → eat-me", "m"),
        "pe2": (1, 5, "PE (mitochondria)", "m"),
        "lps": (2, 0.65, "LysoPS + FA", "m"), "gps": (2, 2.4, "GPS + FA", "m"), "gl": (2, 3.6, "Glycerol-3-P + Serine", "m"),
        "r1": (3, 0.65, "↺ FA → acyl-CoA", "r"), "r2": (3, 2.4, "↺ Re-acylated by LPSAT\n(MBOAT5) → PS", "r"),
        "r3": (3, 3.6, "↺ Serine → PTDSS / SPT", "r"), "r4": (3, 5, "↺ PE → PC (PEMT)\nor Kennedy", "r"),
    },
    edges=[("pc", "ps", "PTDSS1\n(base exchange, −choline)"), ("pe", "ps", "PTDSS2\n(−ethanolamine)"), ("ps", "flip", ""),
           ("ps", "scr", ""), ("ps", "pe2", "PISD\n(−CO₂)", 0.0),
           ("ps", "lps", "ABHD16A / PLA1A\n(PLA1)"), ("lps", "gps", "ABHD12"), ("gps", "gl", "GDE / GDPD"),
           ("lps", "r1", ""), ("lps", "r2", "", 0.0), ("gl", "r3", ""), ("pe2", "r4", "")],
)

P["DG"] = dict(
    title="DG – Diacylglycerol", sub="Metabolic hub and second messenger · PKC / Munc13 / RasGRP signalling",
    nodes={
        "pa": (0, 0, "Phosphatidic acid (PA)", "m"), "pip2": (0, 1.2, "PI(4,5)P₂ → IP₃", "m"), "pc": (0, 2.4, "PC / SM (PLC, SMS)", "m"),
        "tg": (0, 3.6, "TG (lipolysis)", "m"), "mg": (0, 4.8, "2-MG + acyl-CoA", "x"),
        "dg": (1, 2.2, "DG", "s"), "sig": (1, 0.2, "Signalling: PKC, Munc13, RasGRP", "m"),
        "tg2": (1, 4.2, "TG", "m"),
        "pa2": (2, 0.2, "PA", "m"), "ag": (2, 1.7, "2-AG + FA", "m"), "mg2": (2, 3.2, "MG + FA", "m"),
        "pcpe": (2, 4.7, "PC / PE", "m"),
        "r1": (3, 0.2, "↺ PA → CDP-DAG → PI\n(PI cycle) / PG / CL", "r"), "r2": (3, 1.7, "↺ AA + glycerol\n(MAGL; eicosanoids)", "r"),
        "r3": (3, 3.2, "↺ Glycerol → G3P (GK)\nFA → acyl-CoA", "r"), "r4": (3, 4.7, "↺ Membrane phospholipids", "r"),
    },
    edges=[("pa", "dg", "LIPIN1-3 (PAP)"), ("pip2", "dg", "PLCβ/γ/δ"), ("pc", "dg", "PLC / SMS"), ("tg", "dg", "ATGL"), ("mg", "dg", "MGAT"),
           ("dg", "sig", ""), ("dg", "tg2", "DGAT1/2"),
           ("dg", "pa2", "DGKα–ι\n(+ ATP)"), ("dg", "ag", "DAGLα/β\n(sn-1 lipase)"), ("dg", "mg2", "HSL (LIPE)"), ("dg", "pcpe", "CEPT1 / CHPT1"),
           ("pa2", "r1", ""), ("ag", "r2", "MAGL"), ("mg2", "r3", "MGL"), ("pcpe", "r4", "")],
)

P["PG"] = dict(
    title="PG – Phosphatidylglycerol", sub="Mitochondrial inner membrane & lung surfactant lipid · precursor of CL and BMP",
    nodes={
        "pa": (0, 0.3, "PA", "m"), "cdp": (0, 1.8, "CDP-DAG", "m"), "pgp": (0, 3.3, "PGP\n(+ glycerol-3-P)", "m"),
        "pg": (1, 1.8, "PG", "s"), "cl": (1, 0.1, "Cardiolipin (CL)", "m"), "bmp": (1, 3.5, "BMP / LBPA\n(late endosome)", "m"),
        "sur": (1, 4.8, "Lung surfactant\n(alveolar PG)", "m"),
        "lpg": (2, 0.1, "LPG + FA", "m"), "pad": (2, 1.8, "PA + glycerol", "m"), "gpg": (2, 3.5, "GPG + FA", "m"),
        "r1": (3, 0.1, "↺ LPG re-acylated\n(LPGAT1) → PG", "r"), "r2": (3, 1.8, "↺ PA → LIPIN → DG\nor CDP-DAG", "r"),
        "r3": (3, 3.5, "↺ Glycerol / G3P;\nFA → acyl-CoA", "r"),
    },
    edges=[("pa", "cdp", "TAMM41 / CDS"), ("cdp", "pgp", "PGS1"), ("pgp", "pg", "PTPMT1"),
           ("pg", "cl", "CRLS1\n(+ CDP-DAG)"), ("pg", "bmp", "PLA2G15 → CLN5 /\nPLD3/4"), ("pg", "sur", "secretion (type II cells)"),
           ("pg", "lpg", "PLA2\n(mito iPLA2γ)", 0.0), ("pg", "pad", "PLD"), ("lpg", "gpg", "lysophospholipase"),
           ("lpg", "r1", ""), ("pad", "r2", ""), ("gpg", "r3", "")],
)

P["ChE"] = dict(
    title="ChE – Cholesteryl Ester", sub="Neutral storage / transport form of cholesterol · ER, lipid droplets, lipoproteins",
    nodes={
        "chol": (0, 0, "Free cholesterol\n(ER / lysosome exit)", "m"), "acoa": (0, 1.3, "Acyl-CoA (FA-CoA)", "m"),
        "pc": (0, 2.6, "HDL PC (plasma)", "x"), "ldl": (0, 3.9, "LDL / HDL uptake\n(LDLR, SR-B1)", "x"),
        "ce": (1, 1.2, "ChE", "s"), "ld": (1, 0, "Lipid droplets", "m"), "lp": (1, 2.6, "Lipoprotein cores\n(HDL ↔ VLDL/LDL via CETP)", "m"),
        "foam": (1, 4, "Foam cells (excess ChE\nin macrophages)", "m"),
        "ch2": (2, 0, "Cholesterol + FA\n(cytosolic)", "m"), "ch3": (2, 1.6, "Cholesterol + FA\n(lysosome)", "m"),
        "r1": (3, 0, "↺ Re-esterified (ACAT1/2)", "r"), "r2": (3, 1.6, "↺ ABCA1/G1 efflux → HDL\n→ liver (RCT)", "r"),
        "r3": (3, 3.2, "↺ Bile acids (CYP7A1);\nsteroid hormones", "r"), "r4": (3, 4.6, "↺ FA → acyl-CoA\n→ β-oxidation", "r"),
    },
    edges=[("chol", "ce", "ACAT1/SOAT1 (ER)\nACAT2 (liver/gut)"), ("acoa", "ce", ""), ("pc", "ce", "LCAT (+ cholesterol\n→ LPC)"),
           ("ldl", "ce", "", 0.0), ("ce", "ld", ""), ("ce", "lp", "CETP transfer"), ("ce", "foam", "", 0.0),
           ("ld", "ch2", "NCEH1 / CES1 / HSL"), ("ce", "ch3", "LIPA (LAL)\nlysosomal acid lipase"),
           ("ch2", "r1", ""), ("ch3", "r2", ""), ("ch3", "r3", ""), ("ch3", "r4", "")],
)

P["AcCa"] = dict(
    title="AcCa – Acylcarnitine", sub="Carnitine shuttle for fatty acid entry into mitochondria · peroxisomal / short-chain buffering",
    nodes={
        "fa": (0, 0, "Fatty acid", "m"), "acoa": (0, 1.3, "Acyl-CoA", "m"), "car": (0, 2.6, "L-Carnitine\n(diet; lys/met; OCTN2)", "m"),
        "acca": (1, 1.3, "Acylcarnitine (cytosol)", "s"), "mito": (1, 2.8, "Acylcarnitine (matrix)", "s"),
        "perox": (1, 4.2, "Peroxisome: CROT (medium)\nCRAT (short chain)", "m"), "export": (1, 0, "Efflux to plasma/\nurine (OCTN, OCTN2)", "x"),
        "acoa2": (2, 2.8, "Acyl-CoA + carnitine", "m"), "boxn": (2, 1.3, "β-oxidation", "m"), "ac": (2, 4.2, "Acetyl-CoA\n(TCA, ketones)", "m"),
        "r1": (3, 0.1, "↺ Urinary / bile excretion\n(diagnostic biomarker)", "r"), "r2": (3, 1.4, "↺ Carnitine → CACT → cytosol\n→ CPT1 again", "r"),
        "r3": (3, 2.8, "↺ Acetyl-CoA → ATP", "r"), "r4": (3, 4.2, "↺ Acetylcarnitine buffer\n(CRAT, CoA pool)", "r"),
    },
    edges=[("fa", "acoa", "ACSL1/ACSVL"), ("acoa", "acca", "CPT1A/B/C\n(outer mito. membr.)"), ("car", "acca", "(+ carnitine)"),
           ("acca", "mito", "CACT (SLC25A20)\ntranslocase"), ("mito", "perox", "", 0.0), ("acca", "export", "", 0.0),
           ("mito", "acoa2", "CPT2\n(inner mito. membr.)"), ("acoa2", "boxn", "", 0.0), ("boxn", "ac", "", 0.0),
           ("acoa2", "r2", "carnitine return"), ("boxn", "r3", ""), ("perox", "r4", "CRAT"), ("export", "r1", "")],
)

P["LPC"] = dict(
    title="LPC – Lysophosphatidylcholine", sub="Plasma transport lipid (albumin) · signalling · LPC-DHA brain uptake via MFSD2A",
    nodes={
        "pc": (0, 0.8, "PC", "m"), "hdl": (0, 2.4, "HDL: PC + cholesterol", "x"), "el": (0, 3.8, "Dietary / membrane PC", "x"),
        "lpc": (1, 1.6, "LPC", "s"), "alb": (1, 0.1, "Albumin-bound transport\n(plasma)", "m"), "mf": (1, 3.2, "MFSD2A: LPC-DHA\nacross BBB", "m"),
        "lpa": (1, 4.4, "LPA + choline", "m"),
        "gpc": (2, 1.6, "GPC + FA", "m"), "gpcc": (2, 3.2, "Glycerol-3-P + choline", "m"),
        "pc2": (2, 0.1, "PC (re-acylation)", "m"),
        "r1": (3, 1.6, "↺ FA → acyl-CoA", "r"), "r2": (3, 3.2, "↺ Choline → CHKA → Kennedy;\nG3P → glycolysis", "r"),
        "r3": (3, 0.1, "↺ PC re-enters membranes", "r"), "r4": (3, 4.4, "↺ LPA → receptors / LPP →\nMG → glycerolipids", "r"),
    },
    edges=[("pc", "lpc", "PLA2 (cPLA2α, sPLA2)"), ("hdl", "lpc", "LCAT (+ ChE)", 0.0), ("el", "lpc", "endothelial lipase\n(PLA1)", 0.0),
           ("lpc", "alb", ""), ("lpc", "mf", ""), ("lpc", "lpa", "Autotaxin (ENPP2)"),
           ("lpc", "pc2", "LPCAT1–4\n+acyl-CoA", 0.0), ("lpc", "gpc", "Lysophospholipase\n(PNPLA6, LYPLA1/2)"), ("gpc", "gpcc", "GDPD5"),
           ("gpc", "r1", ""), ("gpcc", "r2", ""), ("pc2", "r3", ""), ("lpa", "r4", "", 0.0)],
)

P["Hex1Cer"] = dict(
    title="Hex1Cer – Hexosylceramide (GlcCer / GalCer)", sub="Golgi-synthesised glycosphingolipid · precursor of complex gangliosides · lysosomal GBA1 degradation",
    nodes={
        "cer": (0, 0.5, "Ceramide (ER)", "m"), "udp": (0, 1.8, "UDP-glucose\nUDP-galactose", "m"), "cerg": (0, 3.2, "Ceramide (Golgi)\nvia CERT / vesicles", "m"),
        "hex": (1, 1.5, "Hex1Cer\n(GlcCer / GalCer)", "s"), "lac": (1, 0, "LacCer → GM3 → gangliosides", "m"),
        "sulf": (1, 3, "Sulfatide (from GalCer)", "m"), "flip": (1, 4.4, "FAPP2 / ABCB1 flip\ncytosolic ⇄ lumen", "x"),
        "ceg": (2, 0, "Ceramide + glucose\n(lysosome, GBA1+SapC)", "m"), "ceg2": (2, 1.7, "Ceramide + galactose\n(GALC + SapA)", "m"),
        "gbp": (2, 3.4, "GlcSph (Gaucher)\nvia GBA2 / ASAH1", "m"),
        "r1": (3, 0, "↺ Glucose → glycolysis", "r"), "r2": (3, 1.7, "↺ Cer → ASAH1 →\nSph + FA → salvage", "r"),
        "r3": (3, 3.4, "↺ Galactose → Leloir\n→ glycolysis", "r"),
    },
    edges=[("cer", "hex", "UGCG (Glc) – Golgi\nUGT8 (Gal) – ER"), ("udp", "hex", ""), ("cerg", "hex", ""),
           ("hex", "lac", "B4GALT5/6"), ("hex", "sulf", "CST (UGT8 → CST)"), ("hex", "flip", ""),
           ("hex", "ceg", "GBA1 / GBA2"), ("hex", "ceg2", "GALC"), ("hex", "gbp", "", 0.0),
           ("ceg", "r1", ""), ("ceg2", "r2", ""), ("ceg2", "r3", "", 0.0)],
)

P["PA"] = dict(
    title="PA – Phosphatidic Acid", sub="Central phospholipid precursor and signalling lipid (mTOR, Raf) · ER / membranes",
    nodes={
        "g3p": (0, 0, "Glycerol-3-P + acyl-CoA", "m"), "lpa": (0, 1.2, "LPA", "m"), "dg": (0, 2.4, "DG", "m"),
        "pc": (0, 3.6, "PC / PE", "m"), "dhap": (0, 4.8, "DHAP (peroxisome)", "x"),
        "pa": (1, 2, "PA", "s"), "cdp": (1, 0.3, "CDP-DAG → PI, PG, CL", "m"), "sig": (1, 3.7, "Signalling: mTOR, Raf, Opi1", "m"),
        "dg2": (2, 0.3, "DG", "m"), "lpa2": (2, 2, "LPA + FA", "m"), "mg": (2, 3.7, "MG + Pi", "m"),
        "r1": (3, 0.3, "↺ DG → TG, PC, PE", "r"), "r2": (3, 2, "↺ LPA → AGPAT → PA\n(receptors LPAR1–6)", "r"),
        "r3": (3, 3.7, "↺ MG → glycerolipids;\nPi released", "r"),
    },
    edges=[("g3p", "lpa", "GPAT1-4"), ("lpa", "pa", "AGPAT1-5"), ("dg", "pa", "DGKα–ι"), ("pc", "pa", "PLD1/2"),
           ("dhap", "pa", "GNPAT → AGPS …", 0.0),
           ("pa", "cdp", "CDS1/2\n(CTP)"), ("pa", "sig", ""),
           ("pa", "dg2", "LIPIN1-3 (PAP)\nLPP1-3"), ("pa", "lpa2", "PLA1/PLA2\n(DDHD1, PLA2G4)"), ("lpa2", "mg", "LPP1-3", 0.0),
           ("dg2", "r1", ""), ("lpa2", "r2", ""), ("mg", "r3", "")],
)

P["PI"] = dict(
    title="PI – Phosphatidylinositol", sub="ER synthesis · phosphoinositide signalling (PIPs) · PI cycle",
    nodes={
        "g6p": (0, 0, "Glucose-6-P → myo-inositol\n(ISYNA1, IMPA1)", "m"), "pa": (0, 1.5, "PA", "m"), "cdp": (0, 2.7, "CDP-DAG", "m"),
        "pi": (1, 1.5, "PI", "s"), "pi4": (1, 0, "PI4P → PI(4,5)P₂ → PI(3,4,5)P₃", "m"),
        "pi35": (1, 3, "PI3P / PI(3,5)P₂\n(endosome / lysosome)", "m"), "lpi": (1, 4.3, "LPI (GPR55 ligand)", "m"),
        "dg": (2, 0, "DG + IP₃", "m"), "ip": (2, 1.5, "IP₃ → IP₂ → IP₁ → inositol", "m"),
        "lpi2": (2, 3, "GPI + FA", "m"), "pt": (2, 4.3, "Phosphoinositide\nphosphatases\n(PTEN, SHIP, MTM)", "m"),
        "r1": (3, 0, "↺ DG → DGK → PA\n→ CDS → CDP-DAG", "r"), "r2": (3, 1.5, "↺ myo-inositol →\nPIS (CDIPT) → PI", "r"),
        "r3": (3, 3, "↺ Glycerol / G3P;\nFA → acyl-CoA", "r"), "r4": (3, 4.3, "↺ Dephosphorylated →\nPI (PI cycle)", "r"),
    },
    edges=[("g6p", "pi", "myo-inositol", 0.0), ("pa", "cdp", "CDS1/2"), ("cdp", "pi", "PIS / CDIPT\n(+ inositol)"),
           ("pi", "pi4", "PI4KA/B, PIP5K"), ("pi", "pi35", "VPS34, PIKfyve"), ("pi", "lpi", "PLA1/PLA2 (DDHD1)"),
           ("pi4", "dg", "PLCβ/γ/δ"), ("dg", "ip", "", 0.0), ("pi35", "pt", "", 0.0), ("lpi", "lpi2", "lysophospholipase"),
           ("dg", "r1", ""), ("ip", "r2", "IMPA1 / IPPase"), ("lpi2", "r3", ""), ("pt", "r4", "")],
)

P["GM1"] = dict(
    title="GM1 – Monosialotetrahexosylganglioside", sub="Golgi synthesis → plasma membrane rafts → endolysosomal degradation (GLB1, HEXA/B) → salvage",
    nodes={
        "cer": (0, 0, "Ceramide → GlcCer\n(UGCG)", "m"), "lac": (0, 1.3, "LacCer", "m"), "gm3": (0, 2.6, "GM3", "m"), "gm2": (0, 3.9, "GM2", "m"),
        "gm1": (1, 3.2, "GM1", "s"), "gd1a": (1, 1.7, "GD1a / GD1b / GT1b", "m"), "pm": (1, 0.2, "PM raft / neuronal membranes\n(cholera toxin receptor)", "m"),
        "gm2d": (2, 0.2, "GM2 + galactose", "m"), "gm3d": (2, 1.7, "GM3 + GalNAc", "m"), "lacd": (2, 3.0, "LacCer + sialic acid", "m"),
        "glc": (2, 4.3, "GlcCer + gal →\nCer + glc →\nSph + FA", "m"),
        "r1": (3, 0.2, "↺ Galactose → Leloir pathway", "r"), "r2": (3, 1.7, "↺ GalNAc → hexosamine pathway", "r"),
        "r3": (3, 3.0, "↺ Sialic acid → CMP-sia\n(SLC17A5/sialin)", "r"), "r4": (3, 4.3, "↺ Sph → CerS (salvage)\nSph → S1P → lyase", "r"),
    },
    edges=[("cer", "lac", "B4GALT5/6"), ("lac", "gm3", "ST3GAL5"), ("gm3", "gm2", "B4GALNT1"), ("gm2", "gm1", "B3GALT4\n(+UDP-Gal)"),
           ("gm1", "gd1a", "ST3GAL2 / ST3GAL1\n(sialylation)", 0.0), ("gm1", "pm", "", 0.0),
           ("gm1", "gm2d", "GLB1 (β-gal) + SapB\n+ GM2A", 0.0), ("gm2d", "gm3d", "HEXA/HEXB\n+ GM2A"), ("gm3d", "lacd", "NEU1/3/4\n(sialidase)"),
           ("lacd", "glc", "GLB1 / GALC, GBA1"),
           ("gm2d", "r1", ""), ("gm3d", "r2", ""), ("lacd", "r3", ""), ("glc", "r4", "")],
)

P["SPH"] = dict(
    title="SPH – Sphingosine", sub="Sphingolipid backbone · salvage-pathway hub · S1P signalling (S1PR1–5) · irreversible exit via S1P lyase",
    nodes={
        "cer": (0, 0.4, "Ceramide", "m"), "s1p": (0, 1.9, "S1P (phosphatase route)", "m"), "gc": (0, 3.4, "Lysosomal complex\nsphingolipids", "x"),
        "sph": (1, 1.6, "Sphingosine", "s"), "s1pa": (1, 0, "S1P\n(secreted: SPNS2, MFSD2B)", "m"), "cer2": (1, 3.2, "Ceramide (salvage)", "m"),
        "sig": (1, 4.4, "N-methylsphingosine,\nSph-1-… derivatives", "x"),
        "ly": (2, 0, "Hexadecenal + P-Etn", "m"), "pal": (2, 1.6, "Hexadecenoic acid\n→ palmitoyl-CoA", "m"),
        "pe": (2, 3.2, "CDP-ethanolamine / PE", "m"),
        "r1": (3, 0, "↺ P-Etn → PE (Kennedy)", "r"), "r2": (3, 1.6, "↺ Palmitoyl-CoA → SPT\n(new sphingolipid)", "r"),
        "r3": (3, 3.2, "↺ Sph → CerS → SM / GlcCer\n(≥50 % of Sph salvaged)", "r"),
    },
    edges=[("cer", "sph", "Ceramidases\nASAH1/ASAH2/ACER"), ("s1p", "sph", "SGPP1/2"), ("gc", "sph", "lysosome\nGBA1/ASAH1"),
           ("sph", "s1pa", "SphK1/2 (+ATP)"), ("sph", "cer2", "CerS1-6\n(salvage)"), ("sph", "sig", "", 0.0),
           ("s1pa", "ly", "SGPL1\nS1P lyase"), ("ly", "pal", "ALDH3A2"), ("ly", "pe", "", 0.0),
           ("ly", "r1", "", 0.0), ("pal", "r2", ""), ("cer2", "r3", "")],
)

P["LPE"] = dict(
    title="LPE – Lysophosphatidylethanolamine", sub="Lyso-phospholipid intermediate of the Lands cycle · plasma / signalling pool",
    nodes={
        "pe": (0, 0.6, "PE", "m"), "pls": (0, 2.2, "Plasmalogen PE\n(ether-PE)", "m"), "plasma": (0, 3.8, "Plasma / dietary PE", "x"),
        "lpe": (1, 1.6, "LPE", "s"), "sig": (1, 0.1, "Signalling: GPR119, Ca²⁺ responses\n(pharmacological ligand)", "m"),
        "pe2": (1, 3.3, "PE (re-acylated)", "m"),
        "gpe": (2, 1.6, "GPE + FA", "m"), "eg": (2, 3.3, "Ethanolamine + G3P", "m"),
        "r1": (3, 1.6, "↺ FA → acyl-CoA\n→ LPEAT / β-oxidation", "r"), "r2": (3, 3.3, "↺ Ethanolamine → ETNK\n→ CDP-Etn → PE", "r"),
        "r3": (3, 0.1, "↺ G3P → glycolysis / GPAT", "r"),
    },
    edges=[("pe", "lpe", "PLA2 (cPLA2, sPLA2)"), ("pls", "lpe", "plasmalogen PLA2\n(iPLA2)", 0.0), ("plasma", "lpe", "PLA1 / LCAT-like", 0.0),
           ("lpe", "sig", ""), ("lpe", "pe2", "LPEAT\nMBOAT1/2 + acyl-CoA"), ("lpe", "gpe", "Lysophospholipase\n(LYPLA1/2, PNPLA6)"),
           ("gpe", "eg", "GDPD1/5"), ("gpe", "r1", ""), ("eg", "r2", ""), ("eg", "r3", "", 0.0)],
)

P["GM2"] = dict(
    title="GM2 – Ganglioside GM2", sub="Golgi-synthesised ganglioside · lysosomal degradation by β-hexosaminidase A + GM2 activator (Tay-Sachs / Sandhoff / AB variant)",
    nodes={
        "cer": (0, 0, "Ceramide → GlcCer", "m"), "lac": (0, 1.3, "LacCer", "m"), "gm3": (0, 2.6, "GM3", "m"),
        "gm2": (1, 1.5, "GM2", "s"), "gm1": (1, 0, "GM1 / GD2 / GD1a", "m"), "pm": (1, 3.2, "Plasma membrane →\nendocytosis → late endosome", "m"),
        "gm3d": (2, 1.5, "GM3 + GalNAc", "m"), "lacd": (2, 3.2, "LacCer + sialic acid", "m"),
        "cerd": (2, 4.6, "GlcCer → Cer + Glc\n→ Sph + FA", "m"),
        "r1": (3, 1.5, "↺ GalNAc → hexosamine\nsalvage", "r"), "r2": (3, 3.2, "↺ Sialic acid → CMP-sialic\nacid → new gangliosides", "r"),
        "r3": (3, 4.6, "↺ Sph → CerS (salvage)\nGlc → glycolysis", "r"),
    },
    edges=[("cer", "lac", "B4GALT5/6"), ("lac", "gm3", "ST3GAL5"), ("gm3", "gm2", "B4GALNT1\n(+UDP-GalNAc)"),
           ("gm2", "gm1", "B3GALT4 (GM1)\nST8SIA1 (GD2)", 0.0), ("gm2", "pm", "", 0.0),
           ("gm2", "gm3d", "HEXA/HEXB\n+ GM2A activator"), ("gm3d", "lacd", "NEU sialidase", 0.0), ("lacd", "cerd", "GLB1/GALC, GBA1"),
           ("gm3d", "r1", ""), ("lacd", "r2", ""), ("cerd", "r3", "")],
)

P["Co"] = dict(
    title="Co – Coenzyme Q (CoQ₁₀, ubiquinone)", sub="Mitochondrial electron carrier & lipid antioxidant · isoprenoid + aromatic ring synthesis (note: “Co” read as coenzyme Q)",
    nodes={
        "ac": (0, 0, "Acetyl-CoA →\nmevalonate (HMGCR)", "m"), "fpp": (0, 1.3, "FPP → decaprenyl-PP\n(PDSS1/2)", "m"), "hb": (0, 2.6, "Tyrosine / Phe → 4-HB", "m"),
        "pre": (0, 3.9, "Polyprenyl-4-HB\n(COQ2 prenylation)", "m"),
        "coq": (1, 2.2, "CoQ₁₀ (ubiquinone)", "s"), "cplx": (1, 0.6, "Complex Q: COQ3–7, COQ9\n(inner mito. membrane)", "m"),
        "qh2": (1, 3.7, "CoQH₂ (ubiquinol)", "s"),
        "ox": (2, 0.4, "Redox use: Complexes I/II/III,\nETF-QO, DHODH, SQOR", "m"), "fsp": (2, 2.0, "Plasma membrane: FSP1\nrecycles CoQH₂ (anti-ferroptosis)", "m"),
        "cat": (2, 3.6, "Side-chain ω/β-oxidation →\nquinone-acid metabolites", "m"), "gl": (2, 5, "Glucuronide / sulfate\nconjugates", "m"),
        "r1": (3, 0.4, "↺ CoQH₂ ⇄ CoQ cycling\n(Q cycle)", "r"), "r2": (3, 2.0, "↺ CoQH₂ regenerates\nvitamin E / lipids", "r"),
        "r3": (3, 3.6, "↺ Isoprenoid / ring\nunits to pool", "r"), "r4": (3, 5, "↺ Biliary excretion", "r"),
    },
    edges=[("ac", "fpp", "mevalonate / FPS"), ("hb", "pre", "COQ2"), ("fpp", "pre", "", 0.0),
           ("pre", "cplx", "COQ3–7, COQ9", 0.0), ("cplx", "coq", "COQ7 final"), ("coq", "qh2", "Complex I/II →\nreduction (+2e⁻)", 0.0),
           ("coq", "ox", ""), ("qh2", "fsp", "FSP1 / NQO1"), ("qh2", "cat", "ω-oxidation"), ("cat", "gl", "UGT / SULT"),
           ("ox", "r1", ""), ("fsp", "r2", ""), ("cat", "r3", "", 0.0), ("gl", "r4", "")],
)

P["ST"] = dict(
    title="ST – Sterol (Cholesterol)", sub="Membrane sterol · de novo ER synthesis + LDL uptake · SREBP2 regulation · bile acid / hormone fates (note: “ST” read as free sterol)",
    nodes={
        "ac": (0, 0, "Acetyl-CoA → HMG-CoA\n→ mevalonate (HMGCR)", "m"), "fpp": (0, 1.4, "IPP → FPP → squalene\n(SQLE) → lanosterol", "m"),
        "cho": (0, 2.8, "Lanosterol → cholesterol\n(DHCR7 / DHCR24)", "m"), "ldl": (0, 4.2, "LDL uptake\n(LDLR → LIPA, NPC1/2)", "x"),
        "st": (1, 1.9, "Cholesterol", "s"), "ce": (1, 0.3, "Cholesteryl ester\n(ACAT1/2)", "m"), "ox": (1, 3.4, "Oxysterols\n(CH25H, CYP27A1, CYP46A1)", "m"),
        "mem": (1, 5, "Membrane / rafts; SREBP2\nfeedback (INSIG, SCAP)", "m"),
        "ba": (2, 0.3, "Bile acids\n(CYP7A1, CYP27A1)", "m"), "hor": (2, 1.9, "Pregnenolone → steroid\nhormones (CYP11A1)", "m"),
        "vd": (2, 3.4, "7-DHC → vitamin D", "m"), "bile": (2, 5, "Biliary cholesterol\n(ABCG5/8)", "m"),
        "r1": (3, 0.3, "↺ Enterohepatic recycling\n(ASBT, NTCP; ~95 %)", "r"), "r2": (3, 1.9, "↺ Reverse cholesterol transport\nABCA1/G1 → HDL → SR-B1", "r"),
        "r3": (3, 3.4, "↺ ChE hydrolysis\n(LIPA, NCEH1)", "r"), "r4": (3, 5, "↺ Intestinal reabsorption\nNPC1L1 / faecal loss", "r"),
    },
    edges=[("ac", "fpp", "HMGCS1, FDPS"), ("fpp", "cho", "LSS, CYP51 …"), ("cho", "st", ""), ("ldl", "st", "", 0.0),
           ("st", "ce", "ACAT / SOAT"), ("st", "ox", "CYP / ROS"), ("st", "mem", "", 0.0),
           ("st", "ba", "CYP7A1 (classic)"), ("st", "hor", "CYP11A1"), ("st", "vd", "(via 7-DHC)", 0.0), ("st", "bile", "", 0.0),
           ("ba", "r1", ""), ("hor", "r2", "", 0.0), ("ce", "r3", "", 0.0), ("bile", "r4", "")],
)

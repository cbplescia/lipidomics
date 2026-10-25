#!/usr/bin/env python3
"""Curated human gene lists linked to the lipids changed by IFNB/IFNG.

Literature-based, NOT derived from the lipidomics data. Mouse orthologs: same symbol
in title case for most (e.g. Pdcd6ip, Smpd1). Writes grouped and flat lists to this folder.
Run:  python gene_lists/make_gene_lists.py
"""
import os

R = lambda p, a, b: [f"{p}{i}" for i in range(a, b + 1)]

# ------------------------------------------------------------------ 1. effectors / readers / processes
EFFECTORS = {
    "BMP/ILV biogenesis, ESCRT": ["PDCD6IP", "HGS", "STAM", "STAM2", "TSG101", "VPS28", "VPS37A", "MVB12A", "SNF8", "VPS25", "VPS36",
        "CHMP1A", "CHMP1B", "CHMP2A", "CHMP2B", "CHMP3", "CHMP4A", "CHMP4B", "CHMP4C", "CHMP5", "CHMP6", "IST1", "VPS4A", "VPS4B",
        "SDCBP", "VTA1", "CD63", "CD9", "CD81", "TSPAN6"],
    "BMP / lysosomal lipid hydrolase cofactors": ["SMPD1", "ASAH1", "GBA1", "PSAP", "GM2A", "NPC1", "NPC2", "HSPA1A", "HSPA1B", "HSPA8",
        "CTSD", "CTSB", "CTSL", "HEXA", "HEXB", "GLB1", "LIPA", "SCARB2", "LAMP1", "LAMP2"],
    "Ceramide effectors": ["SMPD3", "PPP2CA", "PPP2CB", "PPP1CA", "SET", "PRKCZ", "KSR1", "FAS", "CD40", "MAP1LC3B", "CERT1", "CPTP", "TXNIP", "NLRP3"],
    "Cholesterol sensing / transport at endolysosomes": ["SLC38A9", "GPR155", "OSBP", "OSBPL1A", "VAPA", "VAPB", "STARD3", "STARD3NL", "GRAMD1A", "GRAMD1B",
        "GRAMD1C", "SOAT1", "SOAT2", "SREBF2", "SCAP", "INSIG1", "INSIG2", "LDLR", "ABCA1", "ABCG1", "SCARB1", "NR1H2", "NR1H3",
        "CH25H", "CYP7B1", "RAB7A", "RILP", "ORP1L"],
    "mTORC1 / lysosome signalling": ["MTOR", "RPTOR", "MLST8", "RHEB", "TSC1", "TSC2", "FLCN", "FNIP1", "FNIP2", "RRAGA", "RRAGB", "RRAGC", "RRAGD",
        "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5", "TFEB", "TFE3", "MITF", "MCOLN1", "TPCN1", "TPCN2", "PPP3CA", "CAMKK2"],
    "PA / LPA / DG effectors": ["RAF1", "PLD1", "PLD2", "LPAR1", "LPAR2", "LPAR3", "LPAR4", "LPAR5", "LPAR6", "ENPP2", "PRKCA", "PRKCB", "PRKCD", "PRKCE",
        "PRKCQ", "PRKD1", "PRKD2", "PRKD3", "UNC13A", "UNC13B", "RASGRP1", "RASGRP3", "CHN1", "CHN2", "DAGLA", "DAGLB", "MGLL"],
    "Phosphoinositide effectors (PI3P, PI(3,5)P2, PI4P)": ["PIK3C3", "PIK3R4", "BECN1", "UVRAG", "ATG14", "EEA1", "ZFYVE1", "ZFYVE16", "WDFY2", "SNX1", "SNX2", "SNX3",
        "SNX5", "SNX6", "SNX27", "PIKFYVE", "VAC14", "FIG4", "PI4K2A", "PI4K2B", "PI4KB", "PI4KA", "PIP5K1C", "TLR4", "GPR55", "GPR132", "GPR4"],
    "Autophagy / LC3 lipidation (PE) and IFN-induced GTPases": ["ATG3", "ATG5", "ATG7", "ATG10", "ATG12", "ATG16L1", "MAP1LC3A", "MAP1LC3B", "MAP1LC3C", "GABARAP",
        "GABARAPL1", "GABARAPL2", "RUBCN", "CYBB", "NCF1", "NCF2", "NCF4", "ATP6V1A", "ATP6V1B2", "ATP6V0A1", "ATP6V0D1", "ATP6AP1", "TECPR1", "IRGM",
        "GBP1", "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "GBP7", "LRRK2"],
    "Endosome maturation / fusion": ["RAB5A", "RAB5B", "RAB5C", "RAB4A", "RAB7A", "RAB9A", "RAB11A", "ARF6", "ARL8B", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39", "VPS41",
        "STX7", "STX8", "VTI1B", "VAMP7", "VAMP8", "SYT7"],
    "Gangliosides / rafts / sphingosine": ["CAV1", "FLOT1", "FLOT2", "NTRK1", "SPHK1", "SPHK2", "S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5", "SGPL1"],
    "Interferon-induced membrane / antiviral": ["IFITM1", "IFITM2", "IFITM3", "CH25H", "OAS1", "MX1", "ISG15", "STAT1", "STAT2", "IRF1", "IRF9"],
}

# ------------------------------------------------------------------ 2. lipid-modifying genes (synthesis, remodeling, degradation, transport/flip/flop)
MODIFIERS = {
    "Flippases, P4-ATPases (to cytosolic leaflet) + CDC50 cofactors": ["ATP8A1", "ATP8A2", "ATP8B1", "ATP8B2", "ATP8B3", "ATP8B4", "ATP9A", "ATP9B",
        "ATP10A", "ATP10B", "ATP10D", "ATP11A", "ATP11B", "ATP11C", "TMEM30A", "TMEM30B", "TMEM30C"],
    "Scramblases": ["XKR4", "XKR8", "XKR9", "BSG", "NPTN", "ANO3", "ANO4", "ANO5", "ANO6", "ANO7", "ANO9", "ANO10", "PLSCR1", "PLSCR2", "PLSCR3", "PLSCR4",
        "TMEM41B", "VMP1", "ATG9A", "ATG9B", "TMEM63A", "TMEM63B"],
    "Floppases / ABC lipid exporters and sterol efflux": ["ABCA1", "ABCA3", "ABCA4", "ABCA7", "ABCA12", "ABCB1", "ABCB4", "ABCC1", "ABCG1", "ABCG2", "ABCG4",
        "ABCG5", "ABCG8", "ABCD1", "ABCD2", "ABCD3", "SCARB1"],
    "Lipid transfer / contact-site proteins": ["OSBP", "OSBPL1A", "OSBPL2", "OSBPL3", "OSBPL5", "OSBPL8", "OSBPL9", "OSBPL10", "OSBPL11", "STARD2", "STARD3", "STARD10", "STARD4", "STARD5",
        "CERT1", "CPTP", "GLTP", "PITPNA", "PITPNB", "PITPNM1", "PITPNM2", "PITPNM3", "VPS13A", "VPS13B", "VPS13C", "VPS13D", "ATG2A", "ATG2B", "TMEM24",
        "ESYT1", "ESYT2", "ESYT3", "GRAMD1A", "GRAMD1B", "GRAMD1C", "NPC1", "NPC2", "SCARB2", "LAMP2", "SPNS1", "SPNS2", "MFSD2A", "MFSD2B", "CD36", "FABP4", "FABP5"],
    "PC": ["PCYT1A", "PCYT1B", "CHKA", "CHKB", "CHPT1", "CEPT1", "PEMT", "LPCAT1", "LPCAT2", "LPCAT3", "LPCAT4", "PLA2G4A", "PLA2G2A", "PLA2G5", "PLA2G10",
        "PLA2G6", "PLA2G7", "PLA2G15", "PLD1", "PLD2", "SGMS1", "SGMS2", "PNPLA6", "PNPLA7", "GPCPD1", "LCAT", "PCYT2"],
    "PE": ["PCYT2", "EPT1", "SELENOI", "CEPT1", "PISD", "PTDSS2", "PEMT", "LPCAT3", "MBOAT1", "MBOAT2", "PLA2G4A", "PLA2G6", "ATG3", "ATG7"],
    "PS": ["PTDSS1", "PTDSS2", "PISD", "ATP8A1", "ATP8A2", "ATP11A", "ATP11B", "ATP11C", "TMEM30A", "XKR8", "ANO6", "PLSCR1", "OSBPL5", "OSBPL8", "PLA1A", "PLA2G4A", "LPCAT3"],
    "PI and phosphoinositides": ["CDIPT", "CDS1", "CDS2", "PIK3C2A", "PIK3C2B", "PIK3C2G", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIKFYVE", "PIP5K1A", "PIP5K1B",
        "PIP5K1C", "PIP4K2A", "PIP4K2B", "PIP4K2C", "PI4KA", "PI4KB", "PI4K2A", "PI4K2B", "PTEN", "MTM1", "MTMR1", "MTMR2", "MTMR3", "MTMR4", "MTMR6",
        "INPP4A", "INPP4B", "INPP5D", "INPP5E", "INPP5F", "INPP5B", "INPPL1", "OCRL", "SACM1L", "FIG4", "SYNJ1", "SYNJ2", "PLCB1", "PLCB2", "PLCB3", "PLCG1", "PLCG2",
        "PLCD1", "PLCE1", "PITPNA", "PITPNB", "MBOAT7", "AGPAT3", "LCLAT1", "DGKE"],
    "LPI": ["MBOAT7", "PLA2G4A", "PLA2G6", "DDHD1", "DDHD2", "ABHD12", "GPR55"],
    "PA, LPA, cPA": ["PLD1", "PLD2", "PLD3", "PLD4", "DGKA", "DGKB", "DGKD", "DGKE", "DGKG", "DGKH", "DGKI", "DGKK", "DGKQ", "DGKZ", "AGPAT1", "AGPAT2", "AGPAT3",
        "AGPAT4", "AGPAT5", "GPAT1", "GPAT2", "GPAT3", "GPAT4", "PLPP1", "PLPP2", "PLPP3", "LPIN1", "LPIN2", "LPIN3", "ENPP2", "DDHD1", "PLA1A", "PLA2G4A"],
    "DG / TG / lipid droplets": ["DGAT1", "DGAT2", "MOGAT1", "MOGAT2", "PNPLA2", "ABHD5", "LIPE", "MGLL", "LIPA", "G0S2", "PLIN1", "PLIN2", "PLIN3", "PLIN4", "PLIN5",
        "BSCL2", "FITM1", "FITM2", "CIDEA", "CIDEB", "CIDEC", "DAGLA", "DAGLB", "ABHD6", "ABHD12", "DGKA", "DGKZ", "LPIN1", "GPAT4", "AGPAT3", "LPL", "ATGL"],
    "BMP / LBPA": ["CLN5", "PLD3", "PLD4", "PLA2G15", "ABHD6", "CLN3", "GRN", "TMEM106B", "LRRK2", "NPC1", "NPC2", "PSAP", "SMPD1", "ASAH1", "GBA1", "CTSA"],
    "Cardiolipin / PG": ["PGS1", "PTPMT1", "CRLS1", "TAMM41", "CDS1", "CDS2", "TAFAZZIN", "TAZ", "LCLAT1", "PLA2G6", "PNPLA8", "AGK", "PRELID1", "PRELID3B", "TRIAP1", "PLSCR3", "CHCHD10"],
    "Ceramide / sphingolipid synthesis and turnover": ["SPTLC1", "SPTLC2", "SPTLC3", "SPTSSA", "SPTSSB", "ORMDL1", "ORMDL2", "ORMDL3", "KDSR", "CERS1", "CERS2", "CERS3", "CERS4", "CERS5", "CERS6",
        "DEGS1", "DEGS2", "SMPD1", "SMPD2", "SMPD3", "SMPD4", "ENPP7", "ASAH1", "ASAH2", "ACER1", "ACER2", "ACER3", "SGMS1", "SGMS2", "SAMD8", "CERK", "CERT1", "CPTP",
        "UGCG", "GBA1", "GBA2", "UGT8", "SPHK1", "SPHK2", "SGPL1", "SGPP1", "SGPP2", "SPNS2", "FA2H", "ELOVL1", "ELOVL3", "ELOVL4", "ELOVL7", "ARV1"],
    "Glycosphingolipids / gangliosides / sulfatide": ["UGCG", "B4GALT5", "B4GALT6", "ST3GAL5", "B4GALNT1", "B3GALT4", "ST3GAL2", "ST8SIA1", "GLB1", "HEXA", "HEXB", "GM2A", "NEU1", "NEU2",
        "NEU3", "NEU4", "GBA1", "GBA2", "GALC", "ARSA", "GAL3ST1", "UGT8", "PSAP", "GLA"],
    "Cholesterol synthesis, esterification, trafficking": ["HMGCS1", "HMGCR", "MVK", "PMVK", "MVD", "FDPS", "FDFT1", "SQLE", "LSS", "CYP51A1", "TM7SF2", "LBR", "MSMO1", "NSDHL", "HSD17B7", "EBP", "SC5D",
        "DHCR7", "DHCR24", "SREBF1", "SREBF2", "SCAP", "INSIG1", "INSIG2", "MBTPS1", "MBTPS2", "NPC1", "NPC2", "LIPA", "SOAT1", "SOAT2", "NCEH1", "LIPE", "CES1", "LCAT",
        "CH25H", "CYP7B1", "CYP27A1", "CYP46A1", "CYP11A1", "NR1H2", "NR1H3", "ABCA1", "ABCG1", "SCARB1", "LDLR", "OSBP", "OSBPL1A", "STARD3", "STARD4", "GRAMD1A", "GRAMD1B", "GRAMD1C"],
    "Ether lipids / plasmalogens": ["GNPAT", "AGPS", "FAR1", "FAR2", "PEDS1", "TMEM189", "PEX3", "PEX5", "PEX7", "PEX19", "PEX1", "PEX2", "PEX12", "PEX14", "SELENOI", "EPT1", "PLA2G6", "PLA2G4A", "ACSL1", "ACSL4"],
    "Acylcarnitines / FA oxidation": ["CPT1A", "CPT1B", "CPT1C", "CPT2", "CRAT", "CROT", "SLC25A20", "SLC22A5", "ACSL1", "ACSL3", "ACSL4", "ACSL5", "ACSL6", "ACADVL", "ACADM", "HADHA", "HADHB", "ECI1", "PPARA", "PPARG", "PPARGC1A"],
    "CoQ": ["PDSS1", "PDSS2", "COQ2", "COQ3", "COQ4", "COQ5", "COQ6", "COQ7", "COQ8A", "COQ8B", "COQ9", "COQ10A", "COQ10B", "PPTC7", "FDX2", "FDXR"],
    "PUFA incorporation / elongation / desaturation": ["FADS1", "FADS2", "FADS3", "ELOVL2", "ELOVL5", "ELOVL6", "ELOVL4", "SCD", "SCD5", "ACSL1", "ACSL3", "ACSL4", "ACSL5", "ACSL6", "LPCAT3", "MBOAT7", "MBOAT1",
        "AGPAT3", "MFSD2A", "PLA2G4A", "PLA2G6", "PLA2G7", "PLA2G15", "GPX4", "AIFM2", "ALOX5", "ALOX12", "ALOX15", "PTGS1", "PTGS2"],
    "Lysophospholipids": ["PLA2G4A", "PLA2G2A", "PLA2G5", "PLA2G10", "PLA2G6", "PLA2G7", "PLA2G15", "LCAT", "LPCAT1", "LPCAT2", "LPCAT3", "LPCAT4", "MBOAT1", "MBOAT2", "MBOAT7", "LPGAT1", "AGK", "ENPP2", "GDPD1", "GDPD3", "MFSD2A", "SPNS1", "PNPLA6"],
}


# aliases / off-target entries removed so every symbol is a current HGNC symbol
DROP = {"ORP1L", "ATGL", "TAZ", "EPT1", "TMEM189", "FDX2", "FDXR", "CHCHD10", "TXNIP", "NLRP3", "GLA", "CYP11A1"}


def uniq(xs):
    return sorted(set(xs) - DROP)


def write(name, groups):
    os.makedirs(os.path.dirname(__file__) or ".", exist_ok=True)
    d = os.path.dirname(os.path.abspath(__file__))
    with open(f"{d}/{name}_grouped.txt", "w") as f:
        for k, v in groups.items():
            f.write(f"# {k}\n{', '.join(uniq(v))}\n\n")
    flat = uniq(g for v in groups.values() for g in v)
    open(f"{d}/{name}_flat.txt", "w").write("\n".join(flat) + "\n")
    return flat


if __name__ == "__main__":
    a = write("1_effector_proteins", EFFECTORS)
    b = write("2_lipid_modifying_genes", MODIFIERS)
    print(len(a), "effector genes;", len(b), "lipid-modifying genes;", len(set(a) | set(b)), "combined")

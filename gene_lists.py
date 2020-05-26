gwas_hits = {
 'rs35749011': ['GBA'],
 'rs823118': ['NUCKS1', 'SLC41A1'],
 'rs10797576': ['SIPA1L2'],
 'rs6430538': ['TMEM163', 'CCNT2'],
 'rs1474055': ['STK39'],
 'rs115185635': ['CHMP2B'],
 'rs12637471': ['MCCC1'],
 'rs34311866': ['TMEM175', 'DGKQ'],
 'rs11724635': ['FAM200B'],
 'rs6812193g': ['FAM47E'],
 'rs356182': ['SNCA'],
 'rs199347': ['KLHL7', 'NUPL2', 'GPNMB'],
 'rs591323': ['MICU3'],
 'rs117896735': ['BAG3'],
 'rs3793947': ['DLG2'],
 'rs76904798h': ['LRRK2'],
 'rs11060180': ['OGFOD2'],
 'rs11158026': ['GCH1'],
 'rs1555399': ['TMEM229B'],
 'rs2414739': ['VPS13C'],
 'rs14235': ['ZNF646', 'KAT8'],
 'rs17649553': ['ARHGAP27', 'CRHR1', 'MAPT', 'KANSL1'],
 'rs12456492': ['SYT4'],
 'rs62120679': ['LSM7'],
 'rs8118008': ['DDRGK1'],
 'familial': ['PARK7', 'PINK1'],
 'rs4653767': ['ITPKB'],
 'rs353116': ['SCN3A'],
 'rs4073221': ['SATB1'],
 'rs12497850': ['NCKIPSD'],
 'rs143918452': ['ALAS1', 'DNAH1', 'BAP1', 'PHF7', 'NISCH', 'ITIH4'],
 'rs78738012': ['ANK2', 'CAMK2D'],
 'rs2694528': ['ELOVL7'],
 'rs9468199': ['ZNF184'],
 'rs2740594c': ['CTSB'],
 'rs2280104': ['SORBS3', 'PDLIM2', 'C8orf58', 'BIN3'],
 'rs13294100': ['SH3GL2'],
 'rs10906923': ['FAM171A1'],
 'rs8005172': ['GALC'],
 'rs11343': ['COQ7'],
 'rs4784227': ['TOX3'],
 'rs601999': ['ATP6V0A1', 'PSMC3IP', 'TUBG2']
}

import itertools

snps_gene_list = list(itertools.chain.from_iterable(gwas_hits.values()))




fig1_model_heatmap = ["TH", "NR4A2", "FOXA2", "LMX1A", "LMX1B", "DDC", "KCNJ6", "PBX1", "LMO3", "CALB1", "SLC18A2", "TPH1", "SLC18A1", 
                      "SNAP25", "MAP2", "SYT1", "DCX", "VIM", "HES1","NFIA","RFX4", "SLC1A3", "SLIT2", "OTX2" ]
sfig1_model_matrix = ["TH","FOXA2","KCNJ6", "LMX1A", "LMX1B", "SLC25A5", "IGF1", "DDC", "ALDH1A1", "PBX1", "NETO2", "TMCC3"]
sfig1_model_umaps_heatmap = [ "TH", "DCX", "VIM", "SNAP25", "TPH1", "PBX1", "KCNJ6", "SLC18A2", "LMO3", "CALB1", "MARCKS", "NEFM"]
sfig1_welch_heatmap = ["TH","PLXDC1", "P2RY12", "MOBP", "OLIG1", "AQP4", "MGP", "SYT1", "GAD1"]
sfig1_model_violinplot = ["SLC17A6", "GRIN2B", "GLUL", "GLS", "GABBR1", "GABBR2", "ACSL1", "CALCA", "HMGA2", "HMGB2", "HMGB3", "GFAP", 
                          "AQP4", "S100B", "PDGFRA"]

# sfig1_model_review = ["SOX2", "RFX4", "HMGA1", "HMGB2", "CORIN", "WNT1", "ZEB2", "TOX", "JADE-1", "LMX1A", "FOXA2", "NEUROG2", "ASCL1", 
#                      "TUBB3", "SNAP25", "NEUROD1", "NEUROG2", "NR4A2", "TH", "NKX6-2","SLC18A2", "SLC6A3", "SOX6", "LMO3", "ALDH1A1"]

sfig1_model_review = ["SOX2", "RFX4", "HMGA1", "HMGB2", "CORIN", "ZEB2", "TOX", "LMX1A", "FOXA2", "ASCL1", 
                      "TUBB3", "SNAP25", "NEUROD1", "NR4A2", "TH", "NKX6-2","SLC18A2", "SOX6", "LMO3", "ALDH1A1"]

sfig1_model_review_short = ["ALDH1A1", "SOX6", "OTX2"]
# sfig1_model_poulin_extra = ["ALDH1A1", "SOX6", "OTX2", "SLC32A1", "VIP"]

fig2_lfc = ["SNCA", "MAPT", "UCHL1", "ATP13A2", "BDNF", "STMN1", "STMN2", "STMN3", "STMN4", "TUBB2A", "TUBB4A", 
        "TUBB2B", "HMP19", "NSG1", "SYT4", "SYT5", "SYT17", "SYT1", "CALM2","CALM3", "DOC2B", "SNAP25", "RAB3C", 
        "VAMP2", "RAB3A", "SYN1", "LIN7B","PPFIA2", "RTN1", "RTN3", "RTN4", "MAP1B", "MAP1LC3A", "ACAT2", "FDFT1",
        "HMGCR", "HMGCS1", "IDI1", "MVD", "MSMO1", "SQLE", "INSIG1", "MT-CO1", "MT-CO2", "MT-ND3", "COX17", 
        "CALR", "CANX", "XBP1", "PDIA4", "PDIA3", "PDIA6", "ANXA5", "ANXA6","ANXA2", "LAMP1", "LAMP2", "CTSB", 
        "CTSD", "CTSA", "CTSH"]


fig2_dan1_violinplots = {
    'PD': ["SNCA", "MAPT", "ATP13A2"],
    'SynapticSignaling' : ["SYT4", "CALM2", "HMP19", "RAB3A"],
    "Cholesterol" : ["ACAT2", "SQLE", "HMGCS1", "HMGCR"],
    "Supplementary" : ["CALR", "PDIA6", "CANX", "ANXA6"]
}

fig2_dan2_violinplots = {
    "OxidativePhosphorylation": ["COX17", "MT-ND3", "MT-CO3", "MT-CO1"],
    "Supplementary" : ["TUBB2B", "TUBA1A", "MAPT", "HMGCS1"],
    "ChromatinRemodelling" : ["EID1","HIST2H2BE"],
    "HeatShock" : ["DNAJB6", "DNAJA1"],
    "Review-extra" :["CST3", "SEPP1", "TXNIP", "PON2"]
}


# With respect to what? Include cell type
fig3_dan1_violinplots = { 'Cholesterol Biosynthesis' : ['SQLE', 'HMGCR'],
              'Oxidative phosphorylation' : ["MT-CO1", "NDUFA5"],
              'PD' : ["SNCA"]
}

fig3_dan2_violinplots = {
    "ChromatinRemodelling" : ["EID1","HIST2H2BE"],
    "HeatShock" : ["DNAJB6", "DNAJA1"]
}

# Legend change
fig3_dan1dan2_lfc_down = ["SNCA", "BDNF", "ATP13A2", "SQLE", "ACAT2", "ACSL3", "HMGCR", "INSIG1", "MSMO1", "CYP46A1", "KIF1A", 
                     "KIF3A","KRAS", "TUBB2A", "TUBB4B", "TUBA1B", "TUBA1A", "MT-CO1", "MT-CO3", "MT-CO2", "MT-ND3", "MT-ND4", 
                     "NDUFB6", "NDUFA5"]


fig3_dan1dan2_lfc_up = ["HIST2H2BE", "EID1", "HIST1H2AC", "HIST3H2A", "FOXB1",  "DUSP6", "BCLAF1", "SERINC3","IP6K2", "HSPA8", 
                     "DNAJA1", "HSPB11", "BCLAF1", "BAD", "BAG1", "DNAJB6", "UCHL1", "UBL3", "UBB", "UBA52", "HSPA1A", 
                     "HSP90AB1", "HSP90AA1", "HSPA1B", "HSPA6", "FAU", "PLCG2"]

fig4_dan1dan2_lfc = ["HMGCS1", "ACAT2", "IDI1", "MSMO1", "HMGCR", "FDFT1", "INSIG1", "SQLE", "MVD", "TUBA1A", "TUBA1B", "TUBB2A", 
                     "TUBB2B", "TUBB", "CALM2","CALM3", "NDUFB2", "NDUFA3", "NDUFA4", "COX6B1", "COX6C", "COX7B", "TPI1", "ALDOA", 
                     "ALDOC", "PGAM1", "LDHA", "ENO1", "ENO2", "ATF3", "ATF4", "DDIT3", "DDIT4", "EIF1", "DNAJB9", "PDIA3", "PDIA4", 
                     "PDIA6", "XBP1", "CALR", "CANX", "EIF4EBP1", "BBC3", 'DNAJC3', "DNAJC10", "DNAJC12", "DNAJB11","HSPA5","HSP90B1"]

fig4_violinplots = {'up':['MANF', 'HSPA5'],
                    'down':["NDUFA4", "SQLE"]
}


fig5_violinplots = {'up': ["MAPT", "BASP1", "ID4"],
                    'down' : ["ELAVL4", "SNCA", "GPX4"]
}

fig5_lfc_down = ["ELAVL4", "SNCA", "SNCG", "TXNDC16", "TXNDC17", "GPX4", "GCH1","GCHFR", "FAAP20", "RNMT", "POLR2K", "GTF2H5", "TMEM258", 
            "DPM3", "KRTCAP2", "INSIG1", "HMGCR", "MSMO1", "HMGCS1", "SQLE", "ACAT2", "NDUFA3", "NDUFA1", "COX6A1", "MT-ND4L", "RPL37A", 
            "RPL38", "RPS27", "STMN1", "STMN4", "STMN2", "CALM1", "CALM2"]

fig5_lfc_up = ["ID1", "ID2", "ID3", "ID4", "H3F3B", "HIST1H4C", "H1FX", "H1F0", "HIST2H2BE", "HIST1H1E", "HIST1H2BD", "KDM5B", "KDM5A", 
               "SMARCD3", "NCOR1", "H2AFJ", "PRMT2", "ARID1A",  "KMT2E", "KMT2E", "KMT5B", "KMT2A", "JMJD1C", "BAZ1B", "TADA3", "HDAC3", 
               "ING2", "NLRP1", "BASP1","HRAS", "JUN", "BCL2L1", "FYN", "UBE2S", "UBE2R2", "SIAH2", "UBE2Q1","KIF1A", "KIF5C", "MAPT", 
               "MAP6", "TUBA1B", "SFPQ" ]


fig6_violinplots_rot = {
    "SynapticSignalling":["RAB1A", "CALM2"],
    "Glycolysis":["LDHA", "ENO1"]
}


fig6_violinplots_tun = {
    "StressResponse" :["UBC", "RNF187", "PARK7", "PDIA4"],
}

fig6_dan1dan2_tun_lfc_up = [ "HSPA8", "HSPA5", "PARK7", "DNAJC4", "DNAJC8", "DNAJB9", "DNAJA1", "DNAJC12", 'DNAJB6', "DNAJB14", "MEAF6", 
                             "MORF4L1", "MORF4L2","XBP1", "PDIA4", "CALR", "DDIT3", "PPP1R15A", "HERPUD1", "UBB", 'UBC', "UCHL1", "SUMO1", 
                             "RNF5", "RNF187", "SIAH2", "SKP1", 'PSMD14', 'PSMD8', 'PSMB3', 'PSMB6', 'PSMC4', 'PSMB1', 'PSMB2', 'PSMD13']


fig6_dan1dan2_tun_lfc_down = ["MAP2", "PCLO", "ITSN1", "PREPL", "KIAA1107", "KIAA1109", "RAB3B", "RAB3C", "RAB30", "RAB6B", "KIF2A", "KIF1B", 
                          "KIF1A", "KIF3A", "KIF5C", "KIF21A", "SNCA"]

fig6_dan1dan2_rot_lfc_up = [ "GAPDH", "LDHA", "ENO1", "TPI1", "ALDOA", "ENO2", "PGK1", "PKM",  "VAMP2", "RAB3A", "RAB3C", "RAB1A", "SNAP25", 
                             "CALM2", "CALM3", "SQLE", "HMGCS1", "MSMO1", "CYP46A1", "TXN","GLRX", "SEPW1", "GPX3"]

fig6_dan1dan2_rot_lfc_down = ["MEST", "ELAVL4", "RPS29", "RPS27",
            "RPL38", "RPL37A", "NDUFA3",
            "MT-ND3", "MT-ND4L", "SNCA", "SNCG", "NSG1", "TUBB", 'STMN4',
            "TUBB2B", "KIF5C", "DCX"]


fig7_umap_gwas = ["SNCA", "MAPT", "SYT4", "SCN3A", "PDLIM2"]


# DBH absent
poulin = ["ALDH1A1", "SOX6", "OTX2", "TH", "CCK","CALB1", "SLC17A6", "SLC18A2","DDC","PITX3"]
#poulin = ["ALDH1A1", "SOX6", "OTX2"]


figure_map = {
    'fig1-model-heatmap': fig1_model_heatmap,
    'sfig1-model-umap-heatmap': sfig1_model_umaps_heatmap,
    'sfig1-model-violinplot': sfig1_model_violinplot,
    'sfig1-model-matrixplot': sfig1_model_matrix,
    'sfig1-welch-heatmap' : sfig1_welch_heatmap,
    'sfig1-violinplot' : sfig1_model_violinplot,
    'sfig1-model-new-review' : sfig1_model_review,
    'sfig1-model-new-review-short' : sfig1_model_review_short,
    'fig2-lfc' : fig2_lfc,
    'fig2-dan1-violinplots': fig2_dan1_violinplots,
    'fig2-dan2-violinplots': fig2_dan2_violinplots,
    'fig3-dan1-violinplots' : fig3_dan1_violinplots,
    'fig3-dan2-violinplots' : fig3_dan2_violinplots,
    'fig3-dan1dan2-lfc-down' : fig3_dan1dan2_lfc_down,
    'fig3-dan1dan2-lfc-up' : fig3_dan1dan2_lfc_up,
    'fig4-violinplots' : fig4_violinplots,
    'fig4-dan1dan2-lfc' : fig4_dan1dan2_lfc, 
    'fig5-violinplots' : fig5_violinplots,
    'fig5-lfc-down' : fig5_lfc_down,
    'fig5-lfc-up' : fig5_lfc_up,
    'fig6-violinplots-rot' : fig6_violinplots_rot,
    'fig6-violinplots-tun' : fig6_violinplots_tun,
    'fig6-dan1dan2-tun-lfc-up' : fig6_dan1dan2_tun_lfc_up,
    'fig6-dan1dan2-tun-lfc-down' : fig6_dan1dan2_tun_lfc_down,
    'fig6-dan1dan2-rot-lfc-up' : fig6_dan1dan2_rot_lfc_up,
    'fig6-dan1dan2-rot-lfc-down' : fig6_dan1dan2_rot_lfc_down,
    'fig7-umaps' : fig7_umap_gwas,
    'gwas-hits':gwas_hits,
    'gwas-list': snps_gene_list,
    'poulin' : poulin
}


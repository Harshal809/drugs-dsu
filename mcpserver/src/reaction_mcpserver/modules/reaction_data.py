reaction_smarts_map = {
    # Basic Organic Reactions (Original)
    'Esterification': {
        'smarts': '[C:1][OH:2].[C:3](=[O])[OH:4]>>[C:1][O:2][C:3](=[O]).[OH2:4]',
        'description': 'Alcohol + Carboxylic Acid -> Ester + Water',
        'priority': 8
    },
    'AmideFormation': {
        'smarts': '[N:1].[C:2](=[O])[OH:3]>>[N:1][C:2](=[O]).[OH2:3]',
        'description': 'Amine + Carboxylic Acid -> Amide + Water',
        'priority': 8
    },
    'Hydrolysis': {
        'smarts': '[C:1][O:2][C:3](=[O])>>[C:1][OH:2].[C:3](=[O])[OH]',
        'description': 'Ester -> Alcohol + Carboxylic Acid',
        'priority': 6
    },
    'AmideHydrolysis': {
        'smarts': '[C:1][C:2](=[O])[N:3]>>[C:1][C:2](=[O])[OH].[N:3]',
        'description': 'Amide -> Carboxylic Acid + Amine',
        'priority': 6
    },
    
    # Nucleophilic Substitution Reactions (Original)
    'SN2_Alcohol': {
        'smarts': '[C:1][Cl,Br,I:2].[OH2:3]>>[C:1][OH:3].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Water -> Alcohol + Halide',
        'priority': 7
    },
    'SN2_Amine': {
        'smarts': '[C:1][Cl,Br,I:2].[N:3]>>[C:1][N:3].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Amine -> Substituted Amine + Halide',
        'priority': 7
    },
    'SN2_Thiol': {
        'smarts': '[C:1][Cl,Br,I:2].[S:3][H]>>[C:1][S:3].[Cl,Br,I:2][H]',
        'description': 'Alkyl Halide + Thiol -> Thioether + Hydrogen Halide',
        'priority': 7
    },
    'SN2_Alkoxide': {
        'smarts': '[C:1][Cl,Br,I:2].[C:3][O-:4]>>[C:1][O:4][C:3].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Alkoxide -> Ether + Halide',
        'priority': 7
    },
    
    # Addition Reactions (Original)
    'HalogenAddition': {
        'smarts': '[C:1]=[C:2].[Cl:3][Cl:4]>>[C:1]([Cl:3])[C:2][Cl:4]',
        'description': 'Alkene + Halogen -> Dihalide',
        'priority': 8
    },
    'HydrogenHalideAddition': {
        'smarts': '[C:1]=[C:2].[H:3][Cl,Br,I:4]>>[C:1]([H:3])[C:2][Cl,Br,I:4]',
        'description': 'Alkene + Hydrogen Halide -> Alkyl Halide',
        'priority': 8
    },
    'WaterAddition': {
        'smarts': '[C:1]=[C:2].[OH2:3]>>[C:1]([OH:3])[C:2][H]',
        'description': 'Alkene + Water -> Alcohol (Hydration)',
        'priority': 7
    },
    'HydrogenAddition': {
        'smarts': '[C:1]=[C:2].[H:3][H:4]>>[C:1]([H:3])[C:2][H:4]',
        'description': 'Alkene + Hydrogen -> Alkane (Hydrogenation)',
        'priority': 8
    },
    
    # Elimination Reactions (Original)
    'Dehydration': {
        'smarts': '[C:1][C:2]([OH:3])[H:4]>>[C:1]=[C:2].[OH2:3]',
        'description': 'Alcohol -> Alkene + Water (Dehydration)',
        'priority': 6
    },
    'Dehydrohalogenation': {
        'smarts': '[C:1][C:2]([Cl,Br,I:3])[H:4]>>[C:1]=[C:2].[H:4][Cl,Br,I:3]',
        'description': 'Alkyl Halide -> Alkene + Hydrogen Halide',
        'priority': 6
    },
    
    # Oxidation Reactions (Original)
    'AlcoholToAldehyde': {
        'smarts': '[C:1][CH2:2][OH:3]>>[C:1][C:2](=[O:3])[H]',
        'description': 'Primary Alcohol -> Aldehyde',
        'priority': 7
    },
    'AlcoholToKetone': {
        'smarts': '[C:1][CH:2]([OH:3])[C:4]>>[C:1][C:2](=[O:3])[C:4]',
        'description': 'Secondary Alcohol -> Ketone',
        'priority': 7
    },
    'AldehydeToAcid': {
        'smarts': '[C:1][C:2](=[O:3])[H:4]>>[C:1][C:2](=[O:3])[OH:4]',
        'description': 'Aldehyde -> Carboxylic Acid',
        'priority': 7
    },
    
    # Reduction Reactions (Original)
    'CarbonylToAlcohol': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[H:5][H:6]>>[C:1][C:2]([OH:3])[H:4]',
        'description': 'Aldehyde -> Primary Alcohol',
        'priority': 7
    },
    'KetoneToAlcohol': {
        'smarts': '[C:1][C:2](=[O:3])[C:4].[H:5][H:6]>>[C:1][C:2]([OH:3])[C:4]',
        'description': 'Ketone -> Secondary Alcohol',
        'priority': 7
    },
    'CarboxylicAcidToAlcohol': {
        'smarts': '[C:1][C:2](=[O:3])[OH:4].[H:5][H:6]>>[C:1][C:2][OH:4]',
        'description': 'Carboxylic Acid -> Primary Alcohol',
        'priority': 6
    },
    
    # Aromatic Substitution (Original)
    'Nitration': {
        'smarts': '[c:1][H:2].[N:3](=[O:4])(=[O:5])[OH:6]>>[c:1][N:3](=[O:4])(=[O:5]).[OH2:6]',
        'description': 'Aromatic -> Nitro Compound',
        'priority': 7
    },
    'Halogenation_Aromatic': {
        'smarts': '[c:1][H:2].[Cl:3][Cl:4]>>[c:1][Cl:3].[H:2][Cl:4]',
        'description': 'Aromatic Halogenation',
        'priority': 7
    },
    'Friedel_Crafts_Alkylation': {
        'smarts': '[c:1][H:2].[C:3][Cl:4]>>[c:1][C:3].[H:2][Cl:4]',
        'description': 'Friedel-Crafts Alkylation',
        'priority': 6
    },
    'Friedel_Crafts_Acylation': {
        'smarts': '[c:1][H:2].[C:3](=[O:4])[Cl:5]>>[c:1][C:3](=[O:4]).[H:2][Cl:5]',
        'description': 'Friedel-Crafts Acylation',
        'priority': 6
    },
    
    # Cross-Coupling Reactions (Original)
    'Suzuki_Coupling': {
        'smarts': '[c:1][Br,I:2].[c:3][B:4]([OH:5])[OH:6]>>[c:1][c:3].[Br,I:2][B:4]([OH:5])[OH:6]',
        'description': 'Suzuki Cross-Coupling',
        'priority': 9
    },
    'Heck_Reaction': {
        'smarts': '[c:1][Br,I:2].[C:3]=[C:4]>>[c:1][C:3]=[C:4].[Br,I:2]',
        'description': 'Heck Reaction',
        'priority': 8
    },
    'Sonogashira_Coupling': {
        'smarts': '[c:1][Br,I:2].[C:3]#[C:4][H:5]>>[c:1][C:3]#[C:4].[Br,I:2][H:5]',
        'description': 'Sonogashira Coupling',
        'priority': 8
    },
    
    # Cycloaddition Reactions (Original)
    'Diels_Alder': {
        'smarts': '[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2][C:5][C:6][C:4][C:3]1',
        'description': 'Diels-Alder Cycloaddition',
        'priority': 9
    },
    'Cyclopropanation': {
        'smarts': '[C:1]=[C:2].[C:3][H:4][Cl:5]>>[C:1]1[C:2][C:3]1.[H:4][Cl:5]',
        'description': 'Cyclopropanation',
        'priority': 7
    },
    
    # Rearrangement Reactions (Original)
    'Claisen_Rearrangement': {
        'smarts': '[C:1]=[C:2][C:3][O:4][C:5]=[C:6]>>[C:1][C:2]=[C:3][C:4](=[O])[C:5]=[C:6]',
        'description': 'Claisen Rearrangement',
        'priority': 6
    },
    'Beckmann_Rearrangement': {
        'smarts': '[C:1][C:2](=[N:3][OH:4])[C:5]>>[C:1][N:3][C:2](=[O])[C:5]',
        'description': 'Beckmann Rearrangement',
        'priority': 6
    },
    
    # Specialized Reactions (Original)
    'ThiocarbonylHydrolysis': {
        'smarts': '[C:1](=[S:2])>>[C:1](=[O:2]).[H][S][H]',
        'description': 'Thiocarbonyl -> Carbonyl + H2S',
        'priority': 6
    },
    'ThioureaFormation': {
        'smarts': '[N:1][H].[C:2]=[S:3]>>[H][N:1][C:2](=[S:3])[N:1][H]',
        'description': 'Amine + Thiocarbonyl -> Thiourea Derivative',
        'priority': 7
    },
    'SecondaryAmineThiourea': {
        'smarts': '[N:1]([H])[C:4].[C:2]=[S:3]>>[C:4][N:1][C:2](=[S:3])[N:1][H]',
        'description': 'Secondary Amine + Thiocarbonyl -> Substituted Thiourea',
        'priority': 7
    },
    'Dehalogenation': {
        'smarts': '[c:1][I:2].[H][O:3][H]>>[c:1][H].[I:2][O:3][H]',
        'description': 'Aryl Iodide + Water -> Aryl Hydrogen + Hypoiodous Acid',
        'priority': 5
    },
    
    # Condensation Reactions (Original)
    'Aldol_Condensation': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[C:5][C:6](=[O:7])[H:8]>>[C:1][C:2](=[O:3])[C:4]([OH:8])[C:5][C:6](=[O:7])',
        'description': 'Aldol Condensation',
        'priority': 7
    },
    'Mannich_Reaction': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[N:5].[C:6]=[O:7]>>[C:1][C:2](=[O:3])[C:4]([N:5])[C:6][OH:7]',
        'description': 'Mannich Reaction',
        'priority': 7
    },
    'Knoevenagel_Condensation': {
        'smarts': '[C:1]=[O:2].[C:3][C:4](=[O:5])[OH:6]>>[C:1]=[C:3][C:4](=[O:5])[OH:6].[OH2:2]',
        'description': 'Knoevenagel Condensation',
        'priority': 7
    },
    
    # Grignard Reactions (Original)
    'Grignard_Carbonyl': {
        'smarts': '[C:1][Mg:2][Br:3].[C:4]=[O:5]>>[C:1][C:4]([OH:5])[Mg:2][Br:3]',
        'description': 'Grignard + Carbonyl -> Alcohol',
        'priority': 8
    },
    'Grignard_CO2': {
        'smarts': '[C:1][Mg:2][Br:3].[C:4]=[O:5]=[O:6]>>[C:1][C:4](=[O:5])[OH:6].[Mg:2][Br:3]',
        'description': 'Grignard + CO2 -> Carboxylic Acid',
        'priority': 8
    },
    
    # Protecting Group Chemistry (Original)
    'Silyl_Protection': {
        'smarts': '[C:1][OH:2].[Si:3]([C:4])([C:5])[Cl:6]>>[C:1][O:2][Si:3]([C:4])([C:5]).[H:6][Cl]',
        'description': 'Alcohol Silyl Protection',
        'priority': 6
    },
    'Acetal_Formation': {
        'smarts': '[C:1]=[O:2].[C:3][OH:4].[C:5][OH:6]>>[C:1]([O:4][C:3])([O:6][C:5])[H:2]',
        'description': 'Aldehyde/Ketone Acetal Protection',
        'priority': 6
    },
    
    # Heterocycle Formation (Original)
    'Pyrrole_Formation': {
        'smarts': '[C:1](=[O:2])[C:3](=[O:4]).[N:5][H:6][H:7]>>[c:1]1[c:3][nH:5][c:1][c:3]1',
        'description': 'Diketone + Ammonia -> Pyrrole',
        'priority': 7
    },
    'Imidazole_Formation': {
        'smarts': '[C:1](=[O:2])[C:3](=[O:4]).[N:5][H:6][H:7].[C:8]=[O:9]>>[c:1]1[n:5][c:8][n:5][c:3]1',
        'description': 'Diketone + NH3 + Aldehyde -> Imidazole',
        'priority': 7
    },
    
    # Multi-component Reactions (Original)
    'Ugi_Reaction': {
        'smarts': '[C:1]=[O:2].[N:3].[C:4](=[O:5])[OH:6].[C:7][N:8]#[C:9]>>[C:1]([N:3][C:4](=[O:5])[N:8]([C:7])[C:9])[OH:2]',
        'description': 'Ugi Four-Component Reaction',
        'priority': 8
    },
    'Biginelli_Reaction': {
        'smarts': '[C:1]=[O:2].[C:3][C:4](=[O:5])[C:6].[N:7][H:8][H:9]>>[C:1]1[N:7][C:4](=[O:5])[C:6][C:3][N:7]1',
        'description': 'Biginelli Multicomponent Reaction',
        'priority': 7
    },
    
    # Previous Additions (from first response)
    'SN2_Cyanide': {
        'smarts': '[C:1][Cl,Br,I:2].[C:3]#[N:4]>>[C:1][C:3]#[N:4].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Cyanide -> Nitrile + Halide',
        'priority': 7
    },
    'SN2_Azide': {
        'smarts': '[C:1][Cl,Br,I:2].[N:3]=[N+:4]=[N-:5]>>[C:1][N:3]=[N+:4]=[N-:5].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Azide -> Alkyl Azide + Halide',
        'priority': 7
    },
    'SNAr_Fluoride': {
        'smarts': '[c:1][F:2].[N:3]>>[c:1][N:3].[F:2]',
        'description': 'Aryl Fluoride + Amine -> Aryl Amine + Fluoride (SNAr)',
        'priority': 6
    },
    'SNAr_Alkoxide': {
        'smarts': '[c:1][Cl,Br:2].[O-:3][C:4]>>[c:1][O:3][C:4].[Cl,Br:2]',
        'description': 'Aryl Halide + Alkoxide -> Aryl Ether + Halide (SNAr)',
        'priority': 6
    },
    'Hydroboration': {
        'smarts': '[C:1]=[C:2].[B:3][H:4]>>[C:1][C:2][B:3][H:4]',
        'description': 'Alkene + Borane -> Alkylborane',
        'priority': 7
    },
    'Epoxidation': {
        'smarts': '[C:1]=[C:2].[O:3]>>[C:1]1[C:2][O:3]1',
        'description': 'Alkene + Peroxide -> Epoxide',
        'priority': 7
    },
    'Dihydroxylation': {
        'smarts': '[C:1]=[C:2].[O:3][O:4]>>[C:1]([OH:3])[C:2][OH:4]',
        'description': 'Alkene + OsO4 -> 1,2-Diol',
        'priority': 7
    },
    'Ozonolysis': {
        'smarts': '[C:1]=[C:2].[O:3][O:4][O:5]>>[C:1]=[O:3].[C:2]=[O:4]',
        'description': 'Alkene + Ozone -> Two Carbonyls',
        'priority': 6
    },
    'E2_Elimination': {
        'smarts': '[C:1][C:2]([Cl,Br,I:3])[H:4].[O-:5]>>[C:1]=[C:2].[Cl,Br,I:3][H:4][O:5]',
        'description': 'Alkyl Halide + Base -> Alkene + HX',
        'priority': 6
    },
    'Dehalogenation_Elimination': {
        'smarts': '[C:1]([Cl,Br:2])[C:3]([Cl,Br:4])>>[C:1]=[C:3].[Cl,Br:2][Cl,Br:4]',
        'description': 'Vicinal Dihalide -> Alkene + Halogen',
        'priority': 6
    },
    'AlcoholToAcid': {
        'smarts': '[C:1][CH2:2][OH:3]>>[C:1][C:2](=[O])[OH:3]',
        'description': 'Primary Alcohol -> Carboxylic Acid',
        'priority': 7
    },
    'AlkeneToEpoxide': {
        'smarts': '[C:1]=[C:2].[O:3]>>[C:1]1[C:2][O:3]1',
        'description': 'Alkene -> Epoxide (mCPBA)',
        'priority': 7
    },
    'SulfideToSulfoxide': {
        'smarts': '[C:1][S:2][C:3].[O:4]>>[C:1][S:2](=[O:4])[C:3]',
        'description': 'Sulfide -> Sulfoxide',
        'priority': 6
    },
    'SulfoxideToSulfone': {
        'smarts': '[C:1][S:2](=[O:3])[C:4].[O:5]>>[C:1][S:2](=[O:3])(=[O:5])[C:4]',
        'description': 'Sulfoxide -> Sulfone',
        'priority': 6
    },
    'NitroToAmine': {
        'smarts': '[c:1][N+:2](=[O:3])[O-:4].[H:5][H:6]>>[c:1][N:2].[O:3][H:5][O:4][H:6]',
        'description': 'Nitro -> Amine',
        'priority': 7
    },
    'NitrileToAmine': {
        'smarts': '[C:1]#[N:2].[H:3][H:4]>>[C:1][CH2:2][NH2:2]',
        'description': 'Nitrile -> Primary Amine',
        'priority': 7
    },
    'ImineToAmine': {
        'smarts': '[C:1]=[N:2].[H:3][H:4]>>[C:1][NH2:2]',
        'description': 'Imine -> Amine',
        'priority': 7
    },
    'Negishi_Coupling': {
        'smarts': '[c:1][Br,I:2].[C:3][Zn:4][Cl:5]>>[c:1][C:3].[Br,I:2][Zn:4][Cl:5]',
        'description': 'Negishi Cross-Coupling',
        'priority': 8
    },
    'Stille_Coupling': {
        'smarts': '[c:1][Br,I:2].[C:3][Sn:4]>>[c:1][C:3].[Br,I:2][Sn:4]',
        'description': 'Stille Cross-Coupling',
        'priority': 8
    },
    'Kumada_Coupling': {
        'smarts': '[c:1][Br,I:2].[C:3][Mg:4][Cl:5]>>[c:1][C:3].[Br,I:2][Mg:4][Cl:5]',
        'description': 'Kumada Cross-Coupling',
        'priority': 8
    },
    '[2+2]_Cycloaddition': {
        'smarts': '[C:1]=[C:2].[C:3]=[C:4]>>[C:1]1[C:2][C:3][C:4]1',
        'description': '[2+2] Cycloaddition',
        'priority': 7
    },
    '1,3-Dipolar_Cycloaddition': {
        'smarts': '[C:1]=[C:2].[C:3]=[N+:4]=[N-:5]>>[C:1]1[C:2][C:3][N+:4]=[N-:5]1',
        'description': '1,3-Dipolar Cycloaddition with Azide',
        'priority': 7
    },
    'Pinacol_Rearrangement': {
        'smarts': '[C:1]([OH:2])([C:3])[C:4]([OH:5])[C:6]>>[C:1](=[O:2])[C:3][C:4][C:6].[OH2:5]',
        'description': 'Pinacol Rearrangement',
        'priority': 6
    },
    'Wagner_Meerwein': {
        'smarts': '[C:1][C:2]([OH:3])[C+:4][C:5]>>[C:1]=[C:2][C:4][C:5].[OH2:3]',
        'description': 'Wagner-Meerwein Rearrangement',
        'priority': 6
    },
    'Pyrazole_Formation': {
        'smarts': '[C:1](=[O:2])[C:3](=[O:4]).[N:5][N:6][H:7]>>[c:1]1[c:3][n:5][n:6][c:1]1',
        'description': 'Diketone + Hydrazine -> Pyrazole',
        'priority': 7
    },
    'Triazole_Formation': {
        'smarts': '[C:1]#[C:2].[N:3]=[N+:4]=[N-:5]>>[c:1]1[c:2][n:3][n+:4]=[n-:5]1',
        'description': 'Alkyne + Azide -> 1,2,3-Triazole (Click Chemistry)',
        'priority': 8
    },
    'Boc_Protection': {
        'smarts': '[N:1][H:2].[C:3](=[O:4])[O:5][C:6]>>[N:1][C:3](=[O:4])[O:5][C:6].[H:2]',
        'description': 'Amine + Boc Anhydride -> Boc-Protected Amine',
        'priority': 6
    },
    'Cbz_Protection': {
        'smarts': '[N:1][H:2].[C:3](=[O:4])[O:5][c:6]>>[N:1][C:3](=[O:4])[O:5][c:6].[H:2]',
        'description': 'Amine + Cbz Chloride -> Cbz-Protected Amine',
        'priority': 6
    },
    'Claisen_Condensation': {
        'smarts': '[C:1][C:2](=[O:3])[O:4][C:5].[C:6][C:7](=[O:8])[O:9][C:10]>>[C:1][C:2](=[O:3])[C:6][C:7](=[O:8])[OH:9].[C:5][O:4][C:10]',
        'description': 'Claisen Condensation of Esters',
        'priority': 7
    },
    'Perkin_Reaction': {
        'smarts': '[c:1][C:2](=[O:3])[H:4].[C:5][C:6](=[O:7])[OH:8]>>[c:1][C:2]=[C:5][C:6](=[O:7])[OH:8].[OH2:3]',
        'description': 'Perkin Reaction',
        'priority': 7
    },
    'Passerini_Reaction': {
        'smarts': '[C:1]=[O:2].[C:3](=[O:4])[OH:5].[C:6][N:7]#[C:8]>>[C:1]([O:4][C:3])[O:5][C:6][C:8][N:7][H:2]',
        'description': 'Passerini Three-Component Reaction',
        'priority': 7
    },
    'Hantzsch_Pyridine': {
        'smarts': '[C:1]=[O:2].[C:3][C:4](=[O:5])[C:6].[N:7][H:8][H:9]>>[c:1]1[c:3][c:4][c:6][n:7][c:1]1',
        'description': 'Hantzsch Pyridine Synthesis',
        'priority': 7
    },
    
    # New Reactions: Additional Nucleophilic Substitutions
    'SN1_Alcohol': {
        'smarts': '[C:1][Cl,Br,I:2].[OH2:3]>>[C:1][OH:3].[Cl,Br,I:2]',
        'description': 'Tertiary Alkyl Halide + Water -> Alcohol + Halide (SN1)',
        'priority': 6
    },
    'SN1_Ether': {
        'smarts': '[C:1][Cl,Br,I:2].[C:3][OH:4]>>[C:1][O:4][C:3].[Cl,Br,I:2]',
        'description': 'Tertiary Alkyl Halide + Alcohol -> Ether + Halide (SN1)',
        'priority': 6
    },
    'SNAr_Nitro': {
        'smarts': '[c:1][Cl:2].[N:3][H:4]>>[c:1][N:3].[Cl:2][H:4]',
        'description': 'Nitroaryl Chloride + Amine -> Nitroaryl Amine (SNAr)',
        'priority': 6
    },
    'SN2_Phosphine': {
        'smarts': '[C:1][Cl,Br,I:2].[P:3][H:4]>>[C:1][P:3].[Cl,Br,I:2][H:4]',
        'description': 'Alkyl Halide + Phosphine -> Phosphonium + Halide',
        'priority': 6
    },
    'SN2_Acetylide': {
        'smarts': '[C:1][Cl,Br,I:2].[C:3]#[C:4]>>[C:1][C:3]#[C:4].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Acetylide -> Alkyne + Halide',
        'priority': 7
    },
    'SN2_Hydrosulfide': {
        'smarts': '[C:1][Cl,Br,I:2].[S-:3][H:4]>>[C:1][S:3][H:4].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Hydrosulfide -> Thiol + Halide',
        'priority': 7
    },
    'SNAr_Thiolate': {
        'smarts': '[c:1][Cl,Br:2].[S-:3][C:4]>>[c:1][S:3][C:4].[Cl,Br:2]',
        'description': 'Aryl Halide + Thiolate -> Aryl Thioether + Halide (SNAr)',
        'priority': 6
    },
    'SN2_Enolate': {
        'smarts': '[C:1][Cl,Br,I:2].[C:3][C:4](=[O:5])>>[C:1][C:3][C:4](=[O:5]).[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Enolate -> Alkylated Carbonyl + Halide',
        'priority': 7
    },
    'SN2_Carboxylate': {
        'smarts': '[C:1][Cl,Br,I:2].[C:3](=[O:4])[O-:5]>>[C:1][O:5][C:3](=[O:4]).[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Carboxylate -> Ester + Halide',
        'priority': 7
    },
    'SN2_Hydroxide': {
        'smarts': '[C:1][Cl,Br,I:2].[OH-:3]>>[C:1][OH:3].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Hydroxide -> Alcohol + Halide',
        'priority': 7
    },
    
    # New Reactions: Additional Addition Reactions
    'Michael_Addition': {
        'smarts': '[C:1]=[C:2][C:3](=[O:4]).[C:5][H:6]>>[C:1][C:2]([C:3](=[O:4]))[C:5][H:6]',
        'description': 'Alkene + Nucleophile -> Michael Addition Product',
        'priority': 7
    },
    'Hydroformylation': {
        'smarts': '[C:1]=[C:2].[C:3]=[O:4].[H:5]>>[C:1][C:2][C:3](=[O:4])[H:5]',
        'description': 'Alkene + CO + H2 -> Aldehyde (Hydroformylation)',
        'priority': 7
    },
    'Hydrosilylation': {
        'smarts': '[C:1]=[C:2].[Si:3][H:4]>>[C:1][C:2][Si:3][H:4]',
        'description': 'Alkene + Silane -> Alkylsilane',
        'priority': 6
    },
    'Oxymercuration': {
        'smarts': '[C:1]=[C:2].[Hg:3][O:4][C:5]>>[C:1][OH:4][C:2][Hg:3][O:5][C:5]',
        'description': 'Alkene + Mercury Acetate -> Alcohol (Oxymercuration)',
        'priority': 6
    },
    'Hydroamination': {
        'smarts': '[C:1]=[C:2].[N:3][H:4]>>[C:1][C:2][N:3][H:4]',
        'description': 'Alkene + Amine -> Alkylamine',
        'priority': 6
    },
    'Hydrocyanation': {
        'smarts': '[C:1]=[C:2].[C:3]#[N:4]>>[C:1][C:2][C:3]#[N:4]',
        'description': 'Alkene + HCN -> Nitrile',
        'priority': 6
    },
    'Alkyne_Hydration': {
        'smarts': '[C:1]#[C:2].[OH2:3]>>[C:1][C:2](=[O:3])[H]',
        'description': 'Alkyne + Water -> Ketone (Hydration)',
        'priority': 7
    },
    'Bromohydrin_Formation': {
        'smarts': '[C:1]=[C:2].[Br:3][OH2:4]>>[C:1][Br:3][C:2][OH:4]',
        'description': 'Alkene + Br2/H2O -> Bromohydrin',
        'priority': 7
    },
    'Carbene_Addition': {
        'smarts': '[C:1]=[C:2].[C:3][Cl:4][Cl:5]>>[C:1]1[C:2][C:3]1.[Cl:4][Cl:5]',
        'description': 'Alkene + Dichlorocarbene -> Cyclopropane',
        'priority': 7
    },
    'Alkyne_Halogenation': {
        'smarts': '[C:1]#[C:2].[Cl:3][Cl:4]>>[C:1]([Cl:3])=[C:2][Cl:4]',
        'description': 'Alkyne + Halogen -> Dihaloalkene',
        'priority': 7
    },
    
    # New Reactions: Additional Elimination Reactions
    'E1_Elimination': {
        'smarts': '[C:1][C:2]([Cl,Br,I:3])[H:4]>>[C:1]=[C:2].[Cl,Br,I:3][H:4]',
        'description': 'Tertiary Alkyl Halide -> Alkene + HX (E1)',
        'priority': 6
    },
    'Dehydrohalogenation_Alkyne': {
        'smarts': '[C:1]([Cl,Br:2])[C:3]([Cl,Br:4])[H:5]>>[C:1]#[C:3].[Cl,Br:2][Cl,Br:4][H:5]',
        'description': 'Vicinal Dihalide -> Alkyne + HX',
        'priority': 6
    },
    'Dehydration_Ketone': {
        'smarts': '[C:1][C:2]([OH:3])([C:4])[H:5]>>[C:1]=[C:2][C:4].[OH2:3]',
        'description': 'Tertiary Alcohol -> Alkene + Water (Dehydration to Ketone)',
        'priority': 6
    },
    'Elimination_AmineOxide': {
        'smarts': '[C:1][C:2][N+:3]([C:4])[O-:5]>>[C:1]=[C:2].[N:3][C:4]',
        'description': 'Amine Oxide -> Alkene + Amine (Cope Elimination)',
        'priority': 6
    },
    'Dehydrofluorination': {
        'smarts': '[C:1][C:2]([F:3])[H:4]>>[C:1]=[C:2].[F:3][H:4]',
        'description': 'Alkyl Fluoride -> Alkene + HF',
        'priority': 6
    },
    
    # New Reactions: Additional Oxidation Reactions
    'Baeyer_Villiger': {
        'smarts': '[C:1][C:2](=[O:3])[C:4]>>[C:1][O:3][C:2](=[O])[C:4]',
        'description': 'Ketone -> Ester (Baeyer-Villiger)',
        'priority': 7
    },
    'Swern_Oxidation': {
        'smarts': '[C:1][CH2:2][OH:3]>>[C:1][C:2](=[O:3])[H]',
        'description': 'Primary Alcohol -> Aldehyde (Swern)',
        'priority': 7
    },
    'Jones_Oxidation': {
        'smarts': '[C:1][CH:2]([OH:3])[C:4]>>[C:1][C:2](=[O:3])[C:4]',
        'description': 'Secondary Alcohol -> Ketone (Jones)',
        'priority': 7
    },
    'Alkene_Cleavage': {
        'smarts': '[C:1]=[C:2].[O:3][O:4]>>[C:1]=[O:3].[C:2]=[O:4]',
        'description': 'Alkene -> Two Carbonyls (Oxidative Cleavage)',
        'priority': 6
    },
    'ThiolToDisulfide': {
        'smarts': '[C:1][S:2][H:3].[C:4][S:5][H:6]>>[C:1][S:2][S:5][C:4].[H:3][H:6]',
        'description': 'Thiol -> Disulfide',
        'priority': 6
    },
    'AmineToNitroso': {
        'smarts': '[C:1][N:2][H:3].[O:4]>>[C:1][N:2]=[O:4].[H:3]',
        'description': 'Amine -> Nitroso Compound',
        'priority': 6
    },
    'AlkyneToDiketone': {
        'smarts': '[C:1]#[C:2].[O:3][O:4]>>[C:1](=[O:3])[C:2](=[O:4])',
        'description': 'Alkyne -> 1,2-Diketone (Oxidation)',
        'priority': 6
    },
    'AlcoholToEster': {
        'smarts': '[C:1][OH:2].[C:3](=[O:4])[Cl:5]>>[C:1][O:2][C:3](=[O:4]).[Cl:5][H]',
        'description': 'Alcohol + Acid Chloride -> Ester + HCl',
        'priority': 7
    },
    'AldehydeToAcetal': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[O:5][H:6]>>[C:1][C:2]([O:3][H:6])[O:5][H:4]',
        'description': 'Aldehyde + Alcohol -> Acetal',
        'priority': 6
    },
    'KetoneToAcetal': {
        'smarts': '[C:1][C:2](=[O:3])[C:4].[O:5][H:6]>>[C:1][C:2]([O:3][H:6])[O:5][C:4]',
        'description': 'Ketone + Alcohol -> Acetal',
        'priority': 6
    },
    
    # New Reactions: Additional Reduction Reactions
    'AzideToAmine': {
        'smarts': '[C:1][N:2]=[N+:3]=[N-:4].[H:5][H:6]>>[C:1][NH2:2].[N:3]#[N:4]',
        'description': 'Azide -> Amine + Nitrogen',
        'priority': 7
    },
    'AlkyneToAlkene': {
        'smarts': '[C:1]#[C:2].[H:3][H:4]>>[C:1]=[C:2][H:3][H:4]',
        'description': 'Alkyne -> Alkene (Partial Hydrogenation)',
        'priority': 7
    },
    'EsterToAlcohol': {
        'smarts': '[C:1][O:2][C:3](=[O:4]).[H:5][H:6]>>[C:1][OH:2].[C:3][OH:4]',
        'description': 'Ester -> Two Alcohols',
        'priority': 7
    },
    'AmideToAmine': {
        'smarts': '[C:1][C:2](=[O:3])[N:4].[H:5][H:6]>>[C:1][CH2:2][NH2:4]',
        'description': 'Amide -> Amine (Reduction)',
        'priority': 7
    },
    'NitrosoToAmine': {
        'smarts': '[C:1][N:2]=[O:3].[H:4][H:5]>>[C:1][NH2:2].[OH2:3]',
        'description': 'Nitroso -> Amine',
        'priority': 6
    },
    'AlkyneToAlkane': {
        'smarts': '[C:1]#[C:2].[H:3][H:4][H:5][H:6]>>[C:1][C:2][H:3][H:4][H:5][H:6]',
        'description': 'Alkyne -> Alkane (Full Hydrogenation)',
        'priority': 7
    },
    'EpoxideToAlcohol': {
        'smarts': '[C:1]1[C:2][O:3]1.[H:4][H:5]>>[C:1][C:2][OH:3][H:4]',
        'description': 'Epoxide -> Alcohol (Reduction)',
        'priority': 6
    },
    'DisulfideToThiol': {
        'smarts': '[C:1][S:2][S:3][C:4].[H:5][H:6]>>[C:1][S:2][H:5].[C:4][S:3][H:6]',
        'description': 'Disulfide -> Two Thiols',
        'priority': 6
    },
    'NitrileToAldehyde': {
        'smarts': '[C:1]#[N:2].[H:3][H:4]>>[C:1][C:2](=[O:3])[H:4]',
        'description': 'Nitrile -> Aldehyde (Partial Reduction)',
        'priority': 6
    },
    'KetoneToAlkane': {
        'smarts': '[C:1][C:2](=[O:3])[C:4].[H:5][H:6]>>[C:1][CH2:2][C:4]',
        'description': 'Ketone -> Alkane (Clemmensen/Wolff-Kishner)',
        'priority': 6
    },
    
    # New Reactions: Additional Cross-Coupling Reactions
    'Buchwald_Hartwig': {
        'smarts': '[c:1][Br,I:2].[N:3][H:4]>>[c:1][N:3].[Br,I:2][H:4]',
        'description': 'Aryl Halide + Amine -> Aryl Amine (Buchwald-Hartwig)',
        'priority': 8
    },
    'Hiyama_Coupling': {
        'smarts': '[c:1][Br,I:2].[C:3][Si:4]>>[c:1][C:3].[Br,I:2][Si:4]',
        'description': 'Aryl Halide + Organosilane -> Coupled Product (Hiyama)',
        'priority': 8
    },
    'Miyaura_Borylation': {
        'smarts': '[c:1][Br,I:2].[B:3]([O:4])[O:5]>>[c:1][B:3]([O:4])[O:5].[Br,I:2]',
        'description': 'Aryl Halide -> Aryl Boronate (Miyaura Borylation)',
        'priority': 8
    },
    'Cross_Coupling_Alkyne': {
        'smarts': '[c:1][Br,I:2].[C:3]#[C:4]>>[c:1][C:3]#[C:4].[Br,I:2]',
        'description': 'Aryl Halide + Alkyne -> Aryl Alkyne',
        'priority': 8
    },
    'Suzuki_Coupling_Alkyl': {
        'smarts': '[C:1][Br,I:2].[C:3][B:4]([OH:5])[OH:6]>>[C:1][C:3].[Br,I:2][B:4]([OH:5])[OH:6]',
        'description': 'Alkyl Halide + Boronic Acid -> Coupled Product',
        'priority': 7
    },
    
    # New Reactions: Additional Cycloadditions
    '1,3-Dipolar_NitrileOxide': {
        'smarts': '[C:1]=[C:2].[C:3]=[N+:4][O-:5]>>[C:1]1[C:2][C:3][N+:4][O-:5]1',
        'description': 'Alkene + Nitrile Oxide -> Isoxazoline',
        'priority': 7
    },
    '[4+2]_Hetero_Diels_Alder': {
        'smarts': '[C:1]=[C:2][C:3]=[C:4].[O:5]=[C:6]>>[C:1]1[C:2][C:6][O:5][C:4][C:3]1',
        'description': 'Hetero Diels-Alder with Carbonyl',
        'priority': 7
    },
    'Paterno_Buchi': {
        'smarts': '[C:1]=[C:2].[C:3]=[O:4]>>[C:1]1[C:2][O:4][C:3]1',
        'description': 'Alkene + Carbonyl -> Oxetane (Paterno-Buchi)',
        'priority': 7
    },
    'Azide_Alkyne_Cycloaddition': {
        'smarts': '[C:1]#[C:2].[N:3]=[N+:4]=[N-:5]>>[c:1]1[c:2][n:3][n+:4]=[n-:5]1',
        'description': 'Alkyne + Azide -> Triazole (Click Chemistry)',
        'priority': 8
    },
    '[3+2]_Cycloaddition': {
        'smarts': '[C:1]=[C:2].[C:3][C:4][C:5]>>[C:1]1[C:2][C:3][C:4][C:5]1',
        'description': 'Alkene + 1,3-Dipole -> Five-Membered Ring',
        'priority': 7
    },
    
    # New Reactions: Additional Rearrangements
    'Curtius_Rearrangement': {
        'smarts': '[C:1](=[O:2])[N:3][N:4][N:5]>>[C:1][N:3]#[C:4].[N:5]=[O:2]',
        'description': 'Acyl Azide -> Isocyanate (Curtius)',
        'priority': 6
    },
    'Hofmann_Rearrangement': {
        'smarts': '[C:1][C:2](=[O:3])[N:4][H:5].[Br:6][Br:7]>>[C:1][N:4][H:5].[C:2](=[O:3])[O:6][Br:7]',
        'description': 'Amide -> Amine (Hofmann)',
        'priority': 6
    },
    'Cope_Rearrangement': {
        'smarts': '[C:1]=[C:2][C:3][C:4]=[C:5][C:6]>>[C:1][C:2]=[C:3][C:4][C:5]=[C:6]',
        'description': 'Cope Rearrangement',
        'priority': 6
    },
    'Overman_Rearrangement': {
        'smarts': '[C:1]=[C:2][N:3][C:4](=[O:5])>>[C:1][C:2]=[N:3][C:4](=[O:5])',
        'description': 'Allylic Amide -> Allylic Amine (Overman)',
        'priority': 6
    },
    'Stevens_Rearrangement': {
        'smarts': '[C:1][N+:2]([C:3])[C:4][H:5]>>[C:1][C:4][N:2][C:3][H:5]',
        'description': 'Quaternary Ammonium -> Rearranged Amine (Stevens)',
        'priority': 6
    },
    
    # New Reactions: Additional Heterocycle Formations
    'Oxazole_Formation': {
        'smarts': '[C:1](=[O:2])[C:3][N:4][H:5]>>[c:1]1[o:2][c:3][n:4][c:1]1',
        'description': 'Amide + α-Haloketone -> Oxazole',
        'priority': 7
    },
    'Thiazole_Formation': {
        'smarts': '[C:1](=[O:2])[C:3][S:4][H:5]>>[c:1]1[s:4][c:3][n:2][c:1]1',
        'description': 'Thioamide + α-Haloketone -> Thiazole',
        'priority': 7
    },
    'Pyrimidine_Formation': {
        'smarts': '[C:1](=[O:2])[C:3](=[O:4]).[N:5][H:6][H:7]>>[c:1]1[n:5][c:3][n:5][c:1]1',
        'description': 'Diketone + Urea -> Pyrimidine',
        'priority': 7
    },
    'Isoxazole_Formation': {
        'smarts': '[C:1]=[C:2].[C:3]=[N+:4][O-:5]>>[c:1]1[c:2][c:3][n+:4][o:5]1',
        'description': 'Alkene + Nitrile Oxide -> Isoxazole',
        'priority': 7
    },
    'Indazole_Formation': {
        'smarts': '[c:1][C:2](=[O:3]).[N:4][N:5][H:6]>>[c:1]1[c:2][n:4][n:5][c:1]1',
        'description': 'Aryl Ketone + Hydrazine -> Indazole',
        'priority': 7
    },
    
    # New Reactions: Additional Protecting Group Chemistry
    'Fmoc_Protection': {
        'smarts': '[N:1][H:2].[C:3](=[O:4])[O:5][c:6]>>[N:1][C:3](=[O:4])[O:5][c:6].[H:2]',
        'description': 'Amine + Fmoc Chloride -> Fmoc-Protected Amine',
        'priority': 6
    },
    'Tosylate_Formation': {
        'smarts': '[C:1][OH:2].[S:3](=[O:4])(=[O:5])[Cl:6]>>[C:1][O:2][S:3](=[O:4])(=[O:5]).[Cl:6][H]',
        'description': 'Alcohol + Tosyl Chloride -> Tosylate',
        'priority': 6
    },
    'THP_Protection': {
        'smarts': '[C:1][OH:2].[O:3]1[C:4][C:5][C:6][C:7]1>>[C:1][O:2][O:3]1[C:4][C:5][C:6][C:7]1',
        'description': 'Alcohol + Dihydropyran -> THP-Protected Alcohol',
        'priority': 6
    },
    'MOM_Protection': {
        'smarts': '[C:1][OH:2].[C:3][O:4][Cl:5]>>[C:1][O:2][O:4][C:3].[Cl:5][H]',
        'description': 'Alcohol + MOM Chloride -> MOM-Protected Alcohol',
        'priority': 6
    },
    'Acetyl_Protection': {
        'smarts': '[C:1][OH:2].[C:3](=[O:4])[Cl:5]>>[C:1][O:2][C:3](=[O:4]).[Cl:5][H]',
        'description': 'Alcohol + Acetyl Chloride -> Acetate',
        'priority': 6
    },
    
    # New Reactions: Additional Condensation Reactions
    'Dieckmann_Condensation': {
        'smarts': '[C:1][C:2](=[O:3])[O:4][C:5][C:6](=[O:7])[O:8]>>[C:1][C:2]1[C:5][C:6](=[O:7])[O:8]1.[O:4][H]',
        'description': 'Dieckmann Condensation of Diesters',
        'priority': 7
    },
    'Benzoin_Condensation': {
        'smarts': '[c:1][C:2](=[O:3])[H:4].[c:5][C:6](=[O:7])[H:8]>>[c:1][C:2](=[O:3])[C:6]([OH:8])[c:5]',
        'description': 'Benzoin Condensation of Aromatic Aldehydes',
        'priority': 7
    },
    'Knoevenagel_Malonic': {
        'smarts': '[C:1]=[O:2].[C:3]([C:4](=[O:5])[OH:6])[C:7](=[O:8])[OH:9]>>[C:1]=[C:3][C:4](=[O:5])[OH:6].[C:7](=[O:8])[OH:9]',
        'description': 'Knoevenagel with Malonic Acid',
        'priority': 7
    },
    'Robinson_Annulation': {
        'smarts': '[C:1][C:2](=[O:3]).[C:4]=[C:5][C:6](=[O:7])>>[C:1]1[C:2](=[O:3])[C:4][C:5][C:6]1',
        'description': 'Robinson Annulation',
        'priority': 7
    },
    'Condensation_Imine': {
        'smarts': '[C:1]=[O:2].[N:3][H:4]>>[C:1]=[N:3].[OH2:2]',
        'description': 'Carbonyl + Amine -> Imine + Water',
        'priority': 7
    },
    
    # New Reactions: Additional Multi-component Reactions
    'Strecker_Synthesis': {
        'smarts': '[C:1]=[O:2].[N:3][H:4][H:5].[C:6]#[N:7]>>[C:1]([N:3][H:5])[C:6]#[N:7].[OH2:2]',
        'description': 'Strecker Amino Acid Synthesis',
        'priority': 7
    },
    'Gewald_Reaction': {
        'smarts': '[C:1]=[O:2].[C:3][C:4](=[O:5]).[S:6]>>[c:1]1[c:3][c:4][s:6][c:1]1',
        'description': 'Gewald Thiophene Synthesis',
        'priority': 7
    },
    'Paal_Knorr_Pyrrole': {
        'smarts': '[C:1](=[O:2])[C:3][C:4](=[O:5]).[N:6][H:7]>>[c:1]1[c:3][c:4][n:6][c:1]1',
        'description': 'Paal-Knorr Pyrrole Synthesis',
        'priority': 7
    },
    'Hantzsch_Dihydropyridine': {
        'smarts': '[C:1]=[O:2].[C:3][C:4](=[O:5])[C:6].[N:7][H:8]>>[c:1]1[c:3][c:4][c:6][n:7][c:1]1',
        'description': 'Hantzsch Dihydropyridine Synthesis',
        'priority': 7
    },
    'Biginelli_Pyrimidine': {
        'smarts': '[C:1]=[O:2].[C:3][C:4](=[O:5])[C:6].[N:7][C:8](=[O:9])[N:10]>>[c:1]1[n:7][c:4][c:6][n:10][c:1]1',
        'description': 'Biginelli Pyrimidine Synthesis with Urea',
        'priority': 7
    },
    
    # New Reactions: Specialized Reactions
    'Olefin_Metathesis': {
        'smarts': '[C:1]=[C:2].[C:3]=[C:4]>>[C:1]=[C:3].[C:2]=[C:4]',
        'description': 'Olefin Metathesis',
        'priority': 8
    },
    'CH_Activation': {
        'smarts': '[c:1][H:2].[C:3][Br:4]>>[c:1][C:3].[H:2][Br:4]',
        'description': 'Aryl C-H Activation with Alkyl Halide',
        'priority': 7
    },
    'Wacker_Oxidation': {
        'smarts': '[C:1]=[C:2][C:3].[O:4]>>[C:1][C:2](=[O:4])[C:3]',
        'description': 'Alkene -> Ketone (Wacker Oxidation)',
        'priority': 7
    },
    'Pauson_Khand': {
        'smarts': '[C:1]#[C:2].[C:3]=[C:4].[C:5]=[O:6]>>[C:1]1[C:2][C:3][C:4][C:5](=[O:6])1',
        'description': 'Pauson-Khand Cycloaddition',
        'priority': 7
    },
    'Grubbs_Metathesis': {
        'smarts': '[C:1]=[C:2][C:3].[C:4]=[C:5]>>[C:1]=[C:4].[C:2][C:3]=[C:5]',
        'description': 'Grubbs Olefin Metathesis',
        'priority': 8
    },
     # Additional Nucleophilic Substitution Reactions (10)
    'SN2_Alkylthiolate': {
        'smarts': '[C:1][Cl,Br,I:2].[S-:3][C:4]>>[C:1][S:3][C:4].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Thiolate -> Thioether + Halide',
        'priority': 7
    },
    'SN2_Alkynide': {
        'smarts': '[C:1][Cl,Br,I:2].[C:3]#[C-:4]>>[C:1][C:3]#[C:4].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Alkynide -> Terminal Alkyne + Halide',
        'priority': 7
    },
    'SNAr_Amine': {
        'smarts': '[c:1][Cl,Br:2].[N:3][H:4]>>[c:1][N:3].[Cl,Br:2][H:4]',
        'description': 'Aryl Halide + Amine -> Aryl Amine + Halide (SNAr)',
        'priority': 6
    },
    'SN2_Phosphite': {
        'smarts': '[C:1][Cl,Br,I:2].[P:3]([O:4])[O:5]>>[C:1][P:3]([O:4])[O:5].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Phosphite -> Phosphonate + Halide',
        'priority': 6
    },
    'SN2_Selenol': {
        'smarts': '[C:1][Cl,Br,I:2].[Se:3][H:4]>>[C:1][Se:3].[Cl,Br,I:2][H:4]',
        'description': 'Alkyl Halide + Selenol -> Selenoether + Hydrogen Halide',
        'priority': 7
    },
    'SN2_Alkylamine_Secondary': {
        'smarts': '[C:1][Cl,Br,I:2].[N:3]([C:4])[H:5]>>[C:1][N:3][C:4].[Cl,Br,I:2][H:5]',
        'description': 'Alkyl Halide + Secondary Amine -> Tertiary Amine + Halide',
        'priority': 7
    },
    'SN2_Carboxylate_Ester': {
        'smarts': '[C:1][Cl,Br,I:2].[O-:3][C:4](=[O:5])>>[C:1][O:3][C:4](=[O:5]).[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Carboxylate -> Ester + Halide',
        'priority': 7
    },
    'SN1_Alkene': {
        'smarts': '[C:1][C:2]([Cl,Br,I:3])[C:4]>>[C:1]=[C:2][C:4].[Cl,Br,I:3]',
        'description': 'Tertiary Alkyl Halide -> Alkene + Halide (SN1 via Carbocation)',
        'priority': 6
    },
    'SNAr_Phenoxide': {
        'smarts': '[c:1][Cl,Br:2].[O-:3][c:4]>>[c:1][O:3][c:4].[Cl,Br:2]',
        'description': 'Aryl Halide + Phenoxide -> Aryl Ether + Halide (SNAr)',
        'priority': 6
    },
    'SN2_Nitrile': {
        'smarts': '[C:1][Cl,Br,I:2].[C:3]#[N-:4]>>[C:1][C:3]#[N:4].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Cyanide Ion -> Nitrile + Halide',
        'priority': 7
    },

    # Additional Addition Reactions (10)
    'Hydroalkoxylation': {
        'smarts': '[C:1]=[C:2].[O:3][C:4]>>[C:1][C:2][O:3][C:4]',
        'description': 'Alkene + Alcohol -> Ether (Hydroalkoxylation)',
        'priority': 6
    },
    'Hydrothiolation': {
        'smarts': '[C:1]=[C:2].[S:3][C:4]>>[C:1][C:2][S:3][C:4]',
        'description': 'Alkene + Thiol -> Thioether (Hydrothiolation)',
        'priority': 6
    },
    'Alkyne_Hydroboration': {
        'smarts': '[C:1]#[C:2].[B:3][H:4]>>[C:1]=[C:2][B:3][H:4]',
        'description': 'Alkyne + Borane -> Vinylborane',
        'priority': 7
    },
    'Alkene_Hydrophosphination': {
        'smarts': '[C:1]=[C:2].[P:3][H:4]>>[C:1][C:2][P:3][H:4]',
        'description': 'Alkene + Phosphine -> Alkylphosphine',
        'priority': 6
    },
    'Alkyne_Hydrosilylation': {
        'smarts': '[C:1]#[C:2].[Si:3][H:4]>>[C:1]=[C:2][Si:3][H:4]',
        'description': 'Alkyne + Silane -> Vinylsilane',
        'priority': 7
    },
    'Alkene_Hydrocarboxylation': {
        'smarts': '[C:1]=[C:2].[C:3](=[O:4])[OH:5]>>[C:1][C:2][C:3](=[O:4])[OH:5]',
        'description': 'Alkene + Carboxylic Acid -> Alkylated Carboxylic Acid',
        'priority': 6
    },
    'Alkyne_Hydroamination': {
        'smarts': '[C:1]#[C:2].[N:3][H:4]>>[C:1]=[C:2][N:3][H:4]',
        'description': 'Alkyne + Amine -> Enamine',
        'priority': 7
    },
    'Dihydroxylation_Alkyne': {
        'smarts': '[C:1]#[C:2].[O:3][O:4]>>[C:1]([OH:3])[C:2]([OH:4])',
        'description': 'Alkyne + OsO4 -> 1,2-Diol',
        'priority': 7
    },
    'Halohydrin_Formation': {
        'smarts': '[C:1]=[C:2].[Cl:3][OH2:4]>>[C:1][Cl:3][C:2][OH:4]',
        'description': 'Alkene + Cl2/H2O -> Chlorohydrin',
        'priority': 7
    },
    'Alkene_Cyclopropanation': {
        'smarts': '[C:1]=[C:2].[C:3][H:4][I:5]>>[C:1]1[C:2][C:3]1.[H:4][I:5]',
        'description': 'Alkene + Carbenoid -> Cyclopropane (Simmons-Smith)',
        'priority': 7
    },

    # Additional Elimination Reactions (10)
    'E1_Alcohol': {
        'smarts': '[C:1][C:2]([OH:3])[C:4]>>[C:1]=[C:2][C:4].[OH2:3]',
        'description': 'Tertiary Alcohol -> Alkene + Water (E1 Dehydration)',
        'priority': 6
    },
    'E2_Alcohol': {
        'smarts': '[C:1][C:2]([OH:3])[H:4].[O-:5]>>[C:1]=[C:2].[OH2:3][O:5]',
        'description': 'Alcohol + Base -> Alkene + Water (E2)',
        'priority': 6
    },
    'Dehydrochlorination': {
        'smarts': '[C:1][C:2]([Cl:3])[H:4]>>[C:1]=[C:2].[Cl:3][H:4]',
        'description': 'Alkyl Chloride -> Alkene + HCl',
        'priority': 6
    },
    'Dehydrosulfurization': {
        'smarts': '[C:1][C:2]([S:3][H:4])[H:5]>>[C:1]=[C:2].[S:3][H:4][H:5]',
        'description': 'Thiol -> Alkene + H2S',
        'priority': 6
    },
    'E2_Amine': {
        'smarts': '[C:1][C:2]([N:3][C:4])[H:5]>>[C:1]=[C:2].[N:3][C:4][H:5]',
        'description': 'Tertiary Amine -> Alkene + Amine (Hofmann Elimination)',
        'priority': 6
    },
    'Dehydrobromination_Alkyne': {
        'smarts': '[C:1]([Br:2])[C:3]([H:4])[Br:5]>>[C:1]#[C:3].[Br:2][H:4][Br:5]',
        'description': 'Geminal Dihalide -> Alkyne + 2HX',
        'priority': 6
    },
    'Selenoxide_Elimination': {
        'smarts': '[C:1][C:2][Se+:3]([C:4])[O-:5]>>[C:1]=[C:2].[Se:3][C:4]',
        'description': 'Selenoxide -> Alkene + Selenide (Selenoxide Elimination)',
        'priority': 6
    },
    'Dehydroiodination': {
        'smarts': '[C:1][C:2]([I:3])[H:4]>>[C:1]=[C:2].[I:3][H:4]',
        'description': 'Alkyl Iodide -> Alkene + HI',
        'priority': 6
    },
    'Chugaev_Elimination': {
        'smarts': '[C:1][C:2]([O:3][C:4][S:5])[H:6]>>[C:1]=[C:2].[O:3][C:4][S:5][H:6]',
        'description': 'Alcohol Xanthate -> Alkene + COS + Thiol',
        'priority': 6
    },
    'Corey_Winter': {
        'smarts': '[C:1]([O:2])[C:3]([O:4])[S:5]>>[C:1]=[C:3].[O:2][S:5][O:4]',
        'description': '1,2-Diol -> Alkene + CO2 + S (Corey-Winter)',
        'priority': 6
    },

    # Additional Oxidation Reactions (10)
    'Allylic_Oxidation': {
        'smarts': '[C:1][C:2]=[C:3].[O:4]>>[C:1][C:2](=[O:4])[C:3]',
        'description': 'Alkene -> Allylic Ketone (SeO2 Oxidation)',
        'priority': 7
    },
    'Oxidative_Cleavage_Alkyne': {
        'smarts': '[C:1]#[C:2].[O:3][O:4]>>[C:1](=[O:3])[OH:4].[C:2](=[O:3])[OH:4]',
        'description': 'Alkyne -> Two Carboxylic Acids (Oxidative Cleavage)',
        'priority': 6
    },
    'Tosylate_Oxidation': {
        'smarts': '[C:1][O:2][S:3](=[O:4])(=[O:5])[c:6].[O:7]>>[C:1](=[O:2]).[S:3](=[O:4])(=[O:5])[c:6][O:7]',
        'description': 'Tosylate -> Carbonyl + Tosylate',
        'priority': 6
    },
    'SulfideToSulfone': {
        'smarts': '[C:1][S:2][C:3].[O:4][O:5]>>[C:1][S:2](=[O:4])(=[O:5])[C:3]',
        'description': 'Sulfide -> Sulfone (Double Oxidation)',
        'priority': 6
    },
    'AmineToImine': {
        'smarts': '[C:1][N:2][H:3].[O:4]>>[C:1]=[N:2].[OH2:4]',
        'description': 'Primary Amine -> Imine (Oxidation)',
        'priority': 7
    },
    'AlcoholToAcid_Direct': {
        'smarts': '[C:1][CH2:2][OH:3]>>[C:1][C:2](=[O])[OH:3]',
        'description': 'Primary Alcohol -> Carboxylic Acid (Direct Oxidation)',
        'priority': 7
    },
    'Dess_Martin_Oxidation': {
        'smarts': '[C:1][CH:2]([OH:3])[C:4]>>[C:1][C:2](=[O:3])[C:4]',
        'description': 'Secondary Alcohol -> Ketone (Dess-Martin)',
        'priority': 7
    },
    'Oppenauer_Oxidation': {
        'smarts': '[C:1][CH:2]([OH:3])[C:4].[C:5]=[O:6]>>[C:1][C:2](=[O:3])[C:4].[C:5][OH:6]',
        'description': 'Alcohol + Ketone -> Ketone + Alcohol (Oppenauer)',
        'priority': 7
    },
    'Oxidation_Thioether': {
        'smarts': '[C:1][S:2][C:3].[O:4]>>[C:1][S:2](=[O:4])[C:3]',
        'description': 'Thioether -> Sulfoxide',
        'priority': 6
    },
    'AmineToN_Oxide': {
        'smarts': '[C:1][N:2]([C:3])[C:4].[O:5]>>[C:1][N+:2]([C:3])[C:4][O-:5]',
        'description': 'Tertiary Amine -> Amine N-Oxide',
        'priority': 6
    },

    # Additional Reduction Reactions (10)
    'AlkyneToCisAlkene': {
        'smarts': '[C:1]#[C:2].[H:3][H:4]>>[C:1]=[C:2][H:3][H:4]',
        'description': 'Alkyne -> Cis-Alkene (Lindlar Reduction)',
        'priority': 7
    },
    'NitrileToAmine_Primary': {
        'smarts': '[C:1]#[N:2].[H:3][H:4]>>[C:1][CH2:2][NH2:2]',
        'description': 'Nitrile -> Primary Amine (Reduction)',
        'priority': 7
    },
    'ImineToSecondaryAmine': {
        'smarts': '[C:1]=[N:2][C:3].[H:4][H:5]>>[C:1][NH:2][C:3]',
        'description': 'Imine -> Secondary Amine',
        'priority': 7
    },
    'KetoneToMethylene': {
        'smarts': '[C:1][C:2](=[O:3])[C:4].[H:5][H:6]>>[C:1][CH2:2][C:4]',
        'description': 'Ketone -> Methylene (Clemmensen)',
        'priority': 6
    },
    'AldehydeToPrimaryAlcohol': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[H:5][H:6]>>[C:1][CH2:2][OH:3]',
        'description': 'Aldehyde -> Primary Alcohol (Reduction)',
        'priority': 7
    },
    'NitrosoToHydroxylamine': {
        'smarts': '[C:1][N:2]=[O:3].[H:4][H:5]>>[C:1][N:2][OH:3][H:4]',
        'description': 'Nitroso -> Hydroxylamine',
        'priority': 6
    },
    'AzideToPrimaryAmine': {
        'smarts': '[C:1][N:2]=[N+:3]=[N-:4].[H:5][H:6]>>[C:1][NH2:2].[N:3]#[N:4]',
        'description': 'Azide -> Primary Amine + Nitrogen',
        'priority': 7
    },
    'EsterToAldehyde': {
        'smarts': '[C:1][O:2][C:3](=[O:4]).[H:5][H:6]>>[C:1][OH:2].[C:3](=[O:4])[H:5]',
        'description': 'Ester -> Aldehyde (DIBAL-H Reduction)',
        'priority': 7
    },
    'AmideToAldehyde': {
        'smarts': '[C:1][C:2](=[O:3])[N:4].[H:5][H:6]>>[C:1][C:2](=[O:3])[H:5].[N:4][H:6]',
        'description': 'Amide -> Aldehyde (Reduction)',
        'priority': 7
    },
    'DisulfideToThiol_Reduction': {
        'smarts': '[C:1][S:2][S:3][C:4].[H:5][H:6]>>[C:1][S:2][H:5].[C:4][S:3][H:6]',
        'description': 'Disulfide -> Two Thiols (Reduction)',
        'priority': 6
    },

    # Additional Cross-Coupling Reactions (10)
    'Sonogashira_Alkyl': {
        'smarts': '[C:1][Br,I:2].[C:3]#[C:4][H:5]>>[C:1][C:3]#[C:4].[Br,I:2][H:5]',
        'description': 'Alkyl Halide + Terminal Alkyne -> Coupled Alkyne (Sonogashira)',
        'priority': 8
    },
    'Heck_Alkyl': {
        'smarts': '[C:1][Br,I:2].[C:3]=[C:4]>>[C:1][C:3]=[C:4].[Br,I:2]',
        'description': 'Alkyl Halide + Alkene -> Coupled Alkene (Heck)',
        'priority': 8
    },
    'Suzuki_Alkyl_Boronate': {
        'smarts': '[C:1][Br,I:2].[C:3][B:4]([O:5])[O:6]>>[C:1][C:3].[Br,I:2][B:4]([O:5])[O:6]',
        'description': 'Alkyl Halide + Boronate Ester -> Coupled Product',
        'priority': 7
    },
    'Stille_Alkyl': {
        'smarts': '[C:1][Br,I:2].[C:3][Sn:4]>>[C:1][C:3].[Br,I:2][Sn:4]',
        'description': 'Alkyl Halide + Organotin -> Coupled Product (Stille)',
        'priority': 8
    },
    'Negishi_Alkyl': {
        'smarts': '[C:1][Br,I:2].[C:3][Zn:4][Cl:5]>>[C:1][C:3].[Br,I:2][Zn:4][Cl:5]',
        'description': 'Alkyl Halide + Organozinc -> Coupled Product (Negishi)',
        'priority': 8
    },
    'Kumada_Alkyl': {
        'smarts': '[C:1][Br,I:2].[C:3][Mg:4][Cl:5]>>[C:1][C:3].[Br,I:2][Mg:4][Cl:5]',
        'description': 'Alkyl Halide + Grignard -> Coupled Product (Kumada)',
        'priority': 8
    },
    'Hiyama_Alkyl': {
        'smarts': '[C:1][Br,I:2].[C:3][Si:4]>>[C:1][C:3].[Br,I:2][Si:4]',
        'description': 'Alkyl Halide + Organosilane -> Coupled Product (Hiyama)',
        'priority': 8
    },
    'Buchwald_Hartwig_Secondary': {
        'smarts': '[c:1][Br,I:2].[N:3]([C:4])[H:5]>>[c:1][N:3][C:4].[Br,I:2][H:5]',
        'description': 'Aryl Halide + Secondary Amine -> Tertiary Aryl Amine (Buchwald-Hartwig)',
        'priority': 8
    },
    'Cross_Coupling_Vinyl': {
        'smarts': '[c:1][Br,I:2].[C:3]=[C:4]>>[c:1][C:3]=[C:4].[Br,I:2]',
        'description': 'Aryl Halide + Vinyl -> Aryl Alkene',
        'priority': 8
    },
    'Miyaura_Borylation_Alkyl': {
        'smarts': '[C:1][Br,I:2].[B:3]([O:4])[O:5]>>[C:1][B:3]([O:4])[O:5].[Br,I:2]',
        'description': 'Alkyl Halide -> Alkyl Boronate (Miyaura Borylation)',
        'priority': 8
    },

    # Additional Cycloaddition Reactions (10)
    '1,3-Dipolar_Nitrone': {
        'smarts': '[C:1]=[C:2].[C:3]=[N+:4][O-:5]>>[C:1]1[C:2][C:3][N+:4][O-:5]1',
        'description': 'Alkene + Nitrone -> Isoxazolidine',
        'priority': 7
    },
    'Diels_Alder_Alkyne': {
        'smarts': '[C:1]=[C:2][C:3]=[C:4].[C:5]#[C:6]>>[C:1]1[C:2][C:5][C:6][C:4][C:3]1',
        'description': 'Diene + Alkyne -> Cyclohexene (Diels-Alder)',
        'priority': 9
    },
    '[2+2]_Photocycloaddition': {
        'smarts': '[C:1]=[C:2].[C:3]=[C:4]>>[C:1]1[C:2][C:3][C:4]1',
        'description': 'Alkene + Alkene -> Cyclobutane ([2+2] Photocycloaddition)',
        'priority': 7
    },
    '1,3-Dipolar_Azomethine': {
        'smarts': '[C:1]=[C:2].[C:3]=[N:4][C:5]>>[C:1]1[C:2][C:3][N:4][C:5]1',
        'description': 'Alkene + Azomethine Ylide -> Pyrrolidine',
        'priority': 7
    },
    'Paterno_Buchi_Alkyne': {
        'smarts': '[C:1]#[C:2].[C:3]=[O:4]>>[C:1]1[C:2][O:4][C:3]1',
        'description': 'Alkyne + Carbonyl -> Oxetene (Paterno-Buchi)',
        'priority': 7
    },
    'Cycloaddition_Nitrile': {
        'smarts': '[C:1]=[C:2].[C:3]#[N:4]>>[C:1]1[C:2][C:3][N:4]1',
        'description': 'Alkene + Nitrile -> Four-Membered Heterocycle',
        'priority': 7
    },
    '[4+2]_Hetero_Diels_Alder_Nitroso': {
        'smarts': '[C:1]=[C:2][C:3]=[C:4].[N:5]=[O:6]>>[C:1]1[C:2][N:5][O:6][C:4][C:3]1',
        'description': 'Diene + Nitroso -> Oxazine (Hetero Diels-Alder)',
        'priority': 7
    },
    '1,5-Dipolar_Cycloaddition': {
        'smarts': '[C:1]=[C:2][C:3][N:4]=[N:5]>>[C:1]1[C:2][C:3][N:4][N:5]1',
        'description': '1,5-Dipole -> Five-Membered Heterocycle',
        'priority': 7
    },
    'Cycloaddition_Carbonyl_Ylide': {
        'smarts': '[C:1]=[C:2].[C:3][C:4][O:5]>>[C:1]1[C:2][C:3][C:4][O:5]1',
        'description': 'Alkene + Carbonyl Ylide -> Oxolane',
        'priority': 7
    },
    'Alkyne_Trimerization': {
        'smarts': '[C:1]#[C:2].[C:3]#[C:4].[C:5]#[C:6]>>[c:1]1[c:2][c:3][c:4][c:5][c:6]1',
        'description': 'Three Alkynes -> Benzene (Trimerization)',
        'priority': 8
    },

    # Additional Rearrangement Reactions (10)
    'Sigmatropic_1,3': {
        'smarts': '[C:1]=[C:2][C:3][C:4][H:5]>>[C:1][C:2]=[C:3][C:4][H:5]',
        'description': '1,3-Sigmatropic Rearrangement',
        'priority': 6
    },
    'Wolff_Rearrangement': {
        'smarts': '[C:1][C:2](=[O:3])[N:4]#[N:5]>>[C:1]=[C:2][O:3].[N:4]#[N:5]',
        'description': 'Diazo Ketone -> Ketene (Wolff)',
        'priority': 6
    },
    'Schmidt_Rearrangement': {
        'smarts': '[C:1](=[O:2])[N:3][N:4][H:5]>>[C:1][N:3]#[C:4].[OH2:2]',
        'description': 'Carboxylic Acid + Azide -> Isocyanate (Schmidt)',
        'priority': 6
    },
    'Baeyer_Villiger_Aldehyde': {
        'smarts': '[C:1][C:2](=[O:3])[H:4]>>[C:1][O:3][C:2][OH:4]',
        'description': 'Aldehyde -> Formate Ester (Baeyer-Villiger)',
        'priority': 7
    },
    'Claisen_Allyl_Ether': {
        'smarts': '[C:1]=[C:2][O:3][C:4]=[C:5]>>[C:1][C:2](=[O:3])[C:4]=[C:5]',
        'description': 'Allyl Vinyl Ether -> Carbonyl (Claisen)',
        'priority': 6
    },
    'Cope_Elimination_Oxide': {
        'smarts': '[C:1][C:2][N+:3]([C:4])[O-:5]>>[C:1]=[C:2].[N:3][C:4][OH:5]',
        'description': 'Amine Oxide -> Alkene + Hydroxylamine',
        'priority': 6
    },
    'Fischer_Indole': {
        'smarts': '[c:1][N:2][N:3][C:4](=[O:5])[C:6]>>[c:1]1[c:6][n:2][c:4][c:1]1',
        'description': 'Aryl Hydrazone -> Indole (Fischer)',
        'priority': 7
    },
    'Ireland_Claisen': {
        'smarts': '[C:1][C:2](=[O:3])[O:4][C:5]=[C:6]>>[C:1][C:2](=[O:3])[C:5]=[C:6][O:4][H]',
        'description': 'Silyl Enol Ether -> Rearranged Acid (Ireland-Claisen)',
        'priority': 6
    },
    'Benzilic_Acid': {
        'smarts': '[c:1][C:2](=[O:3])[C:4](=[O:5])[c:6]>>[c:1][C:2]([OH:3])[C:4](=[O:5])[c:6]',
        'description': 'Diketone -> Hydroxy Acid (Benzilic Acid Rearrangement)',
        'priority': 6
    },
    'Favorskii_Rearrangement': {
        'smarts': '[C:1][C:2](=[O:3])[C:4]([Cl,Br:5])[C:6]>>[C:1][C:2](=[O:3])[C:4][C:6].[Cl,Br:5]',
        'description': 'Alpha-Halo Ketone -> Carboxylic Acid (Favorskii)',
        'priority': 6
    },

    # Additional Heterocycle Formation Reactions (10)
    'Pyridine_Formation': {
        'smarts': '[C:1]=[O:2].[C:3]=[O:4].[N:5][H:6]>>[c:1]1[c:3][n:5][c:1][c:3][c:1]1',
        'description': 'Two Carbonyls + Ammonia -> Pyridine (Chichibabin)',
        'priority': 7
    },
    'Thiophene_Formation': {
        'smarts': '[C:1](=[O:2])[C:3][C:4](=[O:5]).[S:6]>>[c:1]1[c:3][c:4][s:6][c:1]1',
        'description': 'Diketone + Sulfur -> Thiophene',
        'priority': 7
    },
    'Furan_Formation': {
        'smarts': '[C:1](=[O:2])[C:3][C:4](=[O:5]).[O:6]>>[c:1]1[c:3][c:4][o:6][c:1]1',
        'description': 'Diketone + Oxygen -> Furan',
        'priority': 7
    },
    'Pyrazoline_Formation': {
        'smarts': '[C:1]=[C:2].[N:3][N:4][H:5]>>[C:1]1[C:2][N:3][N:4][H:5]1',
        'description': 'Alkene + Diazo Compound -> Pyrazoline',
        'priority': 7
    },
    'Oxazoline_Formation': {
        'smarts': '[C:1](=[O:2])[Cl:3].[N:4][C:5][OH:6]>>[c:1]1[o:2][c:5][n:4][c:1]1.[Cl:3][H:6]',
        'description': 'Acid Chloride + Amino Alcohol -> Oxazoline',
        'priority': 7
    },
    'Imidazoline_Formation': {
        'smarts': '[C:1](=[O:2])[Cl:3].[N:4][C:5][N:6][H:7]>>[c:1]1[n:4][c:5][n:6][c:1]1.[Cl:3][H:7]',
        'description': 'Acid Chloride + Diamine -> Imidazoline',
        'priority': 7
    },
    'Thiazoline_Formation': {
        'smarts': '[C:1](=[O:2])[Cl:3].[N:4][C:5][S:6][H:7]>>[c:1]1[s:6][c:5][n:4][c:1]1.[Cl:3][H:7]',
        'description': 'Acid Chloride + Amino Thiol -> Thiazoline',
        'priority': 7
    },
    'Indole_Formation_Skatole': {
        'smarts': '[c:1][N:2][H:3].[C:4](=[O:5])[C:6]>>[c:1]1[c:6][n:2][c:4][c:1]1',
        'description': 'Aniline + Ketone -> Indole (Skatole Synthesis)',
        'priority': 7
    },
    'Quinoline_Formation': {
        'smarts': '[c:1][N:2][H:3].[C:4](=[O:5])[C:6]=[C:7]>>[c:1]1[c:6][c:7][n:2][c:4][c:1]1',
        'description': 'Aniline + Enone -> Quinoline (Skraup Synthesis)',
        'priority': 7
    },
    'Isoquinoline_Formation': {
        'smarts': '[c:1][C:2][N:3][H:4].[C:5](=[O:6])>>[c:1]1[c:2][n:3][c:5][c:1]1.[OH2:6]',
        'description': 'Benzylamine + Carbonyl -> Isoquinoline (Pictet-Spengler)',
        'priority': 7
    },

    # Additional Protecting Group Chemistry (10)
    'Benzyl_Protection': {
        'smarts': '[C:1][OH:2].[C:3][Cl:4][c:5]>>[C:1][O:2][C:3][c:5].[Cl:4][H]',
        'description': 'Alcohol + Benzyl Chloride -> Benzyl Ether',
        'priority': 6
    },
    'TBS_Protection': {
        'smarts': '[C:1][OH:2].[Si:3]([C:4])([C:5])[Cl:6]>>[C:1][O:2][Si:3]([C:4])([C:5]).[Cl:6][H]',
        'description': 'Alcohol + TBS Chloride -> TBS-Protected Alcohol',
        'priority': 6
    },
    'Trityl_Protection': {
        'smarts': '[C:1][OH:2].[C:3]([c:4])([c:5])[Cl:6]>>[C:1][O:2][C:3]([c:4])([c:5]).[Cl:6][H]',
        'description': 'Alcohol + Trityl Chloride -> Trityl Ether',
        'priority': 6
    },
    'Boc_Protection_Secondary': {
        'smarts': '[N:1]([C:2])[H:3].[C:4](=[O:5])[O:6][C:7]>>[N:1]([C:2])[C:4](=[O:5])[O:6][C:7].[H:3]',
        'description': 'Secondary Amine + Boc Anhydride -> Boc-Protected Amine',
        'priority': 6
    },
    'Fmoc_Protection_Secondary': {
        'smarts': '[N:1]([C:2])[H:3].[C:4](=[O:5])[O:6][c:7]>>[N:1]([C:2])[C:4](=[O:5])[O:6][c:7].[H:3]',
        'description': 'Secondary Amine + Fmoc Chloride -> Fmoc-Protected Amine',
        'priority': 6
    },
    'Methyl_Ether_Protection': {
        'smarts': '[C:1][OH:2].[C:3][I:4]>>[C:1][O:2][C:3].[I:4][H]',
        'description': 'Alcohol + Methyl Iodide -> Methyl Ether',
        'priority': 6
    },
    'SEM_Protection': {
        'smarts': '[C:1][OH:2].[C:3][O:4][Si:5][Cl:6]>>[C:1][O:2][C:3][O:4][Si:5].[Cl:6][H]',
        'description': 'Alcohol + SEM Chloride -> SEM-Protected Alcohol',
        'priority': 6
    },
    'Acetal_Protection_Ketone': {
        'smarts': '[C:1][C:2](=[O:3])[C:4].[O:5][C:6][O:7]>>[C:1][C:2]([O:5][C:6][O:7])[C:4].[OH2:3]',
        'description': 'Ketone + Diol -> Cyclic Acetal',
        'priority': 6
    },
    'Cbz_Protection_Secondary': {
        'smarts': '[N:1]([C:2])[H:3].[C:4](=[O:5])[O:6][c:7]>>[N:1]([C:2])[C:4](=[O:5])[O:6][c:7].[H:3]',
        'description': 'Secondary Amine + Cbz Chloride -> Cbz-Protected Amine',
        'priority': 6
    },
    'PMB_Protection': {
        'smarts': '[C:1][OH:2].[C:3][Cl:4][c:5][O:6][C:7]>>[C:1][O:2][C:3][c:5][O:6][C:7].[Cl:4][H]',
        'description': 'Alcohol + PMB Chloride -> PMB Ether',
        'priority': 6
    },

    # Additional Condensation Reactions (10)
    'Aldol_Addition': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[C:5][C:6](=[O:7])[H:8]>>[C:1][C:2]([OH:3])[C:5][C:6](=[O:7])[H:8]',
        'description': 'Aldehyde + Aldehyde -> Beta-Hydroxy Aldehyde (Aldol Addition)',
        'priority': 7
    },
    'Cannizzaro_Reaction': {
        'smarts': '[c:1][C:2](=[O:3])[H:4].[c:5][C:6](=[O:7])[H:8]>>[c:1][C:2][OH:3].[c:5][C:6](=[O:7])[OH:8]',
        'description': 'Aromatic Aldehyde -> Alcohol + Carboxylic Acid (Cannizzaro)',
        'priority': 7
    },
    'Perkin_Condensation': {
        'smarts': '[c:1][C:2](=[O:3])[H:4].[C:5](=[O:6])[O:7][C:8]>>[c:1][C:2]=[C:5][C:6](=[O:7])[O:8].[OH2:3]',
        'description': 'Aromatic Aldehyde + Anhydride -> Cinnamic Acid (Perkin)',
        'priority': 7
    },
    'Knoevenagel_Nitrile': {
        'smarts': '[C:1]=[O:2].[C:3][C:4]#[N:5]>>[C:1]=[C:3][C:4]#[N:5].[OH2:2]',
        'description': 'Carbonyl + Active Methylene Nitrile -> Unsaturated Nitrile',
        'priority': 7
    },
    'Claisen_Schmidt': {
        'smarts': '[c:1][C:2](=[O:3])[H:4].[C:5][C:6](=[O:7])[H:8]>>[c:1][C:2]=[C:5][C:6](=[O:7]).[OH2:3]',
        'description': 'Aromatic Aldehyde + Ketone -> Enone (Claisen-Schmidt)',
        'priority': 7
    },
    'Stobbe_Condensation': {
        'smarts': '[C:1]=[O:2].[C:3]([C:4](=[O:5])[O:6])[C:7](=[O:8])[O:9]>>[C:1]=[C:3][C:4](=[O:5])[O:6].[C:7](=[O:8])[O:9]',
        'description': 'Carbonyl + Succinic Ester -> Unsaturated Acid (Stobbe)',
        'priority': 7
    },
    'Tishchenko_Reaction': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[C:5][C:6](=[O:7])[H:8]>>[C:1][C:2](=[O:3])[O:7][C:6][C:5]',
        'description': 'Two Aldehydes -> Ester (Tishchenko)',
        'priority': 7
    },
    'Enamine_Formation': {
        'smarts': '[C:1][C:2](=[O:3]).[N:4][H:5]>>[C:1][C:2]=[N:4].[OH2:3]',
        'description': 'Ketone + Amine -> Enamine + Water',
        'priority': 7
    },
    'Acetoacetic_Ester_Condensation': {
        'smarts': '[C:1][C:2](=[O:3])[O:4][C:5].[C:6][C:7](=[O:8])[O:9][C:10]>>[C:1][C:2](=[O:3])[C:6][C:7](=[O:8])[OH:9].[C:5][O:4][C:10]',
        'description': 'Two Esters -> Beta-Ketoester (Acetoacetic Ester Condensation)',
        'priority': 7
    },
    'Knoevenagel_Ester': {
        'smarts': '[C:1]=[O:2].[C:3][C:4](=[O:5])[O:6][C:7]>>[C:1]=[C:3][C:4](=[O:5])[O:6][C:7].[OH2:2]',
        'description': 'Carbonyl + Active Methylene Ester -> Unsaturated Ester',
        'priority': 7
    },

    # Additional Multi-Component Reactions (10)
    'Ugi_Three_Component': {
        'smarts': '[C:1]=[O:2].[N:3][H:4].[C:5][N:6]#[C:7]>>[C:1]([N:3][C:5])[N:6][C:7][OH:2]',
        'description': 'Aldehyde + Amine + Isocyanide -> Amide (Ugi 3-Component)',
        'priority': 8
    },
    'Passerini_Aldehyde': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[C:5](=[O:6])[OH:7].[C:8][N:9]#[C:10]>>[C:1][C:2]([O:6][C:5])[O:7][C:8][C:10][N:9][H:4]',
        'description': 'Aldehyde + Acid + Isocyanide -> Ester Amide (Passerini)',
        'priority': 7
    },
    'Biginelli_Urea': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[C:5][C:6](=[O:7])[C:8].[N:9][C:10](=[O:11])[N:12]>>[c:1]1[n:9][c:6][c:8][n:12][c:1]1.[OH2:3]',
        'description': 'Aldehyde + Beta-Ketoester + Urea -> Dihydropyrimidinone (Biginelli)',
        'priority': 7
    },
    'Hantzsch_Pyrrole': {
        'smarts': '[C:1](=[O:2])[C:3][C:4](=[O:5]).[N:6][H:7][C:8]>>[c:1]1[c:3][c:4][n:6][c:8]1',
        'description': 'Diketone + Amine -> Substituted Pyrrole (Hantzsch)',
        'priority': 7
    },
    'Gewald_Thiazole': {
        'smarts': '[C:1][C:2](=[O:3])[C:4].[S:5].[C:6]#[N:7]>>[c:1]1[c:2][c:4][s:5][c:6]1.[N:7][H]',
        'description': 'Ketone + Sulfur + Nitrile -> Thiazole (Gewald)',
        'priority': 7
    },
    'Povarov_Reaction': {
        'smarts': '[c:1][N:2][H:3].[C:4]=[C:5].[C:6]=[O:7]>>[c:1]1[c:4][c:5][n:2][c:6][c:1]1.[OH2:7]',
        'description': 'Aniline + Alkene + Aldehyde -> Quinoline (Povarov)',
        'priority': 7
    },
    'Mannich_Imine': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[N:5][C:6].[C:7]=[O:8]>>[C:1][C:2](=[O:3])[C:6][N:5][C:7][OH:8]',
        'description': 'Aldehyde + Imine + Carbonyl -> Beta-Amino Carbonyl (Mannich)',
        'priority': 7
    },
    'Strecker_Nitrile': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[N:5][H:6][H:7].[C:8]#[N:9]>>[C:1][C:2]([N:5][H:7])[C:8]#[N:9].[OH2:3]',
        'description': 'Aldehyde + Ammonia + Cyanide -> Amino Nitrile (Strecker)',
        'priority': 7
    },
    'Paal_Knorr_Furan': {
        'smarts': '[C:1](=[O:2])[C:3][C:4](=[O:5]).[O:6][H:7]>>[c:1]1[c:3][c:4][o:6][c:1]1.[OH2:5]',
        'description': 'Diketone + Alcohol -> Furan (Paal-Knorr)',
        'priority': 7
    },
    'Hantzsch_Thiazole': {
        'smarts': '[C:1](=[O:2])[Cl:3].[C:4][C:5](=[S:6])[N:7][H:8]>>[c:1]1[c:4][s:6][c:5][n:7]1.[Cl:3][H:8]',
        'description': 'Alpha-Halo Ketone + Thioamide -> Thiazole (Hantzsch)',
        'priority': 7
    },

    # Additional Specialized Reactions (10)
    'Wacker_Oxidation_Alkyne': {
        'smarts': '[C:1]#[C:2][C:3].[O:4]>>[C:1][C:2](=[O:4])[C:3]',
        'description': 'Terminal Alkyne -> Methyl Ketone (Wacker)',
        'priority': 7
    },
    'Ring_Closing_Metathesis': {
        'smarts': '[C:1]=[C:2][C:3].[C:4][C:5]=[C:6]>>[C:1]1[C:2][C:5][C:6]1.[C:3]=[C:4]',
        'description': 'Diene -> Cyclic Alkene + Ethylene (Ring-Closing Metathesis)',
        'priority': 8
    },
    'Pauson_Khand_Alkyne': {
        'smarts': '[C:1]#[C:2].[C:3]=[C:4].[C:5]=[O:6]>>[C:1]1[C:2][C:3][C:4][C:5](=[O:6])1',
        'description': 'Alkyne + Alkene + CO -> Cyclopentenone (Pauson-Khand)',
        'priority': 7
    },
    'Hydroacylation': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[C:5]=[C:6]>>[C:1][C:2](=[O:3])[C:5][C:6]',
        'description': 'Aldehyde + Alkene -> Ketone (Hydroacylation)',
        'priority': 7
    },
    'C_H_Insertion': {
        'smarts': '[C:1][H:2].[C:3][N:4]#[N:5]>>[C:1][C:3].[N:4]#[N:5][H:2]',
        'description': 'Alkane + Diazo -> Alkyl Insertion (C-H Insertion)',
        'priority': 7
    },
    'Metathesis_Cross': {
        'smarts': '[C:1]=[C:2][C:3].[C:4]=[C:5][C:6]>>[C:1]=[C:4].[C:2][C:3]=[C:5][C:6]',
        'description': 'Two Alkenes -> Crossed Alkenes (Cross Metathesis)',
        'priority': 8
    },
    'Hydroarylation': {
        'smarts': '[c:1][H:2].[C:3]=[C:4]>>[c:1][C:3][C:4]',
        'description': 'Aromatic C-H + Alkene -> Alkyl Aromatic (Hydroarylation)',
        'priority': 7
    },
    'Click_Reaction_Thiol': {
        'smarts': '[C:1]=[C:2].[S:3][H:4]>>[C:1][C:2][S:3][H:4]',
        'description': 'Alkene + Thiol -> Thioether (Thiol-Ene Click)',
        'priority': 7
    },
    'Oxidative_Coupling': {
        'smarts': '[c:1][H:2].[c:3][H:4]>>[c:1][c:3].[H:2][H:4]',
        'description': 'Two Aryl C-H -> Biaryl (Oxidative Coupling)',
        'priority': 7
    },
    'Hydroboration_Oxidation': {
        'smarts': '[C:1]=[C:2].[B:3][H:4].[O:5]>>[C:1][C:2][OH:5].[B:3][H:4]',
        'description': 'Alkene + Borane + Oxidation -> Alcohol',
        'priority': 7
    }
}

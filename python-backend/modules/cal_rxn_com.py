

def calculate_reaction_compatibility(groups1, groups2, reaction_type):
    """Calculate compatibility score between reactants for a specific reaction."""
    compatibility_rules = {
        'Esterification': {
            'required': [(['alcohol'], ['carboxylic_acid']), (['carboxylic_acid'], ['alcohol'])],
            'score': 0.9
        },
        'AmideFormation': {
            'required': [(['amine', 'primary_amine', 'secondary_amine'], ['carboxylic_acid']),
                         (['carboxylic_acid'], ['amine', 'primary_amine', 'secondary_amine'])],
            'score': 0.9
        },
        'Hydrolysis': {
            'required': [(['ester'], []), ([], ['ester'])],
            'score': 0.7
        },
        'AmideHydrolysis': {
            'required': [(['amide'], []), ([], ['amide'])],
            'score': 0.7
        },
        'SN2_Alcohol': {
            'required': [(['alkyl_halide'], ['alcohol']), (['alcohol'], ['alkyl_halide'])],
            'score': 0.7
        },
        'SN2_Amine': {
            'required': [(['alkyl_halide'], ['amine']), (['amine'], ['alkyl_halide'])],
            'score': 0.7
        },
        'SN2_Thiol': {
            'required': [(['alkyl_halide'], ['thiol']), (['thiol'], ['alkyl_halide'])],
            'score': 0.7
        },
        'SN2_Alkoxide': {
            'required': [(['alkyl_halide'], ['ether']), (['ether'], ['alkyl_halide'])],
            'score': 0.7
        },
        'HalogenAddition': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.8
        },
        'HydrogenHalideAddition': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.8
        },
        'WaterAddition': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.7
        },
        'HydrogenAddition': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.8
        },
        'Dehydration': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.6
        },
        'Dehydrohalogenation': {
            'required': [(['alkyl_halide'], []), ([], ['alkyl_halide'])],
            'score': 0.6
        },
        'AlcoholToAldehyde': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.7
        },
        'AlcoholToKetone': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.7
        },
        'AldehydeToAcid': {
            'required': [(['aldehyde'], []), ([], ['aldehyde'])],
            'score': 0.7
        },
        'CarbonylToAlcohol': {
            'required': [(['aldehyde'], []), ([], ['aldehyde'])],
            'score': 0.7
        },
        'KetoneToAlcohol': {
            'required': [(['ketone'], []), ([], ['ketone'])],
            'score': 0.7
        },
        'CarboxylicAcidToAlcohol': {
            'required': [(['carboxylic_acid'], []), ([], ['carboxylic_acid'])],
            'score': 0.6
        },
        'Nitration': {
            'required': [(['aromatic'], []), ([], ['aromatic'])],
            'score': 0.7
        },
        'Halogenation_Aromatic': {
            'required': [(['aromatic'], []), ([], ['aromatic'])],
            'score': 0.7
        },
        'Friedel_Crafts_Alkylation': {
            'required': [(['aromatic'], ['alkyl_halide']), (['alkyl_halide'], ['aromatic'])],
            'score': 0.6
        },
        'Friedel_Crafts_Acylation': {
            'required': [(['aromatic'], []), ([], ['aromatic'])],
            'score': 0.6
        },
        'Suzuki_Coupling': {
            'required': [(['aryl_halide'], ['boronic_acid', 'boronate_ester']),
                         (['boronic_acid', 'boronate_ester'], ['aryl_halide'])],
            'score': 0.95
        },
        'Heck_Reaction': {
            'required': [(['aryl_halide'], ['alkene']), (['alkene'], ['aryl_halide'])],
            'score': 0.9
        },
        'Sonogashira_Coupling': {
            'required': [(['aryl_halide'], ['alkyne']), (['alkyne'], ['aryl_halide'])],
            'score': 0.8
        },
        'Diels_Alder': {
            'required': [(['alkene'], ['alkene'])],
            'score': 0.85
        },
        'Cyclopropanation': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.7
        },
        'Claisen_Rearrangement': {
            'required': [(['ether'], []), ([], ['ether'])],
            'score': 0.6
        },
        'Beckmann_Rearrangement': {
            'required': [(['oxime'], []), ([], ['oxime'])],
            'score': 0.6
        },
        'ThiocarbonylHydrolysis': {
            'required': [(['thiocarbonyl'], []), ([], ['thiocarbonyl'])],
            'score': 0.6
        },
        'ThioureaFormation': {
            'required': [(['amine'], []), ([], ['amine'])],
            'score': 0.7
        },
        'SecondaryAmineThiourea': {
            'required': [(['secondary_amine'], []), ([], ['secondary_amine'])],
            'score': 0.7
        },
        'Dehalogenation': {
            'required': [(['aryl_halide'], []), ([], ['aryl_halide'])],
            'score': 0.5
        },
        'Aldol_Condensation': {
            'required': [(['aldehyde'], ['aldehyde'])],
            'score': 0.7
        },
        'Mannich_Reaction': {
            'required': [(['aldehyde'], ['amine'])],
            'score': 0.7
        },
        'Knoevenagel_Condensation': {
            'required': [(['aldehyde'], ['carboxylic_acid'])],
            'score': 0.7
        },
        'Grignard_Carbonyl': {
            'required': [(['organomagnesium'], ['carbonyl', 'aldehyde', 'ketone'])],
            'score': 0.9
        },
        'Grignard_CO2': {
            'required': [(['organomagnesium'], []), ([], ['organomagnesium'])],
            'score': 0.8
        },
        'Silyl_Protection': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.6
        },
        'Acetal_Formation': {
            'required': [(['carbonyl'], ['alcohol'])],
            'score': 0.6
        },
        'Pyrrole_Formation': {
            'required': [(['ketone'], ['amine'])],
            'score': 0.7
        },
        'Imidazole_Formation': {
            'required': [(['ketone'], ['amine', 'aldehyde'])],
            'score': 0.7
        },
        'Ugi_Reaction': {
            'required': [(['carbonyl'], ['amine', 'carboxylic_acid', 'isocyanate'])],
            'score': 0.8
        },
        'Biginelli_Reaction': {
            'required': [(['aldehyde'], ['urea', 'ketone'])],
            'score': 0.7
        },
        'SN2_Cyanide': {
            'required': [(['alkyl_halide'], ['nitrile']), (['nitrile'], ['alkyl_halide'])],
            'score': 0.7
        },
        'SN2_Azide': {
            'required': [(['alkyl_halide'], ['azide']), (['azide'], ['alkyl_halide'])],
            'score': 0.7
        },
        'SNAr_Fluoride': {
            'required': [(['aryl_halide'], ['amine']), (['amine'], ['aryl_halide'])],
            'score': 0.6
        },
        'SNAr_Alkoxide': {
            'required': [(['aryl_halide'], ['ether']), (['ether'], ['aryl_halide'])],
            'score': 0.6
        },
        'Hydroboration': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.7
        },
        'Epoxidation': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.7
        },
        'Dihydroxylation': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.7
        },
        'Ozonolysis': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.6
        },
        'E2_Elimination': {
            'required': [(['alkyl_halide'], []), ([], ['alkyl_halide'])],
            'score': 0.6
        },
        'Dehalogenation_Elimination': {
            'required': [(['alkyl_halide'], []), ([], ['alkyl_halide'])],
            'score': 0.6
        },
        'AlcoholToAcid': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.7
        },
        'AlkeneToEpoxide': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.7
        },
        'SulfideToSulfoxide': {
            'required': [(['thioether'], []), ([], ['thioether'])],
            'score': 0.6
        },
        'SulfoxideToSulfone': {
            'required': [(['sulfoxide'], []), ([], ['sulfoxide'])],
            'score': 0.6
        },
        'NitroToAmine': {
            'required': [(['nitro'], []), ([], ['nitro'])],
            'score': 0.7
        },
        'NitrileToAmine': {
            'required': [(['nitrile'], []), ([], ['nitrile'])],
            'score': 0.7
        },
        'ImineToAmine': {
            'required': [(['imine'], []), ([], ['imine'])],
            'score': 0.7
        },
        'Negishi_Coupling': {
            'required': [(['aryl_halide'], ['organozinc']), (['organozinc'], ['aryl_halide'])],
            'score': 0.8
        },
        'Stille_Coupling': {
            'required': [(['aryl_halide'], ['organotin']), (['organotin'], ['aryl_halide'])],
            'score': 0.8
        },
        'Kumada_Coupling': {
            'required': [(['aryl_halide'], ['organomagnesium']), (['organomagnesium'], ['aryl_halide'])],
            'score': 0.8
        },
        '[2+2]_Cycloaddition': {
            'required': [(['alkene'], ['alkene'])],
            'score': 0.7
        },
        '1,3-Dipolar_Cycloaddition': {
            'required': [(['alkene'], ['azide']), (['azide'], ['alkene'])],
            'score': 0.7
        },
        'Pinacol_Rearrangement': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.6
        },
        'Wagner_Meerwein': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.6
        },
        'Pyrazole_Formation': {
            'required': [(['ketone'], ['hydrazine']), (['hydrazine'], ['ketone'])],
            'score': 0.7
        },
        'Triazole_Formation': {
            'required': [(['alkyne'], ['azide']), (['azide'], ['alkyne'])],
            'score': 0.8
        },
        'Boc_Protection': {
            'required': [(['amine'], []), ([], ['amine'])],
            'score': 0.6
        },
        'Cbz_Protection': {
            'required': [(['amine'], []), ([], ['amine'])],
            'score': 0.6
        },
        'Claisen_Condensation': {
            'required': [(['ester'], ['ester'])],
            'score': 0.7
        },
        'Perkin_Reaction': {
            'required': [(['aromatic'], ['carboxylic_acid'])],
            'score': 0.7
        },
        'Passerini_Reaction': {
            'required': [(['carbonyl'], ['carboxylic_acid', 'isocyanate'])],
            'score': 0.7
        },
        'Hantzsch_Pyridine': {
            'required': [(['aldehyde'], ['ketone', 'amine'])],
            'score': 0.7
        },
        'SN1_Alcohol': {
            'required': [(['alkyl_halide'], ['alcohol']), (['alcohol'], ['alkyl_halide'])],
            'score': 0.6
        },
        'SN1_Ether': {
            'required': [(['alkyl_halide'], ['alcohol']), (['alcohol'], ['alkyl_halide'])],
            'score': 0.6
        },
        'SNAr_Nitro': {
            'required': [(['aryl_halide'], ['amine']), (['amine'], ['aryl_halide'])],
            'score': 0.6
        },
        'SN2_Phosphine': {
            'required': [(['alkyl_halide'], ['phosphine']), (['phosphine'], ['alkyl_halide'])],
            'score': 0.6
        },
        'SN2_Acetylide': {
            'required': [(['alkyl_halide'], ['acetylide']), (['acetylide'], ['alkyl_halide'])],
            'score': 0.7
        },
        'SN2_Hydrosulfide': {
            'required': [(['alkyl_halide'], ['thiol']), (['thiol'], ['alkyl_halide'])],
            'score': 0.7
        },
        'SNAr_Thiolate': {
            'required': [(['aryl_halide'], ['thiol']), (['thiol'], ['aryl_halide'])],
            'score': 0.6
        },
        'SN2_Enolate': {
            'required': [(['alkyl_halide'], ['enolate']), (['enolate'], ['alkyl_halide'])],
            'score': 0.7
        },
        'SN2_Carboxylate': {
            'required': [(['alkyl_halide'], ['carboxylic_acid']), (['carboxylic_acid'], ['alkyl_halide'])],
            'score': 0.7
        },
        'SN2_Hydroxide': {
            'required': [(['alkyl_halide'], ['alcohol']), (['alcohol'], ['alkyl_halide'])],
            'score': 0.7
        },
        'Michael_Addition': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.7
        },
        'Hydroformylation': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.7
        },
        'Hydrosilylation': {
            'required': [(['alkene'], ['silane']), (['silane'], ['alkene'])],
            'score': 0.6
        },
        'Oxymercuration': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.6
        },
        'Hydroamination': {
            'required': [(['alkene'], ['amine']), (['amine'], ['alkene'])],
            'score': 0.6
        },
        'Hydrocyanation': {
            'required': [(['alkene'], ['nitrile']), (['nitrile'], ['alkene'])],
            'score': 0.6
        },
        'Alkyne_Hydration': {
            'required': [(['alkyne'], []), ([], ['alkyne'])],
            'score': 0.7
        },
        'Bromohydrin_Formation': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.7
        },
        'Carbene_Addition': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.7
        },
        'Alkyne_Halogenation': {
            'required': [(['alkyne'], []), ([], ['alkyne'])],
            'score': 0.7
        },
        'E1_Elimination': {
            'required': [(['alkyl_halide'], []), ([], ['alkyl_halide'])],
            'score': 0.6
        },
        'Dehydrohalogenation_Alkyne': {
            'required': [(['alkyl_halide'], []), ([], ['alkyl_halide'])],
            'score': 0.6
        },
        'Dehydration_Ketone': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.6
        },
        'Elimination_AmineOxide': {
            'required': [(['tertiary_amine'], []), ([], ['tertiary_amine'])],
            'score': 0.6
        },
        'Dehydrofluorination': {
            'required': [(['alkyl_halide'], []), ([], ['alkyl_halide'])],
            'score': 0.6
        },
        'Baeyer_Villiger': {
            'required': [(['ketone'], []), ([], ['ketone'])],
            'score': 0.7
        },
        'Swern_Oxidation': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.7
        },
        'Jones_Oxidation': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.7
        },
        'Alkene_Cleavage': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.6
        },
        'ThiolToDisulfide': {
            'required': [(['thiol'], ['thiol'])],
            'score': 0.6
        },
        'AmineToNitroso': {
            'required': [(['amine'], []), ([], ['amine'])],
            'score': 0.6
        },
        'AlkyneToDiketone': {
            'required': [(['alkyne'], []), ([], ['alkyne'])],
            'score': 0.6
        },
        'AlcoholToEster': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.7
        },
        'AldehydeToAcetal': {
            'required': [(['aldehyde'], ['alcohol']), (['alcohol'], ['aldehyde'])],
            'score': 0.6
        },
        'KetoneToAcetal': {
            'required': [(['ketone'], ['alcohol']), (['alcohol'], ['ketone'])],
            'score': 0.6
        },
        'AzideToAmine': {
            'required': [(['azide'], []), ([], ['azide'])],
            'score': 0.7
        },
        'AlkyneToAlkene': {
            'required': [(['alkyne'], []), ([], ['alkyne'])],
            'score': 0.7
        },
        'EsterToAlcohol': {
            'required': [(['ester'], []), ([], ['ester'])],
            'score': 0.7
        },
        'AmideToAmine': {
            'required': [(['amide'], []), ([], ['amide'])],
            'score': 0.7
        },
        'NitrosoToAmine': {
            'required': [(['nitro'], []), ([], ['nitro'])],
            'score': 0.6
        },
        'AlkyneToAlkane': {
            'required': [(['alkyne'], []), ([], ['alkyne'])],
            'score': 0.7
        },
        'EpoxideToAlcohol': {
            'required': [(['epoxide'], []), ([], ['epoxide'])],
            'score': 0.6
        },
        'DisulfideToThiol': {
            'required': [(['disulfide'], []), ([], ['disulfide'])],
            'score': 0.6
        },
        'NitrileToAldehyde': {
            'required': [(['nitrile'], []), ([], ['nitrile'])],
            'score': 0.6
        },
        'KetoneToAlkane': {
            'required': [(['ketone'], []), ([], ['ketone'])],
            'score': 0.6
        },
        'Buchwald_Hartwig': {
            'required': [(['aryl_halide'], ['amine']), (['amine'], ['aryl_halide'])],
            'score': 0.8
        },
        'Hiyama_Coupling': {
            'required': [(['aryl_halide'], ['silane']), (['silane'], ['aryl_halide'])],
            'score': 0.8
        },
        'Miyaura_Borylation': {
            'required': [(['aryl_halide'], []), ([], ['aryl_halide'])],
            'score': 0.8
        },
        'Cross_Coupling_Alkyne': {
            'required': [(['aryl_halide'], ['alkyne']), (['alkyne'], ['aryl_halide'])],
            'score': 0.8
        },
        'Suzuki_Coupling_Alkyl': {
            'required': [(['alkyl_halide'], ['boronic_acid']), (['boronic_acid'], ['alkyl_halide'])],
            'score': 0.7
        },
        '1,3-Dipolar_NitrileOxide': {
            'required': [(['alkene'], ['nitrile_oxide']), (['nitrile_oxide'], ['alkene'])],
            'score': 0.7
        },
        '[4+2]_Hetero_Diels_Alder': {
            'required': [(['alkene'], ['carbonyl']), (['carbonyl'], ['alkene'])],
            'score': 0.7
        },
        'Paterno_Buchi': {
            'required': [(['alkene'], ['carbonyl']), (['carbonyl'], ['alkene'])],
            'score': 0.7
        },
        'Azide_Alkyne_Cycloaddition': {
            'required': [(['alkyne'], ['azide']), (['azide'], ['alkyne'])],
            'score': 0.8
        },
        '[3+2]_Cycloaddition': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.7
        },
        'Curtius_Rearrangement': {
            'required': [(['azide'], []), ([], ['azide'])],
            'score': 0.6
        },
        'Hofmann_Rearrangement': {
            'required': [(['amide'], []), ([], ['amide'])],
            'score': 0.6
        },
        'Cope_Rearrangement': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.6
        },
        'Overman_Rearrangement': {
            'required': [(['amide'], []), ([], ['amide'])],
            'score': 0.6
        },
        'Stevens_Rearrangement': {
            'required': [(['tertiary_amine'], []), ([], ['tertiary_amine'])],
            'score': 0.6
        },
        'Oxazole_Formation': {
            'required': [(['amide'], []), ([], ['amide'])],
            'score': 0.7
        },
        'Thiazole_Formation': {
            'required': [(['thioether'], []), ([], ['thioether'])],
            'score': 0.7
        },
        'Pyrimidine_Formation': {
            'required': [(['ketone'], ['urea']), (['urea'], ['ketone'])],
            'score': 0.7
        },
        'Isoxazole_Formation': {
            'required': [(['alkene'], ['nitrile_oxide']), (['nitrile_oxide'], ['alkene'])],
            'score': 0.7
        },
        'Indazole_Formation': {
            'required': [(['ketone'], ['hydrazine']), (['hydrazine'], ['ketone'])],
            'score': 0.7
        },
        'Fmoc_Protection': {
            'required': [(['amine'], []), ([], ['amine'])],
            'score': 0.6
        },
        'Tosylate_Formation': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.6
        },
        'THP_Protection': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.6
        },
        'MOM_Protection': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.6
        },
        'Acetyl_Protection': {
            'required': [(['alcohol'], []), ([], ['alcohol'])],
            'score': 0.6
        },
        'Dieckmann_Condensation': {
            'required': [(['ester'], ['ester'])],
            'score': 0.7
        },
        'Benzoin_Condensation': {
            'required': [(['aldehyde'], ['aldehyde'])],
            'score': 0.7
        },
        'Knoevenagel_Malonic': {
            'required': [(['aldehyde'], ['carboxylic_acid'])],
            'score': 0.7
        },
        'Robinson_Annulation': {
            'required': [(['ketone'], ['alkene'])],
            'score': 0.7
        },
        'Condensation_Imine': {
            'required': [(['carbonyl'], ['amine']), (['amine'], ['carbonyl'])],
            'score': 0.7
        },
        'Strecker_Synthesis': {
            'required': [(['aldehyde'], ['amine', 'nitrile'])],
            'score': 0.7
        },
        'Gewald_Reaction': {
            'required': [(['ketone'], []), ([], ['ketone'])],
            'score': 0.7
        },
        'Paal_Knorr_Pyrrole': {
            'required': [(['ketone'], ['amine']), (['amine'], ['ketone'])],
            'score': 0.7
        },
        'Hantzsch_Dihydropyridine': {
            'required': [(['aldehyde'], ['ketone', 'amine'])],
            'score': 0.7
        },
        'Biginelli_Pyrimidine': {
            'required': [(['aldehyde'], ['urea', 'ketone'])],
            'score': 0.7
        },
        'Olefin_Metathesis': {
            'required': [(['alkene'], ['alkene'])],
            'score': 0.8
        },
        'CH_Activation': {
            'required': [(['aromatic'], ['alkyl_halide']), (['alkyl_halide'], ['aromatic'])],
            'score': 0.7
        },
        'Wacker_Oxidation': {
            'required': [(['alkene'], []), ([], ['alkene'])],
            'score': 0.7
        },
        'Pauson_Khand': {
            'required': [(['alkyne'], ['alkene', 'carbonyl'])],
            'score': 0.7
        },
        'Grubbs_Metathesis': {
            'required': [(['alkene'], ['alkene'])],
            'score': 0.8
        },
           # Additional Nucleophilic Substitution Reactions (10)
    'SN2_Alkylthiolate': {
        'required': [(['alkyl_halide'], ['thioether']), (['thioether'], ['alkyl_halide'])],
        'score': 0.7
    },
    'SN2_Alkynide': {
        'required': [(['alkyl_halide'], ['acetylide']), (['acetylide'], ['alkyl_halide'])],
        'score': 0.7
    },
    'SNAr_Amine': {
        'required': [(['aryl_halide'], ['amine', 'primary_amine', 'secondary_amine']), 
                     (['amine', 'primary_amine', 'secondary_amine'], ['aryl_halide'])],
        'score': 0.6
    },
    'SN2_Phosphite': {
        'required': [(['alkyl_halide'], ['phosphine']), (['phosphine'], ['alkyl_halide'])],
        'score': 0.6
    },
    'SN2_Selenol': {
        'required': [(['alkyl_halide'], ['thiol']), (['thiol'], ['alkyl_halide'])],
        'score': 0.7
    },
    'SN2_Alkylamine_Secondary': {
        'required': [(['alkyl_halide'], ['secondary_amine']), (['secondary_amine'], ['alkyl_halide'])],
        'score': 0.7
    },
    'SN2_Carboxylate_Ester': {
        'required': [(['alkyl_halide'], ['carboxylic_acid']), (['carboxylic_acid'], ['alkyl_halide'])],
        'score': 0.7
    },
    'SN1_Alkene': {
        'required': [(['alkyl_halide'], []), ([], ['alkyl_halide'])],
        'score': 0.6
    },
    'SNAr_Phenoxide': {
        'required': [(['aryl_halide'], ['phenol']), (['phenol'], ['aryl_halide'])],
        'score': 0.6
    },
    'SN2_Nitrile': {
        'required': [(['alkyl_halide'], ['nitrile']), (['nitrile'], ['alkyl_halide'])],
        'score': 0.7
    },

    # Additional Addition Reactions (10)
    'Hydroalkoxylation': {
        'required': [(['alkene'], ['alcohol']), (['alcohol'], ['alkene'])],
        'score': 0.6
    },
    'Hydrothiolation': {
        'required': [(['alkene'], ['thiol']), (['thiol'], ['alkene'])],
        'score': 0.6
    },
    'Alkyne_Hydroboration': {
        'required': [(['alkyne'], []), ([], ['alkyne'])],
        'score': 0.7
    },
    'Alkene_Hydrophosphination': {
        'required': [(['alkene'], ['phosphine']), (['phosphine'], ['alkene'])],
        'score': 0.6
    },
    'Alkyne_Hydrosilylation': {
        'required': [(['alkyne'], ['silane']), (['silane'], ['alkyne'])],
        'score': 0.7
    },
    'Alkene_Hydrocarboxylation': {
        'required': [(['alkene'], ['carboxylic_acid']), (['carboxylic_acid'], ['alkene'])],
        'score': 0.6
    },
    'Alkyne_Hydroamination': {
        'required': [(['alkyne'], ['amine', 'primary_amine', 'secondary_amine']), 
                     (['amine', 'primary_amine', 'secondary_amine'], ['alkyne'])],
        'score': 0.7
    },
    'Dihydroxylation_Alkyne': {
        'required': [(['alkyne'], []), ([], ['alkyne'])],
        'score': 0.7
    },
    'Halohydrin_Formation': {
        'required': [(['alkene'], []), ([], ['alkene'])],
        'score': 0.7
    },
    'Alkene_Cyclopropanation': {
        'required': [(['alkene'], []), ([], ['alkene'])],
        'score': 0.7
    },

    # Additional Elimination Reactions (10)
    'E1_Alcohol': {
        'required': [(['alcohol'], []), ([], ['alcohol'])],
        'score': 0.6
    },
    'E2_Alcohol': {
        'required': [(['alcohol'], []), ([], ['alcohol'])],
        'score': 0.6
    },
    'Dehydrochlorination': {
        'required': [(['alkyl_halide'], []), ([], ['alkyl_halide'])],
        'score': 0.6
    },
    'Dehydrosulfurization': {
        'required': [(['thiol'], []), ([], ['thiol'])],
        'score': 0.6
    },
    'E2_Amine': {
        'required': [(['tertiary_amine'], []), ([], ['tertiary_amine'])],
        'score': 0.6
    },
    'Dehydrobromination_Alkyne': {
        'required': [(['alkyl_halide'], []), ([], ['alkyl_halide'])],
        'score': 0.6
    },
    'Selenoxide_Elimination': {
        'required': [(['sulfoxide'], []), ([], ['sulfoxide'])],
        'score': 0.6
    },
    'Dehydroiodination': {
        'required': [(['alkyl_halide'], []), ([], ['alkyl_halide'])],
        'score': 0.6
    },
    'Chugaev_Elimination': {
        'required': [(['alcohol'], []), ([], ['alcohol'])],
        'score': 0.6
    },
    'Corey_Winter': {
        'required': [(['alcohol'], []), ([], ['alcohol'])],
        'score': 0.6
    },

    # Additional Oxidation Reactions (10)
    'Allylic_Oxidation': {
        'required': [(['alkene'], []), ([], ['alkene'])],
        'score': 0.7
    },
    'Oxidative_Cleavage_Alkyne': {
        'required': [(['alkyne'], []), ([], ['alkyne'])],
        'score': 0.6
    },
    'Tosylate_Oxidation': {
        'required': [(['tosylate'], []), ([], ['tosylate'])],
        'score': 0.6
    },
    'SulfideToSulfone': {
        'required': [(['thioether'], []), ([], ['thioether'])],
        'score': 0.6
    },
    'AmineToImine': {
        'required': [(['primary_amine'], []), ([], ['primary_amine'])],
        'score': 0.7
    },
    'AlcoholToAcid_Direct': {
        'required': [(['alcohol'], []), ([], ['alcohol'])],
        'score': 0.7
    },
    'Dess_Martin_Oxidation': {
        'required': [(['alcohol'], []), ([], ['alcohol'])],
        'score': 0.7
    },
    'Oppenauer_Oxidation': {
        'required': [(['alcohol'], ['ketone']), (['ketone'], ['alcohol'])],
        'score': 0.7
    },
    'Oxidation_Thioether': {
        'required': [(['thioether'], []), ([], ['thioether'])],
        'score': 0.6
    },
    'AmineToN_Oxide': {
        'required': [(['tertiary_amine'], []), ([], ['tertiary_amine'])],
        'score': 0.6
    },

    # Additional Reduction Reactions (10)
    'AlkyneToCisAlkene': {
        'required': [(['alkyne'], []), ([], ['alkyne'])],
        'score': 0.7
    },
    'NitrileToAmine_Primary': {
        'required': [(['nitrile'], []), ([], ['nitrile'])],
        'score': 0.7
    },
    'ImineToSecondaryAmine': {
        'required': [(['imine'], []), ([], ['imine'])],
        'score': 0.7
    },
    'KetoneToMethylene': {
        'required': [(['ketone'], []), ([], ['ketone'])],
        'score': 0.6
    },
    'AldehydeToPrimaryAlcohol': {
        'required': [(['aldehyde'], []), ([], ['aldehyde'])],
        'score': 0.7
    },
    'NitrosoToHydroxylamine': {
        'required': [(['nitro'], []), ([], ['nitro'])],
        'score': 0.6
    },
    'AzideToPrimaryAmine': {
        'required': [(['azide'], []), ([], ['azide'])],
        'score': 0.7
    },
    'EsterToAldehyde': {
        'required': [(['ester'], []), ([], ['ester'])],
        'score': 0.7
    },
    'AmideToAldehyde': {
        'required': [(['amide'], []), ([], ['amide'])],
        'score': 0.7
    },
    'DisulfideToThiol_Reduction': {
        'required': [(['disulfide'], []), ([], ['disulfide'])],
        'score': 0.6
    },

    # Additional Cross-Coupling Reactions (10)
    'Sonogashira_Alkyl': {
        'required': [(['alkyl_halide'], ['alkyne']), (['alkyne'], ['alkyl_halide'])],
        'score': 0.8
    },
    'Heck_Alkyl': {
        'required': [(['alkyl_halide'], ['alkene']), (['alkene'], ['alkyl_halide'])],
        'score': 0.8
    },
    'Suzuki_Alkyl_Boronate': {
        'required': [(['alkyl_halide'], ['boronate_ester']), (['boronate_ester'], ['alkyl_halide'])],
        'score': 0.7
    },
    'Stille_Alkyl': {
        'required': [(['alkyl_halide'], ['organotin']), (['organotin'], ['alkyl_halide'])],
        'score': 0.8
    },
    'Negishi_Alkyl': {
        'required': [(['alkyl_halide'], ['organozinc']), (['organozinc'], ['alkyl_halide'])],
        'score': 0.8
    },
    'Kumada_Alkyl': {
        'required': [(['alkyl_halide'], ['organomagnesium']), (['organomagnesium'], ['alkyl_halide'])],
        'score': 0.8
    },
    'Hiyama_Alkyl': {
        'required': [(['alkyl_halide'], ['silane']), (['silane'], ['alkyl_halide'])],
        'score': 0.8
    },
    'Buchwald_Hartwig_Secondary': {
        'required': [(['aryl_halide'], ['secondary_amine']), (['secondary_amine'], ['aryl_halide'])],
        'score': 0.8
    },
    'Cross_Coupling_Vinyl': {
        'required': [(['aryl_halide'], ['alkene']), (['alkene'], ['aryl_halide'])],
        'score': 0.8
    },
    'Miyaura_Borylation_Alkyl': {
        'required': [(['alkyl_halide'], []), ([], ['alkyl_halide'])],
        'score': 0.8
    },

    # Additional Cycloaddition Reactions (10)
    '1,3-Dipolar_Nitrone': {
        'required': [(['alkene'], ['nitrile_oxide']), (['nitrile_oxide'], ['alkene'])],
        'score': 0.7
    },
    'Diels_Alder_Alkyne': {
        'required': [(['alkene'], ['alkyne']), (['alkyne'], ['alkene'])],
        'score': 0.85
    },
    '[2+2]_Photocycloaddition': {
        'required': [(['alkene'], ['alkene'])],
        'score': 0.7
    },
    '1,3-Dipolar_Azomethine': {
        'required': [(['alkene'], []), ([], ['alkene'])],
        'score': 0.7
    },
    'Paterno_Buchi_Alkyne': {
        'required': [(['alkyne'], ['carbonyl']), (['carbonyl'], ['alkyne'])],
        'score': 0.7
    },
    'Cycloaddition_Nitrile': {
        'required': [(['alkene'], ['nitrile']), (['nitrile'], ['alkene'])],
        'score': 0.7
    },
    '[4+2]_Hetero_Diels_Alder_Nitroso': {
        'required': [(['alkene'], ['nitro']), (['nitro'], ['alkene'])],
        'score': 0.7
    },
    '1,5-Dipolar_Cycloaddition': {
        'required': [(['alkene'], []), ([], ['alkene'])],
        'score': 0.7
    },
    'Cycloaddition_Carbonyl_Ylide': {
        'required': [(['alkene'], []), ([], ['alkene'])],
        'score': 0.7
    },
    'Alkyne_Trimerization': {
        'required': [(['alkyne'], ['alkyne'])],
        'score': 0.8
    },

    # Additional Rearrangement Reactions (10)
    'Sigmatropic_1,3': {
        'required': [(['alkene'], []), ([], ['alkene'])],
        'score': 0.6
    },
    'Wolff_Rearrangement': {
        'required': [(['azide'], []), ([], ['azide'])],
        'score': 0.6
    },
    'Schmidt_Rearrangement': {
        'required': [(['azide'], ['carboxylic_acid']), (['carboxylic_acid'], ['azide'])],
        'score': 0.6
    },
    'Baeyer_Villiger_Aldehyde': {
        'required': [(['aldehyde'], []), ([], ['aldehyde'])],
        'score': 0.7
    },
    'Claisen_Allyl_Ether': {
        'required': [(['ether'], []), ([], ['ether'])],
        'score': 0.6
    },
    'Cope_Elimination_Oxide': {
        'required': [(['tertiary_amine'], []), ([], ['tertiary_amine'])],
        'score': 0.6
    },
    'Fischer_Indole': {
        'required': [(['aromatic'], ['hydrazine', 'ketone']), (['hydrazine', 'ketone'], ['aromatic'])],
        'score': 0.7
    },
    'Ireland_Claisen': {
        'required': [(['ester'], []), ([], ['ester'])],
        'score': 0.6
    },
    'Benzilic_Acid': {
        'required': [(['ketone'], []), ([], ['ketone'])],
        'score': 0.6
    },
    'Favorskii_Rearrangement': {
        'required': [(['ketone'], ['alkyl_halide']), (['alkyl_halide'], ['ketone'])],
        'score': 0.6
    },

    # Additional Heterocycle Formation Reactions (10)
    'Pyridine_Formation': {
        'required': [(['carbonyl'], ['amine']), (['amine'], ['carbonyl'])],
        'score': 0.7
    },
    'Thiophene_Formation': {
        'required': [(['ketone'], []), ([], ['ketone'])],
        'score': 0.7
    },
    'Furan_Formation': {
        'required': [(['ketone'], []), ([], ['ketone'])],
        'score': 0.7
    },
    'Pyrazoline_Formation': {
        'required': [(['alkene'], ['hydrazine']), (['hydrazine'], ['alkene'])],
        'score': 0.7
    },
    'Oxazoline_Formation': {
        'required': [(['carboxylic_acid'], ['alcohol', 'amine']), (['alcohol', 'amine'], ['carboxylic_acid'])],
        'score': 0.7
    },
    'Imidazoline_Formation': {
        'required': [(['carboxylic_acid'], ['amine']), (['amine'], ['carboxylic_acid'])],
        'score': 0.7
    },
    'Thiazoline_Formation': {
        'required': [(['carboxylic_acid'], ['thiol', 'amine']), (['thiol', 'amine'], ['carboxylic_acid'])],
        'score': 0.7
    },
    'Indole_Formation_Skatole': {
        'required': [(['aromatic'], ['ketone']), (['ketone'], ['aromatic'])],
        'score': 0.7
    },
    'Quinoline_Formation': {
        'required': [(['aromatic'], ['ketone']), (['ketone'], ['aromatic'])],
        'score': 0.7
    },
    'Isoquinoline_Formation': {
        'required': [(['aromatic'], ['carbonyl']), (['carbonyl'], ['aromatic'])],
        'score': 0.7
    },

    # Additional Protecting Group Chemistry (10)
    'Benzyl_Protection': {
        'required': [(['alcohol'], []), ([], ['alcohol'])],
        'score': 0.6
    },
    'TBS_Protection': {
        'required': [(['alcohol'], []), ([], ['alcohol'])],
        'score': 0.6
    },
    'Trityl_Protection': {
        'required': [(['alcohol'], []), ([], ['alcohol'])],
        'score': 0.6
    },
    'Boc_Protection_Secondary': {
        'required': [(['secondary_amine'], []), ([], ['secondary_amine'])],
        'score': 0.6
    },
    'Fmoc_Protection_Secondary': {
        'required': [(['secondary_amine'], []), ([], ['secondary_amine'])],
        'score': 0.6
    },
    'Methyl_Ether_Protection': {
        'required': [(['alcohol'], []), ([], ['alcohol'])],
        'score': 0.6
    },
    'SEM_Protection': {
        'required': [(['alcohol'], []), ([], ['alcohol'])],
        'score': 0.6
    },
    'Acetal_Protection_Ketone': {
        'required': [(['ketone'], ['alcohol']), (['alcohol'], ['ketone'])],
        'score': 0.6
    },
    'Cbz_Protection_Secondary': {
        'required': [(['secondary_amine'], []), ([], ['secondary_amine'])],
        'score': 0.6
    },
    'PMB_Protection': {
        'required': [(['alcohol'], []), ([], ['alcohol'])],
        'score': 0.6
    },

    # Additional Condensation Reactions (10)
    'Aldol_Addition': {
        'required': [(['aldehyde'], ['aldehyde']), (['aldehyde'], ['aldehyde'])],
        'score': 0.7
    },
    'Cannizzaro_Reaction': {
        'required': [(['aldehyde'], []), ([], ['aldehyde'])],
        'score': 0.7
    },
    'Perkin_Condensation': {
        'required': [(['aromatic'], ['carboxylic_acid']), (['carboxylic_acid'], ['aromatic'])],
        'score': 0.7
    },
    'Knoevenagel_Nitrile': {
        'required': [(['carbonyl'], ['nitrile']), (['nitrile'], ['carbonyl'])],
        'score': 0.7
    },
    'Claisen_Schmidt': {
        'required': [(['aldehyde'], ['ketone']), (['ketone'], ['aldehyde'])],
        'score': 0.7
    },
    'Stobbe_Condensation': {
        'required': [(['carbonyl'], ['ester']), (['ester'], ['carbonyl'])],
        'score': 0.7
    },
    'Tishchenko_Reaction': {
        'required': [(['aldehyde'], ['aldehyde'])],
        'score': 0.7
    },
    'Enamine_Formation': {
        'required': [(['ketone'], ['amine']), (['amine'], ['ketone'])],
        'score': 0.7
    },
    'Acetoacetic_Ester_Condensation': {
        'required': [(['ester'], ['ester'])],
        'score': 0.7
    },
    'Knoevenagel_Ester': {
        'required': [(['carbonyl'], ['ester']), (['ester'], ['carbonyl'])],
        'score': 0.7
    },

    # Additional Multi-Component Reactions (10)
    'Ugi_Three_Component': {
        'required': [(['carbonyl'], ['amine', 'isocyanate']), (['amine', 'isocyanate'], ['carbonyl'])],
        'score': 0.8
    },
    'Passerini_Aldehyde': {
        'required': [(['aldehyde'], ['carboxylic_acid', 'isocyanate']), 
                     (['carboxylic_acid', 'isocyanate'], ['aldehyde'])],
        'score': 0.7
    },
    'Biginelli_Urea': {
        'required': [(['aldehyde'], ['urea', 'ketone']), (['urea', 'ketone'], ['aldehyde'])],
        'score': 0.7
    },
    'Hantzsch_Pyrrole': {
        'required': [(['ketone'], ['amine']), (['amine'], ['ketone'])],
        'score': 0.7
    },
    'Gewald_Thiazole': {
        'required': [(['ketone'], ['nitrile']), (['nitrile'], ['ketone'])],
        'score': 0.7
    },
    'Povarov_Reaction': {
        'required': [(['aromatic'], ['alkene', 'aldehyde']), (['alkene', 'aldehyde'], ['aromatic'])],
        'score': 0.7
    },
    'Mannich_Imine': {
        'required': [(['aldehyde'], ['imine', 'carbonyl']), (['imine', 'carbonyl'], ['aldehyde'])],
        'score': 0.7
    },
    'Strecker_Nitrile': {
        'required': [(['aldehyde'], ['amine', 'nitrile']), (['amine', 'nitrile'], ['aldehyde'])],
        'score': 0.7
    },
    'Paal_Knorr_Furan': {
        'required': [(['ketone'], ['alcohol']), (['alcohol'], ['ketone'])],
        'score': 0.7
    },
    'Hantzsch_Thiazole': {
        'required': [(['ketone'], ['thiourea']), (['thiourea'], ['ketone'])],
        'score': 0.7
    },

    # Additional Specialized Reactions (10)
    'Wacker_Oxidation_Alkyne': {
        'required': [(['alkyne'], []), ([], ['alkyne'])],
        'score': 0.7
    },
    'Ring_Closing_Metathesis': {
        'required': [(['alkene'], ['alkene'])],
        'score': 0.8
    },
    'Pauson_Khand_Alkyne': {
        'required': [(['alkyne'], ['alkene', 'carbonyl']), (['alkene', 'carbonyl'], ['alkyne'])],
        'score': 0.7
    },
    'Hydroacylation': {
        'required': [(['aldehyde'], ['alkene']), (['alkene'], ['aldehyde'])],
        'score': 0.7
    },
    'C_H_Insertion': {
        'required': [(['aromatic'], []), ([], ['aromatic'])],
        'score': 0.7
    },
    'Metathesis_Cross': {
        'required': [(['alkene'], ['alkene'])],
        'score': 0.8
    },
    'Hydroarylation': {
        'required': [(['aromatic'], ['alkene']), (['alkene'], ['aromatic'])],
        'score': 0.7
    },
    'Click_Reaction_Thiol': {
        'required': [(['alkene'], ['thiol']), (['thiol'], ['alkene'])],
        'score': 0.7
    },
    'Oxidative_Coupling': {
        'required': [(['aromatic'], ['aromatic'])],
        'score': 0.7
    },
    'Hydroboration_Oxidation': {
        'required': [(['alkene'], []), ([], ['alkene'])],
        'score': 0.7
    }
    }

    if reaction_type not in compatibility_rules:
        return 0.5
    rule = compatibility_rules[reaction_type]
    for req_pair in rule['required']:
        if (any(g in groups1 for g in req_pair[0]) and any(g in groups2 for g in req_pair[1])) or \
           (any(g in groups2 for g in req_pair[0]) and any(g in groups1 for g in req_pair[1])):
            return rule['score']
    return 0.1
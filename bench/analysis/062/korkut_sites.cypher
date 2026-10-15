// specs/062: phosphorylated residues Reactome curates for the korkut2015
// phospho-readout genes (human EWAS), to fix each antibody's coordinate.
MATCH (e:EntityWithAccessionedSequence)-[:referenceEntity]->(r:ReferenceGeneProduct),
      (e)-[:hasModifiedResidue]->(m)-[:psiMod]->(p)
WHERE e.speciesName = 'Homo sapiens' AND p.name[0] IN ['O-phospho-L-serine','O-phospho-L-threonine',"O4'-phospho-L-tyrosine"]
  AND r.geneName[0] IN ['EIF4EBP1','RB1','MAPK1','MAPK3','MAP2K1','MAP2K2','RPS6','AKT1','AKT2','AKT3','PRKAA1','PRKAA2','CTNNB1','GSK3A','GSK3B','MAPK14','RPS6KB1','SRC','STAT3','STAT5A','STAT5B','STAT6','TSC2','YAP1','YBX1','CHEK2','ATM','JUN','CHEK1','ELK1','MTOR','WWTR1','PDPK1','ACACA','IRS1','SMAD3']
RETURN DISTINCT r.geneName[0] AS gene, r.identifier AS uniprot, substring(p.name[0], 8, 1) AS res, m.coordinate AS coord
ORDER BY gene, coord;

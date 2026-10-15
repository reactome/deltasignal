:param ids => ["R-HSA-8938854","R-HSA-912690","R-HSA-449798","R-HSA-2473196","R-HSA-2473219","R-HSA-1181341","R-HSA-1181340","R-HSA-1181346","R-HSA-1911561","R-HSA-1470001","R-HSA-1562562","R-HSA-2187341","R-HSA-2187343","R-HSA-201806","R-HSA-5656063","R-HSA-508158","R-HSA-508159","R-HSA-937333","R-HSA-182945","R-HSA-201841","R-HSA-202644","R-HSA-201446","R-HSA-1433442","R-HSA-5635085"];
WITH $ids AS ids UNWIND ids AS id
MATCH (c:PhysicalEntity {stId:id})
OPTIONAL MATCH (rx:ReactionLikeEvent)-[:input]->(c)
OPTIONAL MATCH (cat:CatalystActivity)-[:physicalEntity]->(c)
OPTIONAL MATCH (reg:Regulation)-[:regulator]->(c)
OPTIONAL MATCH (prod:ReactionLikeEvent)-[:output]->(c)
OPTIONAL MATCH (prod)-[:goBiologicalProcess]->(go)
RETURN id, c.displayName AS name, count(DISTINCT rx) AS used_as_input, count(DISTINCT cat) AS catalyst,
 count(DISTINCT reg) AS regulator, [x IN collect(DISTINCT reg.schemaClass) WHERE x IS NOT NULL] AS regclass,
 [x IN collect(DISTINCT go.name) WHERE x IS NOT NULL][..3] AS go;
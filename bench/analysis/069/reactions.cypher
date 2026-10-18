UNWIND $pids AS pid
MATCH (p:Pathway {stId:pid})-[:hasEvent*]->(r:ReactionLikeEvent)
WITH DISTINCT pid, r
OPTIONAL MATCH (r)-[:input]->(i)
OPTIONAL MATCH (r)-[:output]->(o)
OPTIONAL MATCH (r)-[:precedingEvent]->(pe:ReactionLikeEvent)
RETURN pid, r.stId, r.schemaClass AS cls, apoc.text.join(collect(DISTINCT i.stId),";") AS ins, apoc.text.join(collect(DISTINCT o.stId),";") AS outs, apoc.text.join(collect(DISTINCT pe.stId),";") AS prev;

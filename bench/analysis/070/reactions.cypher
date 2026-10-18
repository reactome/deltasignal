// specs/070: every reaction under a catalog pathway, with all curated participants.
UNWIND ['R-HSA-5673001','R-HSA-1227986','R-HSA-5693606','R-HSA-5693567','R-HSA-5685939','R-HSA-5693571','R-HSA-5358508','R-HSA-6783310','R-HSA-5696398','R-HSA-168928','R-HSA-73942','R-HSA-73893','R-HSA-73884','R-HSA-194315','R-HSA-2559583','R-HSA-2871796','R-HSA-1236394','R-HSA-450531','R-HSA-195721','R-HSA-8848021','R-HSA-69620','R-HSA-453279','R-HSA-1538133','R-HSA-453274','R-HSA-3700989','R-HSA-109606','R-HSA-170834','R-HSA-5358351','R-HSA-373755','R-HSA-375165','R-HSA-373752','R-HSA-446728','R-HSA-376176','R-HSA-373760','R-HSA-2682334','R-HSA-1257604','R-HSA-9012999','R-HSA-5625740','R-HSA-5627117','R-HSA-5654736','R-HSA-5654738','R-HSA-5654741','R-HSA-5654743','R-HSA-5625900','R-HSA-5693532','R-HSA-1912422','R-HSA-1980143','R-HSA-1980145','R-HSA-1502540','R-HSA-1181150','R-HSA-201451','R-HSA-194138','R-HSA-140534','R-HSA-418889','R-HSA-6806834','R-HSA-75153','R-HSA-74752','R-HSA-3247509','R-HSA-8878159','R-HSA-8878171','R-HSA-8864260','R-HSA-2172127','R-HSA-6785807','R-HSA-428359','R-HSA-5221030','R-HSA-388841','R-HSA-177929','R-HSA-2404192','R-HSA-1433557','R-HSA-186797','R-HSA-1266695','R-HSA-8878166','R-HSA-452723','R-HSA-114604','R-HSA-5626467','R-HSA-68875','R-HSA-373753','R-HSA-69242','R-HSA-451927','R-HSA-8854691','R-HSA-8853659','R-HSA-450294','R-HSA-512988','R-HSA-5617472','R-HSA-983169','R-HSA-198933','R-HSA-877300','R-HSA-909733','R-HSA-427413','R-HSA-73894','R-HSA-6794361','R-HSA-186712'] AS pid
MATCH (p:Pathway {stId:pid})-[:hasEvent*]->(r:ReactionLikeEvent)
WITH DISTINCT pid, r
OPTIONAL MATCH (r)-[:input]->(i)
WITH pid, r, collect(DISTINCT i.stId) AS ins
OPTIONAL MATCH (r)-[:output]->(o)
WITH pid, r, ins, collect(DISTINCT o.stId) AS outs
OPTIONAL MATCH (r)-[:catalystActivity]->(:CatalystActivity)-[:physicalEntity]->(c)
WITH pid, r, ins, outs, collect(DISTINCT c.stId) AS cats
OPTIONAL MATCH (r)-[:regulatedBy]->(g:PositiveRegulation)-[:regulator]->(x)
WITH pid, r, ins, outs, cats, collect(DISTINCT x.stId) AS pos
OPTIONAL MATCH (r)-[:regulatedBy]->(g:NegativeRegulation)-[:regulator]->(x)
WITH pid, r, ins, outs, cats, pos, collect(DISTINCT x.stId) AS neg
OPTIONAL MATCH (r)-[:precedingEvent]->(pe:ReactionLikeEvent)
RETURN pid, r.stId AS rid, r.schemaClass AS cls,
  apoc.text.join(ins,";") AS ins, apoc.text.join(outs,";") AS outs, apoc.text.join(cats,";") AS cats,
  apoc.text.join(pos,";") AS pos, apoc.text.join(neg,";") AS neg,
  apoc.text.join(collect(DISTINCT pe.stId),";") AS prev;

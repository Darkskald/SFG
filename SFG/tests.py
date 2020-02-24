from SFG.orm.orm import PostProcessor, WorkDatabaseWizard, Substances
import json

w = WorkDatabaseWizard()
p = PostProcessor(w, new=False)

#name = "20161111_NA_2ul"
#print(p.refine_regular(name))

with open("newport/substances.json") as infile:
    substances = json.load(infile)
    substances_orm = []
    for substance in substances:
        s = Substances()
        for key in substance:
            setattr(s, key, substance[key])
        substances_orm.append(s)
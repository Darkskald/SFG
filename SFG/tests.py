from SFG.orm import PostProcessor, WorkDatabaseWizard

w = WorkDatabaseWizard()
p = PostProcessor(w, new=False)

name = "20161111_NA_2ul"
print(p.refine_regular(name))



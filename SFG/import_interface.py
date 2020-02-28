from SFG.orm.orm import ImportDatabaseWizard, PostProcessor
from SFG.natural_samples.gasex import GasExManager
from SFG.natural_samples.boknis_eck import BoknisEckExtension

import os

np = "C:/Users/lange/Desktop/CharmingSFG/SFG/newport"
db = "C:/Users/lange/Desktop/CharmingSFG/SFG/orm.db"

if __name__ == "__main__":
    os.remove("C:/Users/lange/Desktop/CharmingSFG/SFG/orm.db")
    D = ImportDatabaseWizard()
    P = PostProcessor(D)
    P.disconnect()
    GasExManager()
    BoknisEckExtension(new=True)
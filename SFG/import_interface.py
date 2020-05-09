import os
import timeit

# from SFG.natural_samples.boknis_eck import BoknisEckExtension
# from SFG.natural_samples.gasex import GasExManager
from SFG.orm.orm import ImportDatabaseWizard, PostProcessor

np = "C:/Users/lange/Desktop/CharmingSFG/SFG/newport"
db = "C:/Users/lange/Desktop/CharmingSFG/SFG/orm.db"

if __name__ == "__main__":
    os.remove(os.getcwd() + "/orm.db")
    # D = ImportDatabaseWizard()
    # P = PostProcessor(D)
    # P.disconnect()
    # GasExManager()
    # BoknisEckExtension(new=True)
    start = timeit.default_timer()
    D = ImportDatabaseWizard()
    end = timeit.default_timer()
    print(end - start)

import os
import timeit

# from SFG.natural_samples.boknis_eck import BoknisEckExtension
# from SFG.natural_samples.gasex import GasExManager
from SFG.orm.import_db_controller import ImportDatabaseWizard
import logging
np = "C:/Users/lange/Desktop/CharmingSFG/SFG/newport"
db = "C:/Users/lange/Desktop/CharmingSFG/SFG/orm.db"

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    os.remove(os.getcwd() + "/orm.db")
    # GasExManager()
    # BoknisEckExtension(new=True)
    start = timeit.default_timer()
    D = ImportDatabaseWizard()
    end = timeit.default_timer()
    print(end - start)


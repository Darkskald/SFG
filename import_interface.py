from orm import ImportDatabaseWizard, PostProcessor
from gasex import GasExManager
from boknis_eck import BoknisEckExtension

import os

if __name__ == "__main__":
    D = ImportDatabaseWizard()
    P = PostProcessor(D)
    P.disconnect()
    os.chdir("..")
    GasExManager()
    BoknisEckExtension(new=True)
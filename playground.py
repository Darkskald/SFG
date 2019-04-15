import datetime
import re
import pandas as pd

def tprint(*args):
    for i in args:
        print(i)
        print(type(i))


def to_date(s):
    year = int(s[5:9])
    month = int(s[10:12])
    day = int(s[13:15])

    return datetime.date(year, month, day)


def is_date(s):
    dre = re.compile('SFG \(\d{4}-\d{2}-\d{2}\)')

    try:
        out = re.match(dre, s)

        if out is not None:
            return True
        else:
            return False

    except TypeError:

        return False


df = pd.read_excel("Wasserproben_uc.xlsx", sheet_name="Samples", header=2)

sel = "Experiment dates"
sel2 = "Date"
sel3 = "Sample"

test = df[sel]

ts = "SFG (2017-03-14)"

# use boolean indexing by a column preprocessed by applying a function to a row
abc = df[sel].apply(is_date)
a = df[df[sel].apply(is_date)]
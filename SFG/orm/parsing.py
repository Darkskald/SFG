import datetime
import re
from typing import Dict


# todo: this module requires documentation

class InvalidNameError(Exception):
    pass


lt_regex = {
    "sample_no": re.compile(r'x\d'),
    "measurement_no": re.compile(r'#\d'),
    "ratio": re.compile(r'\dto\d'),
    "surfactant": re.compile(r'\wA'),
    "sensitizer": re.compile(r'BX\d'),
    "speed": re.compile(r'^\d.\d*$'),
    "spreading_volume": re.compile(r'\d\d(ul)?$'),
    "conc": re.compile(r'\d+(.\d)?mM$')
}

# for parsing of Boknis Eck systematic names

n1 = re.compile("\D\d{1,2}\D")
n2 = re.compile("[ a-zA-Z_-]\d{1,2}$")
n3 = re.compile("\d{1,2}-?#$")
n4 = re.compile("-#\d{1,2}$")

numreg = [n1, n2, n3, n4]
# sampling dates to extract the sampling date of a sample from the filename:
d1 = re.compile("^\d{8}_\d{1,2}\D")
d2 = re.compile("_[a-zA-z -]*\d{8}")
d3 = re.compile("^\d{8}_[a-zA-Z]*[ -]")
d4 = re.compile("^\d{8}_\d{1,2}$")
d5 = re.compile("^\d{8}_[a-zA-Z]{2}\d{1,2}-")

datereg = [d1, d2, d3, d4, d5]


def refine_regular_lt(string):
    out = {}
    for token in string.split("_"):
        for item in lt_regex:
            temp = re.match(lt_regex[item], token)
            if temp is not None:
                out[item] = token
    return out


def refine_regular(namestring, substances, yield_full=False):
    """This function parses the names of the spectra in the sfg table with the type "regular" in
    order to extract additional metainformation. It makes use of regular expressions."""

    process_list = namestring.split("_")
    sample = re.compile('^x\d{1,2}$')
    measurement = re.compile('^#\d{1,2}$')
    photolysis = re.compile('^\d{1,3}p$')
    spread_vol = re.compile('^\d{1,2}(.\d{1,2})?$')
    conc = re.compile('^\d{1,2}mM$')
    ratio_reg = re.compile('^\dto\d$')

    date = datetime.date(int(process_list[0][0:4]), int(process_list[0][4:6]),
                         int(process_list[0][6:]))
    surf = None
    sens = None
    surf_v = None
    sens_v = None
    surf_c = None
    sens_c = None
    sample_nr = None
    measurement_nr = None
    photo = None
    ratio = None
    comment = None

    for item in process_list[1:]:

        if item in substances:
            if surf is None:
                surf = item
            else:
                if substances[item] == "y":
                    sens = item

        elif re.match(sample, item):
            sample_nr = item

        elif re.match(measurement, item):
            measurement_nr = item

        elif re.match(photolysis, item):
            photo = item

        elif re.match(conc, item):

            if surf_c is None:
                surf_c = item
            else:
                sens_c = item

        elif re.match(spread_vol, item):

            if surf_v is None:
                surf_v = item
            else:
                sens_v = item

        elif re.match(ratio_reg, item):
            ratio = item

        else:
            comment = item

    dic = {"sensitizer": sens, "sample_no": sample_nr, "measurement_no": measurement_nr,
           "surfactant_conc": surf_c, "sensitizer_conc": sens_c, "surfactant": surf, "surfactant_vol": surf_v,
           "sensitizer_vol": sens_v, "comment": comment, "photolysis": photo, "ratio": ratio}

    if yield_full:
        dic["date"] = date
        dic["full_name"]: namestring

    return dic


def get_hashes(name: str) -> Dict[str, str]:
    """ Extract the hashes from file names. Remember to strip the file ending
    before passing namestring to this function."""

    try:
        temp = name.split("_")[1:]
        station_hash = temp[0] + temp[1][1]
        if temp[1][0] != "c":
            sample_hash = temp[0] + temp[1] + temp[2] + temp[3]

        else:
            sample_hash = temp[0] + temp[1] + "deep"
        dic = {}
        dic["name"] = name
        dic["station_hash"] = station_hash
        dic["sample_hash"] = sample_hash
        return dic
    except IndexError:
        raise InvalidNameError(f'{name} is not well-formed!')


def process_station_hash(station_hash):
    """Extracts metadate from the station hash"""
    date = station_hash[:4]
    date = datetime.date(2018, int(date[0:2]), int(date[2:]))
    number = station_hash[4]
    return {"date": date, "number": number}


# todo: this hack needs explanation
def get_station_type(sample_hash):
    temp = sample_hash[4]

    if temp in ("c", "r"):
        stype = "big"
    else:
        stype = "small"

    return stype


def process_sample_hash(sample_hash) -> Dict:
    """Extracts metadata from the sample hash"""
    dic = {}
    location = sample_hash[4]
    dic["location"] = location

    if location == "c":
        if "deep" in sample_hash or "low" in sample_hash:
            dic["type"] = "deep"
            dic["number"] = "-"

    elif location == "r":
        dic["type"] = sample_hash[6]
        if dic["type"] in ("p", "P"):
            dic["number"] = sample_hash[7]
            dic["type"] = "p"

        elif dic["type"] == "s":
            dic["number"] = sample_hash[8]

    elif location == "a":
        dic["type"] = sample_hash[6]
        dic["number"] = sample_hash[8]

    else:
        raise ValueError(f"Invalid sample hash {sample_hash}")

    return dic


def sample_hash_to_station_hash(sample_hash: str) -> str:
    """This function calculates the station hash from a given sample hash"""
    temp = sample_hash[:4] + sample_hash[5]
    return temp

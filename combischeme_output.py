from icecream import ic
import json


# cf. https://stackoverflow.com/a/31631711
def readable_bytes(B):
    """Return the given bytes as a human friendly KiB, MiB, GB, or TiB string."""
    B = float(B)
    KiB = float(1024)
    MiB = float(KiB ** 2)  # 1,048,576
    GiB = float(KiB ** 3)  # 1,073,741,824
    TiB = float(KiB ** 4)  # 1,099,511,627,776

    if B < KiB:
        return '{0}_{1}'.format(B, 'Bytes' if 0 == B > 1 else 'Byte')
    elif KiB <= B < MiB:
        return '{0:.2f}_KiB'.format(B / KiB)
    elif MiB <= B < GiB:
        return '{0:.2f}_MiB'.format(B / MiB)
    elif GiB <= B < TiB:
        return '{0:.2f}_GiB'.format(B / GiB)
    elif TiB <= B:
        return '{0:.2f}_TiB'.format(B / TiB)


def read_data_from_json(filename):
    with open(filename, 'r') as f:
        data = json.load(f)
    return data


def write_scheme_dictionary_to_json(dictionary: dict, filename):
    schemeList = []
    for key, value in dictionary.items():
        schemeList += [{"coeff": value, "level": list(key)}]

    # ic(schemeList)
    jsonString = json.dumps(schemeList)  # , indent=0)

    with open(filename, 'w') as f:
        f.write(jsonString)

from icecream import ic
import combischeme_utils
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


def write_scheme_to_json_without_process_group_number(scheme: combischeme_utils.CombinationScheme):
    schemeList = []
    for key, value in scheme.get_combination_dictionary().items():
        schemeList += [{"coeff": value, "level": list(key)}]

    # ic(schemeList)
    jsonString = json.dumps(schemeList)  # , indent=0)

    mem = (scheme.get_total_num_points_combi()*8)
    # ic(readable_bytes(mem))
    with open('scheme_' + readable_bytes(mem) + '.json', 'w') as f:
        f.write(jsonString)

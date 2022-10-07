import combischeme_utils


#cf. https://stackoverflow.com/a/31631711
def readable_bytes(B):
    """Return the given bytes as a human friendly KiB, MiB, GB, or TiB string."""
    B = float(B)
    KiB = float(1024)
    MiB = float(KiB ** 2) # 1,048,576
    GiB = float(KiB ** 3) # 1,073,741,824
    TiB = float(KiB ** 4) # 1,099,511,627,776

    if B < KiB:
        return '{0} {1}'.format(B,'Bytes' if 0 == B > 1 else 'Byte')
    elif KiB <= B < MiB:
        return '{0:.2f} KiB'.format(B / KiB)
    elif MiB <= B < GiB:
        return '{0:.2f} MiB'.format(B / MiB)
    elif GiB <= B < TiB:
        return '{0:.2f} GB'.format(B / GiB)
    elif TiB <= B:
        return '{0:.2f} TiB'.format(B / TiB)

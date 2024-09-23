from collections import defaultdict

from docutils.nodes import header



def assign_bin(size, bins):
    for i, b in enumerate(bins):
        if size < b:
            return i - 1
    return len(bins) - 1
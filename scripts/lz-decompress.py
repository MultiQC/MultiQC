import sys

from multiqc.utils import lzstring


if __name__ == "__main__":
    lz_string = lzstring.LZString()
    with open(sys.argv[1], "rt") as f:
        data = f.read()
    print(lz_string.decompressFromBase64(data), end="")

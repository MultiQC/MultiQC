#!/usr/bin/python

"""
lz-string for python v1.0.4 (https://github.com/gkovacs/lz-string-python)

(c) Geza Kovacs

+ Python 2 dropped, Python 3.12 added

License: https://github.com/gkovacs/lz-string-python/blob/master/LICENSE.md
"""
import base64
import io
import math
import struct
from dataclasses import dataclass

keyStrBase64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/="
keyStrUriSafe = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-$"
baseReverseDic = {}


@dataclass
class Entry:
    val: int
    position: int
    index: int


def reverse_8bit_int(value: int) -> int:
    """Reverses the bit order of an 8-bit int"""
    value &= 0xFF  # Ensure 8-bit
    value = ((value & 0b11110000) >> 4) | ((value & 0b00001111) << 4)
    value = ((value & 0b11001100) >> 2) | ((value & 0b00110011) << 2)
    value = ((value & 0b10101010) >> 1) | ((value & 0b01010101) << 1)
    return value


# Since bitwise operations are expensive create a lookup table
REVERSE_LOOKUP = bytes(reverse_8bit_int(i) for i in range(256))


assert REVERSE_LOOKUP[0b10000000] == 0b00000001
assert REVERSE_LOOKUP[0b10101010] == 0b01010101


class ContextWriter:
    def __init__(self):
        self.bitbuffer = 0
        self.bits_stored = 0
        self.buffer = io.BytesIO()
        self.context_enlarge_in = 2  # Compensate for the first entry which should not count
        self.context_numBits = 2

    def write(self, value: int, bits: int = 0):
        # Ensure no rogue bits are set after the bits of interest
        bits = bits or self.context_numBits
        bitmask = (1 << bits) - 1
        value &= bitmask

        # Python can store integers bigger than 64 bits. So no checks for overflow.
        self.bitbuffer |= value << self.bits_stored
        self.bits_stored += bits

        self.context_enlarge_in -= 1
        if self.context_enlarge_in == 0:
            self.context_enlarge_in = 1 << self.context_numBits
            self.context_numBits += 1

        if self.bits_stored >= 64:
            integer_to_store = self.bitbuffer & 0xFFFF_FFFF_FFFF_FFFF
            self.bitbuffer >>= 64
            self.bits_stored -= 64
            # The integer needs to be stored in reverse. First use little
            # endian order to reverse the byte order. Then use a lookup table
            # to inverse the individual bytes. This is faster than doing
            # a lot of bitwise operations.
            integer_bytes = struct.pack("<Q", integer_to_store)
            self.buffer.write(integer_bytes.translate(REVERSE_LOOKUP))

    def getvalue(self):
        value = self.buffer.getvalue()
        # Flush the last char
        leftover = struct.pack("<Q", self.bitbuffer)
        leftover = leftover.translate(REVERSE_LOOKUP)
        bytes_to_use = (self.bits_stored + 7) // 8
        return value + leftover[:bytes_to_use]


def getBaseValue(alphabet, character):
    if alphabet not in baseReverseDic:
        baseReverseDic[alphabet] = {}
    for i in range(len(alphabet)):
        baseReverseDic[alphabet][alphabet[i]] = i
    return baseReverseDic[alphabet][character]


def _compress(uncompressed):
    if uncompressed is None:
        return ""

    if isinstance(uncompressed, bytes):
        uncompressed = uncompressed.decode()

    context_dictionary = {}
    context_dictionaryToCreate = {}
    context_c = ""
    context_wc = ""
    context_w = ""
    context_dictSize = 3
    context_data = ContextWriter()

    for context_c in uncompressed:
        if context_c not in context_dictionary:
            context_dictionary[context_c] = context_dictSize
            context_dictSize += 1
            context_dictionaryToCreate[context_c] = True

        context_wc = context_w + context_c
        if context_wc in context_dictionary:
            context_w = context_wc
        else:
            if context_w in context_dictionaryToCreate:
                value = ord(context_w[0])
                if value < 256:
                    context_data.write(0)
                    context_data.write(value, 8)
                else:
                    context_data.write(1)
                    context_data.write(value, 16)
                del context_dictionaryToCreate[context_w]
            else:
                value = context_dictionary[context_w]
                context_data.write(value)

            # Add wc to the dictionary.
            context_dictionary[context_wc] = context_dictSize
            context_dictSize += 1
            context_w = str(context_c)

    # Output the code for w.
    if context_w != "":
        if context_w in context_dictionaryToCreate:
            value = ord(context_w[0])
            if value < 256:
                context_data.write(0)
                context_data.write(value, 8)
            else:
                context_data.write(1)
                context_data.write(value, 16)
            del context_dictionaryToCreate[context_w]
        else:
            value = context_dictionary[context_w]
            context_data.write(value)

    # Mark the end of the stream
    context_data.write(2)

    return context_data.getvalue()


def _decompress(length, resetValue, getNextValue):
    dictionary = {}
    enlargeIn = 4
    dictSize = 4
    numBits = 3
    entry = ""
    result = []

    data = Entry(val=getNextValue(0), position=resetValue, index=1)

    for i in range(3):
        dictionary[i] = i

    bits = 0
    maxpower = math.pow(2, 2)
    power = 1

    while power != maxpower:
        resb = data.val & data.position
        data.position >>= 1
        if data.position == 0:
            data.position = resetValue
            data.val = getNextValue(data.index)
            data.index += 1

        bits |= power if resb > 0 else 0
        power <<= 1

    next = bits
    if next == 0:
        bits = 0
        maxpower = math.pow(2, 8)
        power = 1
        while power != maxpower:
            resb = data.val & data.position
            data.position >>= 1
            if data.position == 0:
                data.position = resetValue
                data.val = getNextValue(data.index)
                data.index += 1
            bits |= power if resb > 0 else 0
            power <<= 1
        c = chr(bits)
    elif next == 1:
        bits = 0
        maxpower = math.pow(2, 16)
        power = 1
        while power != maxpower:
            resb = data.val & data.position
            data.position >>= 1
            if data.position == 0:
                data.position = resetValue
                data.val = getNextValue(data.index)
                data.index += 1
            bits |= power if resb > 0 else 0
            power <<= 1
        c = chr(bits)
    elif next == 2:
        return ""

    dictionary[3] = c
    w = c
    result.append(c)
    counter = 0
    while True:
        counter += 1
        if data.index > length:
            return ""

        bits = 0
        maxpower = math.pow(2, numBits)
        power = 1
        while power != maxpower:
            resb = data.val & data.position
            data.position >>= 1
            if data.position == 0:
                data.position = resetValue
                data.val = getNextValue(data.index)
                data.index += 1
            bits |= power if resb > 0 else 0
            power <<= 1

        c = bits
        if c == 0:
            bits = 0
            maxpower = math.pow(2, 8)
            power = 1
            while power != maxpower:
                resb = data.val & data.position
                data.position >>= 1
                if data.position == 0:
                    data.position = resetValue
                    data.val = getNextValue(data.index)
                    data.index += 1
                bits |= power if resb > 0 else 0
                power <<= 1

            dictionary[dictSize] = chr(bits)
            dictSize += 1
            c = dictSize - 1
            enlargeIn -= 1
        elif c == 1:
            bits = 0
            maxpower = math.pow(2, 16)
            power = 1
            while power != maxpower:
                resb = data.val & data.position
                data.position >>= 1
                if data.position == 0:
                    data.position = resetValue
                    data.val = getNextValue(data.index)
                    data.index += 1
                bits |= power if resb > 0 else 0
                power <<= 1
            dictionary[dictSize] = chr(bits)
            dictSize += 1
            c = dictSize - 1
            enlargeIn -= 1
        elif c == 2:
            return "".join(result)

        if enlargeIn == 0:
            enlargeIn = math.pow(2, numBits)
            numBits += 1

        if c in dictionary:
            entry = dictionary[c]
        else:
            if c == dictSize:
                entry = w + w[0]
            else:
                return None
        result.append(entry)

        # Add w+entry[0] to the dictionary.
        dictionary[dictSize] = w + entry[0]
        dictSize += 1
        enlargeIn -= 1

        w = entry
        if enlargeIn == 0:
            enlargeIn = math.pow(2, numBits)
            numBits += 1


class LZString:
    @staticmethod
    def compressToBase64(uncompressed):
        if uncompressed is None:
            return ""
        res = _compress(uncompressed)
        return base64.b64encode(res).decode("ascii")

    @staticmethod
    def decompress(compressed):
        if compressed is None:
            return ""
        if compressed == "":
            return None
        return _decompress(len(compressed), 32768, lambda index: ord(compressed[index]))

    @staticmethod
    def decompressFromUint8Array(compressed):
        if compressed is None:
            return ""
        if compressed == "":
            return None
        return _decompress(len(compressed), 128, lambda index: compressed[index])

    @staticmethod
    def decompressFromUTF16(compressed):
        if compressed is None:
            return ""
        if compressed == "":
            return None
        return _decompress(len(compressed), 16384, lambda index: compressed[index] - 32)

    @staticmethod
    def decompressFromBase64(compressed):
        if compressed is None:
            return ""
        if compressed == "":
            return None
        return _decompress(len(compressed), 32, lambda index: getBaseValue(keyStrBase64, compressed[index]))

    @staticmethod
    def decompressFromEncodedURIComponent(compressed):
        if compressed is None:
            return ""
        if compressed == "":
            return None
        compressed = compressed.replace(" ", "+")
        return _decompress(len(compressed), 32, lambda index: getBaseValue(keyStrUriSafe, compressed[index]))

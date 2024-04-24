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


class BitWriter:
    def __init__(self):
        self.bit_store = 0
        self.bits_remaining = 8
        self.buffer = io.BytesIO()

    def write(self, value: int, bits: int):
        # Retrieve values from self to save later. This saves a lot of
        # dictionary lookups.
        bit_store = self.bit_store
        bits_remaining = self.bits_remaining
        for i in range(bits):
            bit_store <<= 1
            bit_store |= value & 1
            value >>= 1
            bits_remaining -= 1
            if bits_remaining == 0:
                self.buffer.write(struct.pack("B", self.bit_store))
                bit_store = 0
                bits_remaining = 8
        self.bit_store = bit_store
        self.bits_remaining = bits_remaining

    def getvalue(self):
        value = self.buffer.getvalue()
        # Flush the last char
        leftover = self.bit_store << self.bits_remaining
        return value + struct.pack("B", leftover)


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
    context_enlargeIn = 2  # Compensate for the first entry which should not count
    context_dictSize = 3
    context_numBits = 2
    context_data = BitWriter()

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
                if ord(context_w[0]) < 256:
                    context_data.write(0, context_numBits)
                    context_data.write(value, 8)
                else:
                    context_data.write(1, context_numBits)
                    context_data.write(value, 16)
                context_enlargeIn -= 1
                if context_enlargeIn == 0:
                    context_enlargeIn = 1 << context_numBits
                    context_numBits += 1
                del context_dictionaryToCreate[context_w]
            else:
                value = context_dictionary[context_w]
                context_data.write(value, context_numBits)
            context_enlargeIn -= 1
            if context_enlargeIn == 0:
                context_enlargeIn = 1 << context_numBits
                context_numBits += 1

            # Add wc to the dictionary.
            context_dictionary[context_wc] = context_dictSize
            context_dictSize += 1
            context_w = str(context_c)

    # Output the code for w.
    if context_w != "":
        if context_w in context_dictionaryToCreate:
            value = ord(context_w[0])
            if ord(context_w[0]) < 256:
                context_data.write(0, context_numBits)
                context_data.write(value, 8)
            else:
                context_data.write(1, context_numBits)
                context_data.write(value, 16)
            context_enlargeIn -= 1
            if context_enlargeIn == 0:
                context_enlargeIn = 1 << context_numBits
                context_numBits += 1
            del context_dictionaryToCreate[context_w]
        else:
            value = context_dictionary[context_w]
            context_data.write(value, context_numBits)

    context_enlargeIn -= 1
    if context_enlargeIn == 0:
        context_enlargeIn = 1 << context_numBits
        context_numBits += 1

    # Mark the end of the stream
    context_data.write(2, context_numBits)

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

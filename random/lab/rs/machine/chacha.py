import re
import sys

SIGMA = [0x61707865, 0x3320646E, 0x79622D32, 0x6B206574]
MASK = 0xFFFFFFFF
RFC = "e4e7f110 15593bd1 1fdd0f50 c47120a3 c7f4d1c7 0368c033 9aaa2204 4e6cd4c3 466482d2 09aa9f07 05d7c214 a2028bd9 d19c12b5 b94e16de e883d0cb 4e3c50a2"


def rotl(v, c):
    return ((v << c) | (v >> (32 - c))) & MASK


def quarter(s, a, b, c, d):
    s[a] = (s[a] + s[b]) & MASK
    s[d] = rotl(s[d] ^ s[a], 16)
    s[c] = (s[c] + s[d]) & MASK
    s[b] = rotl(s[b] ^ s[c], 12)
    s[a] = (s[a] + s[b]) & MASK
    s[d] = rotl(s[d] ^ s[a], 8)
    s[c] = (s[c] + s[d]) & MASK
    s[b] = rotl(s[b] ^ s[c], 7)


def block(key, counter, nonce, rounds):
    init = SIGMA + key + [counter] + nonce
    s = list(init)
    for _ in range(rounds // 2):
        for a, b, c, d in [
            (0, 4, 8, 12),
            (1, 5, 9, 13),
            (2, 6, 10, 14),
            (3, 7, 11, 15),
        ]:
            quarter(s, a, b, c, d)
        for a, b, c, d in [
            (0, 5, 10, 15),
            (1, 6, 11, 12),
            (2, 7, 8, 13),
            (3, 4, 9, 14),
        ]:
            quarter(s, a, b, c, d)
    return [(x + y) & MASK for x, y in zip(s, init)]


def main():
    text = sys.stdin.read()
    words = lambda s: [int(w, 16) for w in s.split()]
    rfc_key = [
        int.from_bytes(bytes(range(4 * i, 4 * i + 4)), "little") for i in range(8)
    ]
    rfc = block(rfc_key, 1, [0x09000000, 0x4A000000, 0], 20)
    assert rfc == words(RFC), "python chacha20 misses the rfc 8439 2.3.2 vector"
    head = re.search(r"chacha8 key ([0-9a-f ]+), nonce ([0-9a-f ]+)", text)
    key, nonce = words(head.group(1)), words(head.group(2))
    rows = re.findall(r"chacha8 counter ([0-9a-f]+) out ([0-9a-f]+)", text)
    agree = 0
    for ctr, out in rows:
        b = block(key, int(ctr, 16), nonce, 8)
        agree += (b[0] | b[1] << 32) == int(out, 16)
    rust = re.search(r"chacha20 rfc 8439 2.3.2 block: ([0-9a-f ]+)", text)
    print(
        f"python chacha20 matches rfc 8439 2.3.2: True; rust block matches: {words(rust.group(1)) == rfc}"
    )
    print(
        f"python chacha8 agrees with the rust chacha8 on {agree} of {len(rows)} counters"
    )


main()

use std::collections::BTreeMap;

// SPANS

pub struct Span {
    pub bits: usize,
    pub dim: usize,
    pub gens: Vec<u32>,
    member: Vec<u64>,
    elems: Vec<u32>,
}

impl Span {
    pub fn new(bits: usize) -> Span {
        let mut member = vec![0u64; ((1usize << bits) + 63) / 64];
        member[0] = 1;
        Span {
            bits,
            dim: 0,
            gens: Vec::new(),
            member,
            elems: vec![0],
        }
    }

    pub fn has(&self, v: u32) -> bool {
        self.member[v as usize >> 6] >> (v & 63) & 1 == 1
    }

    pub fn insert(&mut self, v: u32) -> bool {
        if self.has(v) {
            return false;
        }
        for k in 0..self.elems.len() {
            let w = self.elems[k] ^ v;
            self.member[w as usize >> 6] |= 1 << (w & 63);
            self.elems.push(w);
        }
        self.gens.push(v);
        self.dim += 1;
        true
    }

    pub fn perp(&self) -> Vec<u32> {
        let mut out = Span::new(self.bits);
        for phi in 1..1u32 << self.bits {
            if self.gens.iter().all(|&g| (g & phi).count_ones() % 2 == 0) {
                out.insert(phi);
            }
        }
        out.gens
    }
}

pub struct Rng(pub u64);

impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

// BLOCKS

pub fn closure(gens: &[&[u32]], v: u32, bits: usize, full: &[bool], rng: &mut Rng) -> Option<Span> {
    let size = 1usize << bits;
    let mut u = Span::new(bits);
    u.insert(v);
    for _ in 0..4 * bits {
        let g = u.gens[rng.next() as usize % u.gens.len()] as usize;
        let x = rng.next() as usize % size;
        let h = gens[rng.next() as usize % gens.len()];
        let w = h[x ^ g] ^ h[x];
        if full[w as usize] {
            return None;
        }
        u.insert(w);
        if u.dim == bits {
            return None;
        }
    }
    let mut next = 0;
    while next < u.gens.len() {
        let g = u.gens[next] as usize;
        next += 1;
        for h in gens {
            for x in 0..size {
                let w = h[x ^ g] ^ h[x];
                if full[w as usize] || u.insert(w) && u.dim == bits {
                    return None;
                }
            }
        }
    }
    Some(u)
}

pub fn first_block(gens: &[&[u32]], bits: usize) -> Option<Span> {
    let size = 1usize << bits;
    let mut full = vec![false; size];
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    for v in 1..size {
        match closure(gens, v as u32, bits, &full, &mut rng) {
            None => full[v] = true,
            Some(u) => return Some(u),
        }
    }
    None
}

pub fn is_block(gens: &[&[u32]], u: &Span) -> bool {
    let size = 1usize << u.bits;
    gens.iter().all(|h| {
        u.gens
            .iter()
            .all(|&g| (0..size).all(|x| u.has(h[x ^ g as usize] ^ h[x])))
    })
}

// CYCLES

pub fn cycles(perm: &[u32]) -> Vec<usize> {
    let mut seen = vec![false; perm.len()];
    let mut out = Vec::new();
    for s in 0..perm.len() {
        if seen[s] {
            continue;
        }
        let (mut x, mut len) = (s, 0);
        while !seen[x] {
            seen[x] = true;
            x = perm[x] as usize;
            len += 1;
        }
        out.push(len);
    }
    out
}

pub fn is_prime(p: usize) -> bool {
    p >= 2 && (2..).take_while(|d| d * d <= p).all(|d| p % d != 0)
}

pub fn long_prime(lengths: &[usize], size: usize) -> Option<usize> {
    lengths
        .iter()
        .copied()
        .find(|&l| 2 * l > size && l + 3 <= size && is_prime(l))
}

pub fn lone_prime(lengths: &[usize], size: usize) -> Option<usize> {
    let mut count: BTreeMap<usize, usize> = BTreeMap::new();
    for &l in lengths {
        *count.entry(l).or_default() += 1;
    }
    count
        .iter()
        .filter(|&(&l, &c)| c == 1 && is_prime(l) && l + 3 <= size)
        .map(|(&l, _)| l)
        .find(|&p| lengths.iter().filter(|&&l| l % p == 0).count() == 1)
}

pub fn wernsdorf(gens: &[&[u32]], bits: usize, tries: usize, seed: u64) -> Option<usize> {
    let size = 1usize << bits;
    let mut rng = Rng(seed);
    for t in 0..tries {
        let len = 1 + t % 4;
        let mut w: Vec<u32> = (0..size as u32).collect();
        for _ in 0..len {
            for h in gens {
                let c = rng.next() as u32 % size as u32;
                for y in w.iter_mut() {
                    *y = h[*y as usize] ^ c;
                }
            }
        }
        if let Some(p) = long_prime(&cycles(&w), size) {
            return Some(p);
        }
    }
    None
}

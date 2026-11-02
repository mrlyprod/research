use mrly_core::colors::{BLUE, LIGHT, ORANGE, WHITE};
use png::{AdaptiveFilterType, BitDepth, ColorType, Compression, Encoder, FilterType};
use std::path::Path;

pub const WIDTH: usize = 1760;
pub const HEIGHT: usize = 736;
const SURF: [u8; 4] = [WHITE.r, WHITE.g, WHITE.b, 255];
const NAVY: [u8; 4] = [BLUE.r, BLUE.g, BLUE.b, 255];
const RUST: [u8; 4] = [ORANGE.r, ORANGE.g, ORANGE.b, 255];
const GREY: [u8; 4] = [LIGHT.dim.r, LIGHT.dim.g, LIGHT.dim.b, 255];
const FAINT: [u8; 4] = [LIGHT.line.r, LIGHT.line.g, LIGHT.line.b, 255];
const INK: [u8; 4] = [LIGHT.fg.r, LIGHT.fg.g, LIGHT.fg.b, 255];

pub struct Series {
    pub spectral: Vec<f64>,
    pub walker: Vec<f64>,
    pub fractal: Vec<f64>,
}

struct Frame {
    left: f64,
    right: f64,
    top: f64,
    bottom: f64,
    x: (f64, f64),
    y: (f64, f64),
}

impl Frame {
    fn place(&self, x: f64, y: f64) -> (f64, f64) {
        let across = (x - self.x.0) / (self.x.1 - self.x.0);
        let up = (y - self.y.0) / (self.y.1 - self.y.0);
        (
            self.left + across * (self.right - self.left),
            self.bottom - up * (self.bottom - self.top),
        )
    }
}

struct Canvas {
    pixels: Vec<[u8; 4]>,
}

impl Canvas {
    fn blend(&mut self, x: f64, y: f64, colour: [u8; 4], alpha: f64) {
        if alpha <= 0.0 || x < 0.0 || y < 0.0 || x >= WIDTH as f64 || y >= HEIGHT as f64 {
            return;
        }
        let seat = y as usize * WIDTH + x as usize;
        for (pixel, &paint) in self.pixels[seat].iter_mut().zip(&colour).take(3) {
            let under = *pixel as f64;
            *pixel = (under + (paint as f64 - under) * alpha).round() as u8;
        }
    }

    fn disc(&mut self, cx: f64, cy: f64, radius: f64, colour: [u8; 4]) {
        let reach = radius + 2.0;
        let mut y = (cy - reach).max(0.0).floor();
        while y <= cy + reach {
            let mut x = (cx - reach).max(0.0).floor();
            while x <= cx + reach {
                let away = (x + 0.5 - cx).hypot(y + 0.5 - cy);
                self.blend(x, y, colour, (radius + 0.5 - away).clamp(0.0, 1.0));
                x += 1.0;
            }
            y += 1.0;
        }
    }

    fn square(&mut self, cx: f64, cy: f64, half: f64, colour: [u8; 4]) {
        let mut y = (cy - half).max(0.0).floor();
        while y <= cy + half {
            let mut x = (cx - half).max(0.0).floor();
            while x <= cx + half {
                self.blend(x, y, colour, 1.0);
                x += 1.0;
            }
            y += 1.0;
        }
    }

    fn line(&mut self, x0: f64, y0: f64, x1: f64, y1: f64, colour: [u8; 4], thick: f64) {
        let span = (x1 - x0).hypot(y1 - y0).max(1.0);
        let steps = (span * 2.0).ceil() as usize;
        for step in 0..=steps {
            let part = step as f64 / steps as f64;
            let x = x0 + (x1 - x0) * part;
            let y = y0 + (y1 - y0) * part;
            self.square(x, y, thick / 2.0, colour);
        }
    }

    fn star(&mut self, cx: f64, cy: f64, radius: f64, colour: [u8; 4]) {
        for spike in 0..5 {
            let angle = -std::f64::consts::FRAC_PI_2 + spike as f64 * std::f64::consts::TAU / 5.0;
            self.line(
                cx,
                cy,
                cx + radius * angle.cos(),
                cy + radius * angle.sin(),
                colour,
                2.0,
            );
        }
    }

    fn axes(&mut self, frame: &Frame, along: &[f64], up: &[f64]) {
        self.line(
            frame.left,
            frame.bottom,
            frame.right,
            frame.bottom,
            INK,
            1.5,
        );
        self.line(frame.left, frame.top, frame.left, frame.bottom, INK, 1.5);
        for value in along {
            let (x, y) = frame.place(*value, frame.y.0);
            self.line(x, y, x, y + 6.0, INK, 1.5);
        }
        for value in up {
            let (x, y) = frame.place(frame.x.0, *value);
            self.line(x - 6.0, y, x, y, INK, 1.5);
            self.line(x, y, frame.right, y, FAINT, 1.0);
        }
    }

    fn bracket(&mut self, frame: &Frame, from: f64, to: f64, height: f64) {
        let (x0, y0) = frame.place(from, height);
        let (x1, _) = frame.place(to, height);
        self.line(x0, y0, x1, y0, INK, 1.0);
        self.line(x0, y0, x0, y0 + 7.0, INK, 1.0);
        self.line(x1, y0, x1, y0 + 7.0, INK, 1.0);
    }
}

pub fn render(series: &Series) -> Vec<[u8; 4]> {
    let mut canvas = Canvas {
        pixels: vec![SURF; WIDTH * HEIGHT],
    };
    let left = Frame {
        left: 110.0,
        right: 830.0,
        top: 70.0,
        bottom: 650.0,
        x: (-0.5, 8.5),
        y: (1.9, 2.9),
    };
    let right = Frame {
        left: 980.0,
        right: 1700.0,
        top: 70.0,
        bottom: 650.0,
        x: (1.3, 2.9),
        y: (0.9, 2.7),
    };
    let ticks: Vec<f64> = (0..9).map(|index| index as f64).collect();
    canvas.axes(&left, &ticks, &[1.9, 2.1, 2.3, 2.5, 2.7, 2.9]);
    let (x0, y0) = left.place(left.x.0, 2.0);
    let (x1, _) = left.place(left.x.1, 2.0);
    canvas.line(x0, y0, x1, y0, GREY, 1.5);
    canvas.bracket(&left, 2.0, 3.0, 2.71);
    canvas.bracket(&left, 4.0, 5.0, 2.28);
    for (index, value) in series.spectral.iter().enumerate() {
        let (x, y) = left.place(index as f64, *value);
        canvas.disc(x, y, 7.0, NAVY);
    }
    for (index, value) in series.walker.iter().enumerate() {
        let (x, y) = left.place(index as f64, *value);
        canvas.square(x, y, 5.0, RUST);
    }
    canvas.axes(
        &right,
        &[1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8],
        &[1.0, 1.4, 1.8, 2.2, 2.6],
    );
    let start = right.x.0.max(right.y.0);
    let stop = right.x.1.min(right.y.1);
    let (dx0, dy0) = right.place(start, start);
    let (dx1, dy1) = right.place(stop, stop);
    canvas.line(dx0, dy0, dx1, dy1, GREY, 1.5);
    for (df, dw) in series.fractal.iter().zip(&series.walker) {
        let (x, y) = right.place(*df, 2.0 * df / dw);
        canvas.disc(x, y, 7.0, NAVY);
    }
    let gasket = right.place(3f64.ln() / 2f64.ln(), 2.0 * 3f64.ln() / 5f64.ln());
    canvas.star(gasket.0, gasket.1, 13.0, RUST);
    canvas.pixels
}

pub fn write(path: &Path, series: &Series) {
    let bytes = png(&render(series));
    std::fs::write(path, &bytes).expect("the figure writes");
}

// PNG

fn table(pixels: &[[u8; 4]]) -> Option<Vec<[u8; 4]>> {
    let mut keys: Vec<u32> = Vec::with_capacity(256);
    for pixel in pixels {
        let key = u32::from_be_bytes(*pixel);
        if let Err(slot) = keys.binary_search(&key) {
            if keys.len() == 256 {
                return None;
            }
            keys.insert(slot, key);
        }
    }
    let mut table: Vec<[u8; 4]> = keys.into_iter().map(u32::to_be_bytes).collect();
    table.sort_by_key(|color| color[3]);
    Some(table)
}

fn depth(count: usize) -> BitDepth {
    match count {
        0..=2 => BitDepth::One,
        3..=4 => BitDepth::Two,
        5..=16 => BitDepth::Four,
        _ => BitDepth::Eight,
    }
}

fn plte(table: &[[u8; 4]]) -> Vec<u8> {
    table
        .iter()
        .flat_map(|color| [color[0], color[1], color[2]])
        .collect()
}

fn trns(table: &[[u8; 4]]) -> Vec<u8> {
    let opaque = table
        .iter()
        .position(|color| color[3] == 255)
        .unwrap_or(table.len());
    table[..opaque].iter().map(|color| color[3]).collect()
}

fn indices(
    pixels: &[[u8; 4]],
    table: &[[u8; 4]],
    width: usize,
    height: usize,
    depth: BitDepth,
) -> Vec<u8> {
    let mut lookup: Vec<(u32, u8)> = table
        .iter()
        .enumerate()
        .map(|(slot, color)| (u32::from_be_bytes(*color), slot as u8))
        .collect();
    lookup.sort_unstable();
    let bits = depth as usize;
    let per_byte = 8 / bits;
    let stride = width.div_ceil(per_byte);
    let mut data = vec![0u8; stride * height];
    if width == 0 {
        return data;
    }
    for (y, row) in pixels.chunks_exact(width).enumerate() {
        for (x, pixel) in row.iter().enumerate() {
            let key = u32::from_be_bytes(*pixel);
            let slot = lookup.binary_search_by_key(&key, |entry| entry.0).unwrap();
            let shift = 8 - bits * (x % per_byte + 1);
            data[y * stride + x / per_byte] |= lookup[slot].1 << shift;
        }
    }
    data
}

fn png(pixels: &[[u8; 4]]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(pixels.len() + 128);
    {
        let mut encoder = Encoder::new(&mut bytes, WIDTH as u32, HEIGHT as u32);
        encoder.set_compression(Compression::Best);
        encoder.set_filter(FilterType::Paeth);
        encoder.set_adaptive_filter(AdaptiveFilterType::Adaptive);
        let data = match table(pixels) {
            Some(table) => {
                let depth = depth(table.len());
                encoder.set_color(ColorType::Indexed);
                encoder.set_depth(depth);
                encoder.set_palette(plte(&table));
                let veils = trns(&table);
                if !veils.is_empty() {
                    encoder.set_trns(veils);
                }
                indices(pixels, &table, WIDTH, HEIGHT, depth)
            }
            None => {
                encoder.set_color(ColorType::Rgba);
                encoder.set_depth(BitDepth::Eight);
                pixels.concat()
            }
        };
        let mut writer = encoder.write_header().expect("the figure encodes");
        writer.write_image_data(&data).expect("the figure encodes");
        writer.finish().expect("the figure encodes");
    }
    bytes
}

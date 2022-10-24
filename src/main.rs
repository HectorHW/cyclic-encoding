use std::{
    fmt::Write,
    iter::{once, repeat},
    ops::Mul,
};

use algebraics::polynomial::Polynomial;
use bits::BitVector;
use colored::*;
use encoding_rs::WINDOWS_1251;

use itertools::Itertools;

use crate::bits::break_one_bit;

mod bits;

type Pol = Polynomial<i32>;

fn into_cp1251_bits(letter: char) -> BitVector<8> {
    let string_buf = String::from(letter);
    let buf = WINDOWS_1251.encode(&string_buf).0.to_vec();
    debug_assert_eq!(buf.len(), 1);
    buf[0].into()
}

fn with_highlight(vec: &[bool], positions: &[usize]) -> String {
    let mut buf = String::new();
    for (i, digit) in vec.iter().cloned().enumerate() {
        if positions.contains(&i) {
            write!(buf, "{}", (digit as u8).to_string().red()).unwrap();
        } else {
            write!(buf, "{}", (digit as u8).to_string().white()).unwrap();
        }
    }
    buf
}

fn into_polynomial<const N: usize>(bits: &BitVector<N>) -> Pol {
    bits.iter().cloned().map(i32::from).collect()
}

fn into_bitvector<const N: usize>(p: Pol) -> BitVector<N> {
    p.into_coefficients()
        .into_iter()
        .chain(repeat(0))
        .take(N)
        .collect_vec()
        .try_into()
        .unwrap()
}

fn shift_polynomial(pol: &Pol, shift: usize) -> Pol {
    pol.mul(Polynomial::from(
        repeat(0).take(shift).chain(once(1)).collect::<Vec<_>>(),
    ))
}

fn wrap_mod_2(p: Pol) -> Pol {
    let mut coefficients: Vec<_> = p.into_coefficients();
    coefficients.iter_mut().for_each(|item| {
        *item = item.rem_euclid(2);
    });
    while let Some(0) = coefficients.last() {
        coefficients.pop();
    }
    coefficients.into()
}

fn add_mod2(x: &Pol, y: &Pol) -> Pol {
    wrap_mod_2(x + y)
}

fn div_mod2(x: Pol, y: &Pol) -> (Pol, Pol) {
    let (d, r) = x.div_rem(y);
    (wrap_mod_2(d), wrap_mod_2(r))
}

#[cfg(feature = "verbose")]
fn pretty_print_poly(p: &Pol) -> String {
    let coeffincients = p.clone().into_coefficients();
    if coeffincients.is_empty() || coeffincients == vec![0] {
        return "0".to_string();
    }
    coeffincients
        .into_iter()
        .enumerate()
        .fold(String::new(), |mut a, (power, coeff)| {
            if coeff == 0 {
                return a;
            }

            if a.is_empty() {
                write!(
                    a,
                    "{}",
                    match power {
                        0 => "1".to_string(),
                        1 => "x".to_string(),
                        n => format!("x^{n}"),
                    }
                )
                .unwrap();
                a
            } else {
                write!(
                    a,
                    " + {}",
                    match power {
                        1 => "x".to_string(),
                        n => format!("x^{n}"),
                    }
                )
                .unwrap();
                a
            }
        })
}

fn encode<const N: usize, const K: usize>(bits: BitVector<K>, base: &Pol) -> BitVector<N> {
    #[cfg(feature = "verbose")]
    println!("bits: {bits}");
    let p = into_polynomial(&bits);
    #[cfg(feature = "verbose")]
    println!("polynomial: {}", pretty_print_poly(&p));
    let shifted = shift_polynomial(&p, N - K);
    #[cfg(feature = "verbose")]
    println!("shifted polynomial: {}", pretty_print_poly(&shifted));
    let (_, rem) = div_mod2(shifted.clone(), base);
    #[cfg(feature = "verbose")]
    println!("remainder: {}", pretty_print_poly(&rem));
    let result = add_mod2(&rem, &shifted);
    #[cfg(feature = "verbose")]
    println!("shifted + remainder: {}", pretty_print_poly(&result));
    let bits: BitVector<N> = into_bitvector(result);
    #[cfg(feature = "verbose")]
    println!("bits: {bits}");
    bits
}

fn syndrome<const N: usize, const N_K: usize>(data: &BitVector<N>, base: &Pol) -> BitVector<N_K> {
    into_bitvector(div_mod2(into_polynomial(data), base).1)
}

#[cfg(feature = "compute-distance")]
fn distance<const N: usize>(v1: &BitVector<N>, v2: &BitVector<N>) -> usize {
    v1.iter().zip(v2.iter()).filter(|(x, y)| x != y).count()
}

#[cfg(feature = "compute-distance")]
fn compute_distance(g: &Pol) -> (usize, BitVector<8>, BitVector<8>) {
    let mut min_distance = 13;
    let mut examples = (None, None);

    for w1 in 0..255u8 {
        for w2 in 0..w1 {
            let (w1, w2) = (BitVector::from(w1), BitVector::from(w2));
            let (c1, c2) = (
                encode::<13, 8>(w1.clone(), g),
                encode::<13, 8>(w2.clone(), g),
            );
            if distance(&c1, &c2) < min_distance {
                min_distance = distance(&c1, &c2);
                let _ = examples.0.replace(w1);
                examples.1.replace(w2);
            }
        }
    }
    (min_distance, examples.0.unwrap(), examples.1.unwrap())
}

fn main() {
    // 5 вариант - x^5 + x^4 + x^3 + x + 1

    const N: usize = 13;
    const K: usize = 8;

    let g = Polynomial::from(vec![1, 1, 0, 1, 1, 1]);

    assert_eq!(g.degree(), Some(N - K));

    let alphabet = ['л', 'м', 'н', 'о', 'п'];

    let mut codes = vec![];

    for (code, &letter) in alphabet
        .iter()
        .cloned()
        .map(into_cp1251_bits)
        .zip(alphabet.iter())
    {
        let encoded = encode::<N, K>(code.clone(), &g);
        println!("{letter} => {} => {}", code, encoded);
        codes.push(encoded);
    }

    for encoded in &codes {
        let s = syndrome::<N, { N - K }>(encoded, &g);
        println!("encoded: {}, syndrome: {}", encoded, s);
    }

    for encoded in &codes {
        let mut with_bit_broken = encoded.clone();
        let idx = break_one_bit(&mut with_bit_broken);
        let s = syndrome::<N, { N - K }>(&with_bit_broken, &g);
        println!(
            "{} syndrome: {}",
            with_highlight(&with_bit_broken, &[idx]),
            s
        );
    }

    #[cfg(feature = "compute-distance")]
    {
        let (dist, w1, w2) = compute_distance(&g);
        let (c1, c2) = (
            encode::<13, 8>(w1.clone(), &g),
            encode::<13, 8>(w2.clone(), &g),
        );
        println!(
            "min distance: {} for codes: {} and {}\nThese words correspond to {} and {}",
            dist, c1, c2, w1, w2
        )
    }
}

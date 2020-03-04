use algebra::{
    fields::{BitIterator, Field},
    prelude::{One, Zero},
    TEModelParameters,
};
use blake2_rfc::blake2b::Blake2b;
use core::{mem, ops::AddAssign};

fn hash_to_scalar<E>(persona: &[u8], a: &[u8], b: &[u8]) -> E::ScalarField
where
    E: TEModelParameters,
    E::ScalarField: One + Zero + Field,
{
    let mut hasher = Blake2b::with_params(64, &[], &[], persona);
    hasher.update(a);
    hasher.update(b);
    let res = hasher.finalize();
    assert_eq!(res.len(), 64);
    let mut slice = res.as_ref();
    let mut repr: [u64; 8] = [0; 8];
    for r in &mut repr {
        let (value, rest) = slice.split_at(mem::size_of::<u64>());
        let mut v = [0u8; 8];
        v.copy_from_slice(value);
        *r = u64::from_le_bytes(v);
        slice = rest;
    }
    let bits = BitIterator::new(repr);

    let one = E::ScalarField::one();

    // mul_bits
    let mut res = E::ScalarField::zero();
    for bit in bits {
        res.double_in_place();

        if bit {
            res.add_assign(one);
        }
    }

    res
}

pub fn h_star<E>(a: &[u8], b: &[u8]) -> E::ScalarField
where
    E: TEModelParameters,
    E::ScalarField: One + Zero + Field,
{
    hash_to_scalar::<E>(b"Zcash_RedJubjubH", a, b)
}

#[cfg(test)]
mod tests {
    use super::h_star;
    use algebra::{
        biginteger::BigInteger256 as BigInteger, curves::jubjub::JubJubParameters, fields::Fp256,
    };

    #[test]
    fn test_h_start_empty() {
        let a = h_star::<JubJubParameters>(b"", b"");
        let expected = Fp256::new(BigInteger([
            0x19fe572ca31d6630,
            0xd10a239b34282920,
            0x5205644c4cd0968f,
            0x06804dcb18325f7c,
        ]));

        assert_eq!(a, expected);
    }

    #[test]
    fn test_h_start_with_params() {
        let a = h_star::<JubJubParameters>(b"ab", b"cde");
        let expected = Fp256::new(BigInteger([
            0x39caf2bba0b7ad6c,
            0x2176ca4852d69ccf,
            0xafb704c4c93bd0db,
            0x0a322f2d738215b2,
        ]));

        assert_eq!(a, expected);
    }
}

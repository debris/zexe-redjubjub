use crate::{
    constants,
    point::{mul_by_cofactor, read_point},
    Point,
};
use algebra::{biginteger::BigInteger256, prelude::Zero, PrimeField, TEModelParameters};
use blake2_rfc::blake2s::Blake2s;

/// Produces a random point in the Jubjub curve.
/// The point is guaranteed to be prime order
/// and not the identity.
pub fn group_hash<E>(tag: &[u8], personalization: &[u8]) -> Option<Point<E>>
where
    E: TEModelParameters,
    E::BaseField: PrimeField + Into<BigInteger256>,
{
    assert_eq!(personalization.len(), 8);

    let mut h = Blake2s::with_params(32, &[], &[], personalization);
    h.update(constants::GH_FIRST_BLOCK);
    h.update(tag);
    let h = h.finalize().as_ref().to_vec();
    assert!(h.len() == 32);

    let mut p = read_point(&h[..])?;

    p = mul_by_cofactor(&p);

    if !p.is_zero() {
        Some(p)
    } else {
        None
    }
}

pub fn find_group_hash<E>(m: &[u8], personalization: &[u8; 8]) -> Point<E>
where
    E: TEModelParameters,
    E::BaseField: PrimeField + Into<BigInteger256>,
{
    let mut tag = m.to_vec();
    let i = tag.len();
    tag.push(0u8);

    loop {
        let gh = group_hash(&tag, personalization);

        // We don't want to overflow and start reusing generators
        assert!(tag[i] != u8::max_value());
        tag[i] += 1;

        if let Some(gh) = gh {
            break gh;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::group_hash;
    use crate::Point;
    use algebra::{
        biginteger::BigInteger256 as BigInteger, curves::jubjub::JubJubParameters, fields::Fp256,
    };

    #[test]
    fn test_group_hash_x() {
        let p = group_hash::<JubJubParameters>(b"abc", b"01234567");

        let x = Fp256::new(BigInteger([
            0x955e8ca79ad70f31,
            0x8bc0225726661e50,
            0xf80045d17b2b5eb8,
            0x02a80185cb9b5fea,
        ]));

        let y = Fp256::new(BigInteger([
            0x65673f3d5ee7e357,
            0xefdc12f784983858,
            0xaa4d01adc4bd1d48,
            0x7356a532a4b76695,
        ]));

        let t = Fp256::new(BigInteger([
            0xb80ffc432eb63a7a,
            0xffdb920ea6d25bc1,
            0xee076e4e5564e7e3,
            0x116f5a0ea64bf238,
        ]));

        let z = Fp256::new(BigInteger([
            0xdcd2a9b715785d16,
            0xb629a1a1c727895f,
            0x5fe596f1edd96753,
            0x2e46f87723ecca72,
        ]));

        let expected = Point::new(x, y, t, z);
        assert_eq!(p, Some(expected));
    }
}

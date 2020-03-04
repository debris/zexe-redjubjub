use crate::{constants, Point};
use algebra::{
    biginteger::{BigInteger, BigInteger256},
    bytes::FromBytes,
    fields::{Field, SquareRootField},
    prelude::{One, Zero},
    PrimeField, TEModelParameters,
};
use blake2_rfc::blake2s::Blake2s;
use core::ops::{AddAssign, MulAssign, Neg, SubAssign};

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

    let mut y_repr = <E::BaseField as PrimeField>::BigInt::read(&h[..]).ok()?;

    let x_sign = (y_repr.as_ref()[3] >> 63) == 1;
    y_repr.as_mut()[3] &= 0x7fffffffffffffff;

    let y = E::BaseField::from_repr(y_repr);
    let mut p = get_for_y(y, x_sign)?;

    // mul_by_cofactor
    p = double(&p);
    p = double(&p);
    p = double(&p);

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

fn get_for_y<E>(y: E::BaseField, sign: bool) -> Option<Point<E>>
where
    E: TEModelParameters,
    E::BaseField: Into<BigInteger256>,
{
    // Given a y on the curve, x^2 = (y^2 - 1) / (dy^2 + 1)
    // This is defined for all valid y-coordinates,
    // as dy^2 + 1 = 0 has no solution in Fr.

    // tmp1 = y^2
    let mut tmp1 = y;
    tmp1.square_in_place();

    // tmp2 = (y^2 * d) + 1
    let mut tmp2 = tmp1;
    //tmp2.mul_assign(params.edwards_d());
    tmp2.mul_assign(&E::COEFF_D);
    tmp2.add_assign(&E::BaseField::one());

    // tmp1 = y^2 - 1
    tmp1.sub_assign(&E::BaseField::one());

    match tmp2.inverse() {
        Some(tmp2) => {
            // tmp1 = (y^2 - 1) / (dy^2 + 1)
            tmp1.mul_assign(&tmp2);

            match tmp1.sqrt() {
                Some(mut x) => {
                    let bi: BigInteger256 = x.into();
                    if bi.is_odd() != sign {
                        x = x.neg();
                    }

                    let mut t = x;
                    t.mul_assign(&y);

                    Some(Point::new(x, y, t, E::BaseField::one()))
                }
                None => None,
            }
        }
        None => None,
    }
}

fn double<E: TEModelParameters>(point: &Point<E>) -> Point<E> {
    // See "Twisted Edwards Curves Revisited"
    //     Huseyin Hisil, Kenneth Koon-Ho Wong, Gary Carter, and Ed Dawson
    //     Section 3.3
    //     http://hyperelliptic.org/EFD/g1p/auto-twisted-extended.html#doubling-dbl-2008-hwcd

    // A = X1^2
    let mut a = point.x;
    a.square_in_place();

    // B = Y1^2
    let mut b = point.y;
    b.square_in_place();

    // C = 2*Z1^2
    let mut c = point.z;
    c.square_in_place();
    c.double_in_place();

    // D = a*A
    //   = -A
    let mut d = a;
    d = d.neg();

    // E = (X1+Y1)^2 - A - B
    let mut e = point.x;
    e.add_assign(&point.y);
    e.square_in_place();
    e.add_assign(&d); // -A = D
    e.sub_assign(&b);

    // G = D+B
    let mut g = d;
    g.add_assign(&b);

    // F = G-C
    let mut f = g;
    f.sub_assign(&c);

    // H = D-B
    let mut h = d;
    h.sub_assign(&b);

    // X3 = E*F
    let mut x3 = e;
    x3.mul_assign(&f);

    // Y3 = G*H
    let mut y3 = g;
    y3.mul_assign(&h);

    // T3 = E*H
    let mut t3 = e;
    t3.mul_assign(&h);

    // Z3 = F*G
    let mut z3 = f;
    z3.mul_assign(&g);

    Point::new(x3, y3, t3, z3)
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

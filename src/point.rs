//! Redjubjub specific Point utilities

use algebra::{
    biginteger::{BigInteger, BigInteger256},
    bytes::{FromBytes, ToBytes},
    curves::{twisted_edwards_extended::GroupProjective, ProjectiveCurve},
    fields::{Field, SquareRootField},
    io::{self, Read, Write},
    prelude::{One, Zero},
    PrimeField, TEModelParameters,
};
use core::ops::{AddAssign, MulAssign, Neg, SubAssign};

pub type Point<E> = GroupProjective<E>;

/// TODO: return algebra::io::Error
pub fn read_point<E, R>(read: R) -> Option<Point<E>>
where
    E: TEModelParameters,
    E::BaseField: PrimeField + Into<BigInteger256>,
    R: Read,
{
    let mut y_repr = <E::BaseField as PrimeField>::BigInt::read(read).ok()?;

    let x_sign = (y_repr.as_ref()[3] >> 63) == 1;
    y_repr.as_mut()[3] &= 0x7fffffffffffffff;

    let y = E::BaseField::from_repr(y_repr);
    let mut p = get_for_y(y, x_sign)?;

    p = mul_by_cofactor(&p);

    if !p.is_zero() {
        Some(p)
    } else {
        None
    }
}

pub fn write_point<E, W>(point: &Point<E>, writer: W) -> io::Result<()>
where
    E: TEModelParameters,
    E::BaseField: PrimeField,
    W: Write,
{
    let affine = point.into_affine();
    let x_repr = affine.x.into_repr();
    let mut y_repr = affine.y.into_repr();

    if x_repr.is_odd() {
        y_repr.as_mut()[3] |= 0x8000000000000000u64;
    }

    y_repr.write(writer)
}

pub fn mul_by_cofactor<E: TEModelParameters>(p: &Point<E>) -> Point<E> {
    double(&double(&double(p)))
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

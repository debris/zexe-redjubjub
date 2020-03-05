#![no_std]

mod constants;
mod generator;
mod group_hash;
mod point;
mod util;

use algebra::{
    biginteger::BigInteger256,
    bytes::{FromBytes, ToBytes},
    io::{self, Read},
    prelude::Zero,
    PrimeField, TEModelParameters,
};
use core::{
    fmt,
    ops::{Add, AddAssign, Mul, MulAssign, Neg},
};
use point::mul_by_cofactor;
use rand::Rng;
use util::h_star;

pub use generator::FixedGenerators;
pub use point::{read_point, write_point, Point};

pub struct PrivateKey<E: TEModelParameters> {
    pub field: E::ScalarField,
}

impl<E: TEModelParameters> fmt::Debug for PrivateKey<E> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("PrivateKey")
            .field("field", &self.field)
            .finish()
    }
}

impl<E> PrivateKey<E>
where
    E: TEModelParameters,
    E::BaseField: PrimeField + Into<BigInteger256>,
{
    pub fn sign<R: Rng>(&self, msg: &[u8], rng: &mut R, generator: FixedGenerators) -> Signature {
        // T = (l_H + 128) bits of randomness
        // For H*, l_H = 512 bits
        let mut t = [0u8; 80];
        rng.fill_bytes(&mut t[..]);

        // r = H*(T || M)
        let r = h_star::<E>(&t[..], msg);

        // R = r . P_G
        //let r_g = params.generator(p_g).mul(r, params);
        let r_g = generator.point::<E>().mul(&r);
        let mut rbar = [0u8; 32];
        write_point(&r_g, &mut rbar[..]).expect("Jubjub points should serialize to 32 bytes");

        // S = r + H*(Rbar || M) . sk
        let mut s = h_star::<E>(&rbar[..], msg);
        s.mul_assign(&self.field);
        s.add_assign(&r);
        let mut sbar = [0u8; 32];
        s.write(&mut sbar[..])
            .expect("Jubjub scalars should serialize to 32 bytes");

        Signature { rbar, sbar }
    }
}

pub struct PublicKey<E: TEModelParameters> {
    pub point: Point<E>,
}

impl<E> FromBytes for PublicKey<E> 
where
    E: TEModelParameters,
    E::BaseField: PrimeField + Into<BigInteger256>,
{
    #[inline]
    fn read<R: Read>(read: R) -> io::Result<Self> {
        // TODO: handle error
        let point = read_point::<E, R>(read).unwrap();
        Ok(PublicKey { point })
    }
}

impl<E> PublicKey<E>
where
    E: TEModelParameters,
    E::BaseField: PrimeField + Into<BigInteger256>,
{
    pub fn new(point: Point<E>) -> Self {
        PublicKey {
            point
        }
    }

    pub fn from_private(privkey: &PrivateKey<E>, generator: FixedGenerators) -> Self {
        PublicKey {
            point: generator.point().mul(&privkey.field),
        }
    }

    pub fn verify(&self, msg: &[u8], sig: &Signature, generator: FixedGenerators) -> bool {
        // c = H*(Rbar || M)
        let c = h_star::<E>(&sig.rbar[..], msg);

        // Signature checks:
        // R != invalid
        let r = match read_point(&sig.rbar[..]) {
            Some(r) => r,
            None => return false,
        };

        // S < order(G)
        // (E::ScalarField guarantees its representation is in the field)
        let s = match E::ScalarField::read(&sig.sbar[..]) {
            Ok(s) => s,
            Err(_) => return false,
        };

        // 0 = h_G(-S . P_G + R + c . vk)
        let p = self
            .point
            .mul(&c)
            .add(&r)
            .add(&generator.point().mul(&s).neg());

        mul_by_cofactor(&p).is_zero()
    }
}

#[derive(Debug, PartialEq)]
pub struct Signature {
    rbar: [u8; 32],
    sbar: [u8; 32],
}

impl FromBytes for Signature {
    fn read<R: Read>(mut reader: R) -> io::Result<Self> {
        let mut rbar = [0u8; 32];
        let mut sbar = [0u8; 32];
        reader.read_exact(&mut rbar)?;
        reader.read_exact(&mut sbar)?;
        Ok(Signature { rbar, sbar })
    }
}

#[cfg(test)]
mod tests {
    use super::{FixedGenerators, PrivateKey, PublicKey, Signature};
    use algebra::curves::jubjub::JubJubParameters;
    use rand::{rngs::mock::StepRng, Rng};

    #[test]
    fn sign_and_verify() {
        let mut rng = StepRng::new(0, 1);
        let generator = FixedGenerators::SpendingKeyGenerator;
        let privkey: PrivateKey<JubJubParameters> = PrivateKey { field: rng.gen() };

        let msg1 = b"Foo bar";
        let sig1 = privkey.sign(msg1, &mut rng, generator);

        let expected = Signature {
            rbar: *b"\
                \xfc\x2a\xea\x10\xcf\xd3\x3d\x75\xbb\x09\x05\xf4\xa7\x08\xb6\x82\
                \xa6\xa0\x33\x80\x18\x53\x0f\x95\x84\xad\x59\x66\x30\x41\x81\x0a\
            ",
            sbar: *b"\
                \xd1\xed\xf2\xf6\xd3\x22\x7b\x49\x83\xe5\x77\x94\x92\x79\x33\x34\
                \x84\xa9\x3a\x35\x58\x7f\xbd\x78\xa5\xe1\x95\x41\x75\x55\x67\x0e\
            ",
        };

        assert_eq!(expected, sig1);
        let pubkey = PublicKey::from_private(&privkey, generator);
        assert!(pubkey.verify(msg1, &sig1, generator));
    }
}

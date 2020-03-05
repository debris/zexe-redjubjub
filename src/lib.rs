//#![no_std]

mod constants;
mod generator;
mod group_hash;
mod point;
mod util;

use algebra::{
    biginteger::BigInteger256,
    bytes::{FromBytes, ToBytes},
    curves::{
        edwards_bls12::EdwardsParameters,
        models::twisted_edwards_extended::{GroupAffine, GroupProjective},
    },
    io::{self, Read},
    prelude::Zero,
    PrimeField, TEModelParameters,
};
use point::mul_by_cofactor;
use rand::Rng;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg};
use util::h_star;

pub use generator::FixedGenerators;
pub use point::{read_point, write_point, Point};

pub struct PrivateKey<E: TEModelParameters> {
    pub field: E::ScalarField,
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

impl<E: TEModelParameters> FromBytes for PublicKey<E> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> io::Result<Self> {
        let point = Point::read(&mut reader)?;
        Ok(PublicKey { point })
    }
}

impl<E> PublicKey<E>
where
    E: TEModelParameters,
    E::BaseField: PrimeField + Into<BigInteger256>,
{
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
    use super::{FixedGenerators, Point, PrivateKey, PublicKey, Signature};
    use algebra::{
        biginteger::BigInteger256 as BigInteger, curves::jubjub::JubJubParameters, fields::Fp256,
        test_rng,
    };
    use rand::Rng;

    #[test]
    fn public_key_verify() {
        let mut rng = test_rng();
        let generator = FixedGenerators::SpendingKeyGenerator;
        let privkey: PrivateKey<JubJubParameters> = PrivateKey { field: rng.gen() };

        let msg1 = b"Foo bar";
        let sig1 = privkey.sign(msg1, &mut rng, generator);
        let pubkey = PublicKey::from_private(&privkey, generator);
        assert!(pubkey.verify(msg1, &sig1, generator));
    }
}

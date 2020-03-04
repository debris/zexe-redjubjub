//#![no_std]

mod util;

use std::ops::{Mul, Add};
use algebra::{
    bytes::FromBytes,
    io::{Read, self},
    curves::{
        edwards_bls12::EdwardsParameters,
        models::{
            TEModelParameters,
            twisted_edwards_extended::{GroupAffine, GroupProjective},
        },
    },
};
use util::h_star;


type Point<E: TEModelParameters> = GroupProjective<E>;

// generic over curve
pub struct PublicKey<E: TEModelParameters> {
    pub point: Point<E>,
}

impl<E: TEModelParameters> FromBytes for PublicKey<E> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> io::Result<Self> {
        let point = Point::read(&mut reader)?;
        Ok(PublicKey { 
            point 
        })
    }
}

impl<E: TEModelParameters> PublicKey<E> {
    // TODO: there are more params: data and so on...
    pub fn verify(&self, msg: &[u8], sig: &Signature, generator: ()) -> bool {
        // c = H*(Rbar || M)
        let c = h_star::<E>(&sig.rbar[..], msg);

        // Signature checks:
        // R != invalid
		let r = match Point::<E>::read(&sig.rbar[..]) {
			Ok(r) => r,
			Err(_) => return false,
		};

        // S < order(G)
        // (E::ScalarField guarantees its representation is in the field)
		let s = match E::ScalarField::read(&sig.sbar[..]) {
            Ok(s) => s,
            Err(_) => return false,
        };

        self.point.mul(&c).add(&r);
        //self.0.mul(c, params).add(&r, params).add(
            //&params.generator(p_g).mul(s, params).negate().into(),
            //params
        //).mul_by_cofactor(params).eq(&Point::zero())



        unimplemented!();
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
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}

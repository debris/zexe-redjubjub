use crate::{constants, group_hash::find_group_hash, Point};
use algebra::{biginteger::BigInteger256, PrimeField, TEModelParameters};

/// Fixed generators of the Jubjub curve of unknown
/// exponent.
#[derive(Copy, Clone)]
pub enum FixedGenerators {
    /// The prover will demonstrate knowledge of discrete log
    /// with respect to this base when they are constructing
    /// a proof, in order to authorize proof construction.
    ProofGenerationKey = 0,

    /// The note commitment is randomized over this generator.
    NoteCommitmentRandomness = 1,

    /// The node commitment is randomized again by the position
    /// in order to supply the nullifier computation with a
    /// unique input w.r.t. the note being spent, to prevent
    /// Faerie gold attacks.
    NullifierPosition = 2,

    /// The value commitment is used to check balance between
    /// inputs and outputs. The value is placed over this
    /// generator.
    ValueCommitmentValue = 3,
    /// The value commitment is randomized over this generator,
    /// for privacy.
    ValueCommitmentRandomness = 4,

    /// The spender proves discrete log with respect to this
    /// base at spend time.
    SpendingKeyGenerator = 5,
}

impl FixedGenerators {
    // TODO: cache value
    pub fn point<E>(&self) -> Point<E>
    where
        E: TEModelParameters,
        E::BaseField: PrimeField + Into<BigInteger256>,
    {
        match *self {
            FixedGenerators::ProofGenerationKey => find_group_hash(
                b"",
                constants::PROOF_GENERATION_KEY_BASE_GENERATOR_PERSONALIZATION,
            ),
            FixedGenerators::NoteCommitmentRandomness => {
                find_group_hash(b"r", constants::PEDERSEN_HASH_GENERATORS_PERSONALIZATION)
            }
            FixedGenerators::NullifierPosition => find_group_hash(
                b"",
                constants::NULLIFIER_POSITION_IN_TREE_GENERATOR_PERSONALIZATION,
            ),
            FixedGenerators::ValueCommitmentValue => {
                find_group_hash(b"v", constants::VALUE_COMMITMENT_GENERATOR_PERSONALIZATION)
            }
            FixedGenerators::ValueCommitmentRandomness => {
                find_group_hash(b"r", constants::VALUE_COMMITMENT_GENERATOR_PERSONALIZATION)
            }
            FixedGenerators::SpendingKeyGenerator => {
                find_group_hash(b"", constants::SPENDING_KEY_GENERATOR_PERSONALIZATION)
            }
        }
    }
}

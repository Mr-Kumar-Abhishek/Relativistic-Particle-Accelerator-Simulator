//! Electromagnetic field container and linear superposition helpers.

use core::ops::Add;

use crate::vec3::Vec3;

/// Combined 3D electric field `E` (V/m) and magnetic flux density `B` (T) at a space-time point.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct EmField {
    pub e: Vec3,
    pub b: Vec3,
}

impl EmField {
    /// Zero electromagnetic field.
    #[inline]
    pub const fn zero() -> Self {
        Self {
            e: Vec3::zero(),
            b: Vec3::zero(),
        }
    }

    /// Constructs a combined electromagnetic field.
    #[inline]
    pub const fn new(e: Vec3, b: Vec3) -> Self {
        Self { e, b }
    }

    /// Constructs a pure magnetic field (`E = 0`).
    #[inline]
    pub const fn magnetic(b: Vec3) -> Self {
        Self {
            e: Vec3::zero(),
            b,
        }
    }

    /// Constructs a pure electric field (`B = 0`).
    #[inline]
    pub const fn electric(e: Vec3) -> Self {
        Self {
            e,
            b: Vec3::zero(),
        }
    }
}

impl Add for EmField {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            e: self.e + rhs.e,
            b: self.b + rhs.b,
        }
    }
}

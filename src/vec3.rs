//! 3D Euclidean vector mathematics in double precision (`f64`).

use core::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

/// A 3-dimensional vector in IEEE-754 double precision (`f64`).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    /// Creates a new `Vec3` with the specified `(x, y, z)` components.
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Returns the zero vector `(0, 0, 0)`.
    #[inline]
    pub const fn zero() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Computes the Euclidean dot product `self · other`.
    #[inline]
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Computes the 3D cross product `self × other`.
    #[inline]
    pub fn cross(self, other: Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// Computes the squared Euclidean norm `|self|^2`.
    #[inline]
    pub fn norm_sq(self) -> f64 {
        self.dot(self)
    }

    /// Computes the Euclidean norm `|self|`.
    #[inline]
    pub fn norm(self) -> f64 {
        self.norm_sq().sqrt()
    }

    /// Returns a unit-length vector along `self`, or `Vec3::zero()` if `|self| == 0`.
    #[inline]
    pub fn normalized(self) -> Self {
        let n = self.norm();
        if n > 0.0 {
            self / n
        } else {
            Self::zero()
        }
    }

    /// Scales the vector by a scalar factor.
    #[inline]
    pub fn scale(self, s: f64) -> Self {
        self * s
    }

    /// Computes the Euclidean distance between two position vectors.
    #[inline]
    pub fn dist(self, other: Self) -> f64 {
        (self - other).norm()
    }

    /// Rotates the vector around the Y axis by `angle_rad` radians (right-handed).
    #[inline]
    pub fn rotate_y(self, angle_rad: f64) -> Self {
        let (sin_a, cos_a) = angle_rad.sin_cos();
        Self {
            x: self.x * cos_a + self.z * sin_a,
            y: self.y,
            z: -self.x * sin_a + self.z * cos_a,
        }
    }

    /// Rotates the vector around the Z axis by `angle_rad` radians (right-handed).
    #[inline]
    pub fn rotate_z(self, angle_rad: f64) -> Self {
        let (sin_a, cos_a) = angle_rad.sin_cos();
        Self {
            x: self.x * cos_a - self.y * sin_a,
            y: self.x * sin_a + self.y * cos_a,
            z: self.z,
        }
    }
}

impl Add for Vec3 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Sub for Vec3 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: f64) -> Self {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;
    #[inline]
    fn mul(self, rhs: Vec3) -> Vec3 {
        Vec3::new(self * rhs.x, self * rhs.y, self * rhs.z)
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;
    #[inline]
    fn div(self, rhs: f64) -> Self {
        let inv = 1.0 / rhs;
        Self::new(self.x * inv, self.y * inv, self.z * inv)
    }
}

impl Neg for Vec3 {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl SubAssign for Vec3 {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::FRAC_PI_2;

    #[test]
    fn cross_and_dot_products_are_orthogonal() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(-4.0, 5.0, -6.0);
        let c = a.cross(b);
        assert!(c.dot(a).abs() < 1.0e-14);
        assert!(c.dot(b).abs() < 1.0e-14);
    }

    #[test]
    fn rotations_preserve_norm() {
        let v = Vec3::new(1.0, 0.0, 0.0);
        let ry = v.rotate_y(FRAC_PI_2);
        assert!((ry.x).abs() < 1.0e-15);
        assert!((ry.z + 1.0).abs() < 1.0e-15);
        let rz = v.rotate_z(FRAC_PI_2);
        assert!((rz.x).abs() < 1.0e-15);
        assert!((rz.y - 1.0).abs() < 1.0e-15);
    }
}

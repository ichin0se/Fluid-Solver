use std::ops::{Add, AddAssign, Div, Mul, Sub, SubAssign};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vector2 {
    pub x: f32,
    pub y: f32,
}

impl Vector2 {
    pub fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }

    pub fn zero() -> Self {
        Self::new(0.0, 0.0)
    }

    pub fn dot(self, other: Self) -> f32 {
        self.x * other.x + self.y * other.y
    }

    pub fn length(self) -> f32 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    pub fn normalized_or_zero(self) -> Self {
        let length = self.length();
        if length > 1e-6 {
            self / length
        } else {
            Self::zero()
        }
    }
}

impl Add for Vector2 {
    type Output = Self;

    fn add(self, right_hand_side: Self) -> Self::Output {
        Self::new(self.x + right_hand_side.x, self.y + right_hand_side.y)
    }
}

impl Sub for Vector2 {
    type Output = Self;

    fn sub(self, right_hand_side: Self) -> Self::Output {
        Self::new(self.x - right_hand_side.x, self.y - right_hand_side.y)
    }
}

impl Mul<f32> for Vector2 {
    type Output = Self;

    fn mul(self, scalar_value: f32) -> Self::Output {
        Self::new(self.x * scalar_value, self.y * scalar_value)
    }
}

impl Div<f32> for Vector2 {
    type Output = Self;

    fn div(self, scalar_value: f32) -> Self::Output {
        Self::new(self.x / scalar_value, self.y / scalar_value)
    }
}

impl AddAssign for Vector2 {
    fn add_assign(&mut self, right_hand_side: Self) {
        self.x += right_hand_side.x;
        self.y += right_hand_side.y;
    }
}

impl SubAssign for Vector2 {
    fn sub_assign(&mut self, right_hand_side: Self) {
        self.x -= right_hand_side.x;
        self.y -= right_hand_side.y;
    }
}

#[cfg(test)]
mod tests {
    use super::Vector2;

    #[test]
    fn dot_product_matches_expected_value() {
        let first_vector = Vector2::new(3.0, 4.0);
        let second_vector = Vector2::new(-2.0, 5.0);

        let dot_product = first_vector.dot(second_vector);

        assert!((dot_product - 14.0).abs() < 1e-6);
    }

    #[test]
    fn normalized_vector_has_unit_length() {
        let vector = Vector2::new(3.0, 4.0);

        let normalized_vector = vector.normalized_or_zero();

        assert!((normalized_vector.length() - 1.0).abs() < 1e-6);
    }
}

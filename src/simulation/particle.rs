use crate::simulation::vector2::Vector2;

#[derive(Clone, Debug)]
pub struct Particle {
    pub position: Vector2,
    pub velocity: Vector2,
    pub affine_velocity_matrix: [[f32; 2]; 2],
    pub mass_kilograms: f32,
}

impl Particle {
    pub fn new(position: Vector2) -> Self {
        Self {
            position,
            velocity: Vector2::zero(),
            affine_velocity_matrix: [[0.0; 2]; 2],
            mass_kilograms: 1.0,
        }
    }
}

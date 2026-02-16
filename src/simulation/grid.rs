use crate::simulation::vector2::Vector2;

#[derive(Clone, Debug)]
pub struct Grid {
    pub node_count_x: usize,
    pub node_count_y: usize,
    pub cell_size_meters: f32,
    pub node_velocity: Vec<Vector2>,
    pub node_mass: Vec<f32>,
    pub node_pressure: Vec<f32>,
    pub node_divergence: Vec<f32>,
}

impl Grid {
    pub fn new(node_count_x: usize, node_count_y: usize, cell_size_meters: f32) -> Self {
        let node_count_total = node_count_x * node_count_y;
        Self {
            node_count_x,
            node_count_y,
            cell_size_meters,
            node_velocity: vec![Vector2::zero(); node_count_total],
            node_mass: vec![0.0; node_count_total],
            node_pressure: vec![0.0; node_count_total],
            node_divergence: vec![0.0; node_count_total],
        }
    }

    pub fn node_index(&self, node_x: usize, node_y: usize) -> usize {
        node_y * self.node_count_x + node_x
    }

    pub fn reset_accumulators(&mut self) {
        self.node_velocity.fill(Vector2::zero());
        self.node_mass.fill(0.0);
        self.node_divergence.fill(0.0);
    }

    pub fn world_position_of_node(&self, node_x: usize, node_y: usize) -> Vector2 {
        Vector2::new(
            node_x as f32 * self.cell_size_meters,
            node_y as f32 * self.cell_size_meters,
        )
    }

    pub fn dimensions_meters(&self) -> Vector2 {
        Vector2::new(
            (self.node_count_x - 1) as f32 * self.cell_size_meters,
            (self.node_count_y - 1) as f32 * self.cell_size_meters,
        )
    }
}

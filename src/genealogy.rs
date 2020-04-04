//! Genealogy of a population.
//!
//! Describes ancestors of a group of individuals until they meet in the first 
//! common ancestor of the group in a distant past. Usually coming from a realization of 
//! the coalescent process. 
//! 

// Structs
use partitions::PartitionVec;
use petgraph::Graph;
use std::collections::HashMap;

// Traits
use std::iter::FromIterator;

/// Genealogic tree. 
/// 
/// This struct is created by the ``sample_genealogy`` method on Coalescent<R>. 
/// See its documentation for more.
#[derive(Debug, Clone)]
pub struct Genealogy {
	path: Vec<PartitionVec<()>>, // including initial state
	steps: Vec<[usize; 2]>,
	time_steps: Vec<f64>, // all positive intervals
	graph: Option<Graph<(usize, usize), f64, petgraph::Undirected, u32>>,
}

impl Genealogy {

	pub(crate) fn new(path: Vec<PartitionVec<()>>, steps: Vec<[usize; 2]>, time_steps: Vec<f64>) -> Self {
		let graph = None;

		Genealogy{path, steps, time_steps, graph}
	}

	/// Total depth of the tree, i.e. the distance from the first common ancestor
	/// of the group. 
	pub fn depth(&self) -> f64 {
		self.time_steps.iter().sum()
	}

	/// Sum of all the time represented in the edges of the genealogy. 
	pub fn length(&self) -> f64 {
		let group_size = self.path[0].len();
		self.time_steps
			.iter()
			.enumerate()
			.map(|(i, time_step)| (group_size - i) as f64 * time_step )
			.sum()
	}

	/// Distance between two individuals in the genealogic tree. 
	/// 
	/// # Remarks
	/// 
	/// To compute mean divergence over all pairs, prefer the method
	/// ``mean_pairwise_divergence``. 
	pub fn divergence(&self, index_1: usize, index_2: usize) -> f64 {
		let mut counter = 0;
		for state in &self.path {
			if state.same_set(index_1, index_2) {
				break;
			} else {
				counter += 1;
			}
		}

		2.0 * self.time_steps.iter().take(counter).sum::<f64>()
	}

	/// Mean distance of all pairs of individual through their first common ancestor, i.e. 
	/// mean distance of all pairs of leaves in the tree. 
	pub fn mean_pairwise_divergence(&self) -> f64 {
		let group_size = self.steps.len() + 1;
		let mut cummulative_time = 0.0;
		let mut cummulative_divergence = 0.0;

		for iteration in 0..(group_size - 1) {
			// Retrieve values

			let value_indexes = self.steps[iteration];
			let state = &self.path[iteration];
			cummulative_time += self.time_steps[iteration];

			// Identify the lengths of sets joints

			let set_size_1 = state.len_of_set(value_indexes[0]);
			let set_size_2 = state.len_of_set(value_indexes[1]); 
			
			// Count number of pairs

			let number_of_pairs = set_size_1 * set_size_2;

			// Add to the counter the respective time

			cummulative_divergence += (2.0 * cummulative_time) * number_of_pairs as f64;
		}

		// Take mean

		cummulative_divergence * 2.0 / (group_size * (group_size - 1)) as f64
	}

	fn compute_graph(&mut self) -> &Graph<(usize, usize), f64, petgraph::Undirected, u32> { 
		let group_size = self.steps.len() + 1;
		let mut graph = Graph::new_undirected();
		let mut node_indexes = HashMap::new();
		let mut representatives_generation = HashMap::new();

		match self.steps.is_empty() {
			true => {
				graph.add_node((0, 0));
			},
			false => {
				let mut state: PartitionVec<()> =
            		PartitionVec::from_iter((0..group_size).map(|_| ()));
				for index in 0..group_size {
					let node_index = graph.add_node((0, index));
					node_indexes.insert((0, index), node_index);
					representatives_generation.insert(index, 0);
				}

				for generation in 0..self.steps.len() {
					let time_step = self.time_steps[generation];
					let value_indexes = self.steps[generation];

					// Retrieve representatives

					let representatives: Vec<usize> = value_indexes.iter()
						.map(|&i| (0..group_size).find(|&j| state.same_set(i, j))
							.expect("Could not retrieve previous state of the genealogy.")
						).collect();

					// Add parent node

					let new_representative = core::cmp::min(representatives[0], representatives[1]);
					let node_index = graph.add_node((generation + 1, new_representative));
					node_indexes.insert((generation + 1, new_representative), node_index);
					
					// Add edges

					graph.add_edge(
						node_indexes[&(generation + 1, new_representative)], 
						node_indexes[&(representatives_generation[&representatives[0]], representatives[0])], 
						time_step
					);
					graph.add_edge(
						node_indexes[&(generation + 1, new_representative)], 
						node_indexes[&(representatives_generation[&representatives[1]], representatives[1])], 
						time_step
					);

					// Update

					representatives_generation.insert(new_representative, generation + 1);
					state.union(value_indexes[0], value_indexes[1]);
				}
			},
		}
		
		self.graph = Some(graph);
		&self.graph.as_ref().unwrap()
	}
}

impl Into<Graph<(usize, usize), f64, petgraph::Undirected, u32>> for Genealogy 
{
	fn into(mut self) -> Graph<(usize, usize), f64, petgraph::Undirected, u32> { 
		match self.graph {
		 	Some(graph) => graph,
		 	None => self.compute_graph().clone(),
		 } 
	}
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn construction() {
		let group_size = 100;
		let rng = rand::thread_rng();
		let coalescent = crate::Coalescent::new(group_size, rng);
		
		let mut rng = rand::thread_rng();
		let genealogy: Genealogy = coalescent.sample_genealogy(&mut rng);

		assert_eq!(genealogy.path.len(), group_size);
		assert_eq!(genealogy.steps.len(), group_size - 1);
		assert_eq!(genealogy.time_steps.len(), group_size - 1);
		assert!(genealogy.graph.is_none());
	}
}
//! Coalescent process.
//!
//! Describes genealogies and is the limit process of various finite population
//! models for ancestry. The n-coalescent process is a Markovian process
//! over the set of partitions of the set {1, 2, ..., n}. It starts with the
//! partition of only singletons and it can go to any partition that
//! uses one set less than the current partition. All options have the same
//! rate one. In particular, the first step is taken after a period that
//! distributes as an exponential with mean 2 / (n (n - 1)), the inverse
//! of the number of possible partitions.
//! 

// Types
use partitions::PartitionVec;
use rand_distr::Exp;
use crate::Genealogy;

// Traits
use markovian::traits::CMarkovChainTrait;
use rand::distributions::Distribution;
use rand::seq::IteratorRandom;
use rand::Rng;
use core::fmt::Debug;
use std::iter::FromIterator;


/// n-Coalescent process in the space of partitions of the set {1, 2, ..., n}.
/// Starts with a finite partition of all singletons and it ends with a single set.
///
/// It has a random number generator associated, R, to be a random iterator.  
///
/// A Coalescent can be seen as:
/// - State-iterator: an iterator with a current state, changing randomly to another
/// state when ``next`` method is called. See
/// [Iterator](https://doc.rust-lang.org/nightly/core/iter/trait.Iterator.html)
/// and [MarkovChainTrait](file:///C:/Users/rasau/projects/markovian/target/doc/markovian/discrete_time/struct.MarkovChain.html)
/// implementation.
/// - Random genealogy generator: random variable over possible genealogies from the
/// current state. See method [sample_genealogy](file:///C:/Users/rasau/projects/coalescence/target/doc/coalescence/coalescent/struct.Coalescent.html#method.sample_genealogy).
#[derive(Debug, Clone)]
pub struct Coalescent<R>
where
    R: Rng + Clone + core::fmt::Debug,
{
    state: PartitionVec<()>, // No selection
    rng: R,
}

impl<R> Coalescent<R>
where
    R: Rng + Clone + core::fmt::Debug,
{
    /// Creates a new Coalescent. 
    /// 
    /// # Examples
    /// 
    /// Construction of a random walk in the integers, using a closure.
    /// ```
    /// # #![allow(unused_mut)]
    /// let group_size = 100;
    /// let rng = rand::thread_rng();
    /// let mut coalescent = coalescence::Coalescent::new(group_size, rng);
    /// ``` 
    pub fn new(group_size: usize, rng: R) -> Self {
        let state: PartitionVec<()> =
            PartitionVec::from_iter((0..group_size).map(|_| ()));

        Coalescent { state, rng }
    }

    /// Mutable reference to the internal random number generator. 
    /// 
    /// # Remarks
    /// 
    /// This method is for debugging mainly. 
    /// 
    pub(crate) fn rng(&mut self) -> &mut R {
        &mut self.rng
    }

    /// Change the internal random number generator for another. 
    /// 
    /// # Remarks
    /// 
    /// This method is for debugging mainly. If the random number generator 
    /// is deterministic, we can restart the random process. 
    /// 
    /// # Examples
    /// 
    /// Using non-deterministic random number generators we do not get 
    /// reproducible results. 
    /// ```
    /// let group_size = 100;
    /// let rng = rand::thread_rng();
    /// let mut coalescent = coalescence::Coalescent::new(group_size, rng.clone());
    /// 
    /// let result_1 = coalescent.next();
    /// coalescent = coalescence::Coalescent::new(group_size, rng.clone());
    /// let result_2 = coalescent.next();
    /// 
    /// assert!(result_1 != result_2);
    /// ``` 
    /// 
    /// Using deterministic random number generators we do get 
    /// reproducible results. 
    /// ```
    /// use rand::SeedableRng;
    /// 
    /// let group_size = 100;
    /// let mut rng = rand_pcg::Pcg32::seed_from_u64(123);
    /// let mut coalescent = coalescence::Coalescent::new(group_size, rng.clone());
    /// 
    /// let result_1 = coalescent.next();
    /// coalescent = coalescence::Coalescent::new(group_size, rng.clone());
    /// let result_2 = coalescent.next();
    /// 
    /// assert_eq!(result_1, result_2);
    /// ```
    pub fn set_rng(&mut self, other_rng: R) -> &mut Self {
        self.rng = other_rng;
        self
    }

    /// Sample a path: from the current state until there is only one set
    /// in the partition. Returns a path represented by a vector
    /// where the first element is ``(0.0, initial_state)`` and 
    /// the other elements are the (n - 1) realizations
    /// the takes to get to a partition with only one set. 
    ///
    /// # Remarks
    ///
    /// No internal state changes, including the internal
    /// random number generator. This is why this methods requires a rng.  
    ///
    /// # Examples
    ///
    /// 
    /// ```
    /// let group_size = 100;
    /// let rng = rand::thread_rng();
    /// let mut coalescent = coalescence::Coalescent::new(group_size, rng);
    ///
    /// let mut rng = rand::thread_rng();
    /// let path = coalescent.sample_path(&mut rng);
    ///
    /// assert_eq!(path.len(), group_size);
    /// ```
    ///
    pub fn sample_path<S>(&self, rng: &mut S) -> Vec<(f64, PartitionVec<()>)> 
    where
        S: Rng + Clone + Debug,
    {
        // Initialize a Coalescent

        let initial_group_size: usize = self.state().len();
        let other_rng = rng.clone();
        let mut coalescent_process = Coalescent::new(initial_group_size, other_rng);

        // Generate a realizations

        let mut realizations = vec![(0.0, coalescent_process.state().clone())];
        loop {
            let (time_step, state) = match coalescent_process.next() {
                Some((t, s)) => (t, s),
                None => break,
            };
            realizations.push((time_step, state));
        }

        // Update rng

        *rng = coalescent_process.rng().clone();

        // Finish

        realizations
    }

    /// Peeks a possible next pair of indices to be joint, chosen 
    /// according to the stochastic process. This does not change the 
    /// state of the ``Coalescent``. 
    /// 
    /// # Examples
    /// 
    /// ```
    /// use markovian::traits::CMarkovChainTrait;
    /// let group_size = 2;
    /// let rng = rand::thread_rng();
    /// let mut coalescent = coalescence::Coalescent::new(group_size, rng);
    ///
    /// // The next state will be a partition with (group_size - 1) sets. 
    /// let (_time, index_pair) = coalescent.peek_next_step().expect("The simulation process failed!");
    /// let current_state = coalescent.state();
    /// 
    /// assert_eq!(2, current_state.amount_of_sets());
    /// assert!(!current_state.same_set(index_pair[0], index_pair[1]));
    /// ``` 
    pub fn peek_next_step(&mut self) -> Option<(f64, [usize; 2])> {
        let current_partition_size = self.state.amount_of_sets();

        if current_partition_size == 1 {
            None
        } else {
            // Simulate time step

            let rate = (current_partition_size * (current_partition_size - 1) / 2) as f64;
            let exp = Exp::new(rate).unwrap();
            let time_step = exp.sample(&mut self.rng());

            // Choose between possible transitions

            let mut set_indexes = [0; 2];
            (0..current_partition_size).choose_multiple_fill(&mut self.rng(), &mut set_indexes);

            // Get values from these sets
            let value_indexes: Vec<usize> = (0..2)
                .map(|i| {
                    self.state
                        .all_sets()
                        .nth(set_indexes[i])
                        .unwrap() // set
                        .nth(0)
                        .unwrap() // (value_index, value)
                        .0
                })
                .collect();
            let value_indexes = [value_indexes[0], value_indexes[1]];

            // Return

            Some((time_step, value_indexes))
        }
    }

    /// Changes to a next state of the ``Coalescent``, chosen 
    /// according to the stochastic process and returning the indexes of 
    /// elements that represent the sets of the partitions that were joint
    /// to produce this next state. 
    /// 
    /// # Examples
    /// 
    /// ```
    /// let group_size = 2;
    /// let rng = rand::thread_rng();
    /// let mut coalescent = coalescence::Coalescent::new(group_size, rng);
    ///
    /// // The next state will be a partition with (group_size - 1) sets. 
    /// let (_time, value_indexes) = coalescent.peek_next_step().expect("The simulation process failed!");
    /// assert!(value_indexes[0] != value_indexes[1]); 
    /// assert!(value_indexes[0] < group_size && value_indexes[1] < group_size ); 
    /// ``` 
    pub fn next_step(&mut self) -> Option<(f64, [usize; 2])> {
        match self.peek_next_step() {
            Some((time_step, value_indexes)) => {
                self.state.union(value_indexes[0], value_indexes[1]);
                Some((time_step, value_indexes))
            },
            None => None,
        }
    }

    /// Sample a genealogy: from the current state until there is only one set
    /// in the partition. Returns a ``Genealogy`` where postprocess is possible. 
    ///
    /// # Remarks
    ///
    /// No internal state changes, including the internal
    /// random number generator. This is why this methods requires a rng.  
    ///
    /// # Examples
    ///
    /// 
    /// ```
    /// let group_size = 100;
    /// let rng = rand::thread_rng();
    /// let mut coalescent = coalescence::Coalescent::new(group_size, rng);
    ///
    /// let mut rng = rand::thread_rng();
    /// let genealogy = coalescent.sample_genealogy(&mut rng);
    /// ```
    ///
    pub fn sample_genealogy<S>(&self, rng: &mut S) -> Genealogy
    where
        S: Rng + Clone + Debug,
    {
       // Initialize a Coalescent

        let group_size: usize = self.state().len();
        let other_rng = rng.clone();
        let mut coalescent_process = Coalescent::new(group_size, other_rng);

        // Generate a transitions

        let mut path: Vec<PartitionVec<()>> = Vec::with_capacity(group_size - 1);
        let mut steps: Vec<[usize; 2]> = Vec::with_capacity(group_size - 1);
        let mut time_steps: Vec<f64> = Vec::with_capacity(group_size - 1);

        let mut state: PartitionVec<()> =
            PartitionVec::from_iter((0..group_size).map(|_| ()));

        path.push(state.clone());
        loop {
            let (time_step, value_indexes) = match coalescent_process.next_step() {
                Some(step) => step,
                None => break,
            };
            state.union(value_indexes[0], value_indexes[1]);

            path.push(state.clone());
            steps.push(value_indexes);
            time_steps.push(time_step);
        }

        // Update rng

        *rng = coalescent_process.rng().clone();

        // Finish

        Genealogy::new(path, steps, time_steps)
    }
}

impl<R> CMarkovChainTrait<PartitionVec<()>> for Coalescent<R>
where
    R: Rng + Clone + core::fmt::Debug,
{
    /// Current state of the process. 
    fn state(&self) -> &PartitionVec<()> {
        &self.state
    }

    /// Change the current state of the process.
    fn set_state(&mut self, state: PartitionVec<()>) -> &mut Self {
        self.state = state;
        self
    }
}

impl<R> Iterator for Coalescent<R>
where
    R: Rng + Clone + core::fmt::Debug,
{
    type Item = (f64, PartitionVec<()>);


    /// Changes the state of the ``Coalescent`` to a new state, chosen 
    /// according to the stochastic process. 
    /// 
    /// # Examples
    /// 
    /// ```
    /// let group_size = 2;
    /// let rng = rand::thread_rng();
    /// let mut coalescent = coalescence::Coalescent::new(group_size, rng);
    ///
    /// // The next state is a partition with (group_size - 1) sets. 
    /// let (_time, new_partition) = coalescent.next().expect("The simulation process failed!");
    /// assert_eq!(group_size - 1, new_partition.amount_of_sets()); 
    /// ``` 
    fn next(&mut self) -> Option<Self::Item> {
        match self.peek_next_step() {
            Some((time_step, value_indexes)) => {
                self.state.union(value_indexes[0], value_indexes[1]);
                Some((time_step, self.state.clone()))
            },
            None => None,
        }
    }
}
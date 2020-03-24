use itertools::Itertools;
use itertools_num::ItertoolsNum;
use partitions::PartitionVec;
use preexplorer::prelude::*;

use markovian::CMarkovChain;
// use numerical_algorithms::markov_chain::traits::CMarkovChainTrait;

fn main() {
    let mut init_state: PartitionVec<usize> = PartitionVec::new();
    for i in 0..20 {
        init_state.push(i);
    }
    fn transition(state: PartitionVec<usize>) -> Vec<(PartitionVec<usize>, f64)> {
        let mut possibilities = Vec::new();
        for (index1, index2) in (0..state.len()).tuple_combinations() {
            let mut possibility = state.clone();
            possibility.union(index1, index2);
            possibilities.push((possibility, 1.0));
        }
        possibilities
    }

    let mc = CMarkovChain::new(init_state.clone(), transition);

    // Realization until extinction
    let mut realization = vec![(0.0, init_state)];
    let mut unfinished = true;
    realization.append(
        &mut mc
            .take_while(|(_, x)| {
                let this_time = unfinished;
                if x.amount_of_sets() == 1 {
                    unfinished = false;
                }
                this_time
            })
            .collect::<Vec<(_, _)>>(),
    );

    if true {
        plot_group_size(realization.clone());
    }

    plot_genealogy(realization);
}

fn plot_genealogy(realization: Vec<(f64, PartitionVec<usize>)>) {
    let (delta_times, values): (Vec<_>, Vec<_>) = realization.into_iter().unzip();
    let times: Vec<f64> = delta_times.iter().cumsum().collect();

    let num_individuals = values[0].amount_of_sets();
    let mut all_trajectories = Vec::with_capacity(num_individuals);

    for i in 0..num_individuals {
        let local_trajectory: Vec<usize> = values
            .iter()
            // Find the set which individual "i" belongs to
            .map(|partition| {
                (0..num_individuals)
                    .find(|&j| partition.same_set(i, j))
                    .unwrap()
            })
            .collect();
        all_trajectories.push((&times, local_trajectory).preexplore());
    }

    pre::process::Comparison::new(all_trajectories)
        .dashtype(0)
        .style("line")
        .plot("genealogy")
        .unwrap();
}

fn plot_group_size(realization: Vec<(f64, PartitionVec<usize>)>) {
    let (delta_times, values): (Vec<_>, Vec<_>) = realization.into_iter().unzip();

    let times: Vec<f64> = delta_times.iter().cumsum().collect();
    let gropu_size = values.iter().map(|x| x.amount_of_sets());

    (times, gropu_size)
        .preexplore()
        .title("Group size of coallescent process")
        .labelx("time")
        .labely("size")
        .plot("test")
        .unwrap();
}

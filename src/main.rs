// Crates
use preexplorer::prelude::*;
use rand::prelude::*;

// Structs
use coalescence::Coalescent;
use partitions::PartitionVec;

// Traits
use itertools_num::ItertoolsNum;
use markovian::traits::CMarkovChainTrait;

fn main() {
    let rng = thread_rng();
    let mut coalescent = Coalescent::new(500, rng);

    // Print realizations

    if false {
        println!("{:?}", coalescent.state());
        println!("{:?}", coalescent.next());
        println!("{:?}", coalescent.next());
        println!("{:?}", coalescent.next());
        println!("{:?}", coalescent.next());
        println!("{:?}", coalescent.next());
    }

    // Plot group size

    if true {
        let realization = coalescent.generate_realization();
        println!("{:?}", coalescent.state());
        plot_group_size(realization);
    }

    // Plot genealogy

    if true {
        let realization = coalescent.generate_realization();
        plot_genealogy(realization);
    }
}

fn plot_genealogy(realization: Vec<(f64, PartitionVec<()>)>) {
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

fn plot_group_size(realization: Vec<(f64, PartitionVec<()>)>) {
    let (delta_times, values): (Vec<_>, Vec<_>) = realization.into_iter().unzip();

    let times: Vec<f64> = delta_times.iter().cumsum().collect();
    let gropu_size = values.iter().map(|x| x.amount_of_sets());

    (times, gropu_size)
        .preexplore()
        .title("Group size of coalescent process")
        .labelx("time")
        .labely("size")
        .plot("test")
        .unwrap();
}

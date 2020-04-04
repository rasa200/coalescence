// Crates
use preexplorer::prelude::*;
use rayon::prelude::*;

// Functions

use rand::thread_rng;

// Structs
use coalescence::Coalescent;


// Traits



fn main() {

    // Compare depth: empirical vs theoretical
    if false {
        let power = 10;
        let group_sizes: Vec<usize> = (1..=power).map(|i| 2_usize.pow(i as u32)).collect();
        let samples = 1000;

        let (empirical, theoretical) = depth_comparison(&group_sizes, samples);
        plot_comparison(&group_sizes, empirical, theoretical, "Depth");
    }

    // Compare length: empirical vs theoretical
    if false {
        let power = 10;
        let group_sizes: Vec<usize> = (1..=power).map(|i| 2_usize.pow(i as u32)).collect();
        let samples = 10;

        let (empirical, theoretical) = length_comparison(&group_sizes, samples);
        plot_comparison(&group_sizes, empirical, theoretical, "Length");
    }

    // Compare pairwise divergence: empirical vs theoretical
    if true {
        let power = 10;
        let group_sizes: Vec<usize> = (1..=power).map(|i| 2_usize.pow(i as u32)).collect();
        let samples = 1000;

        let (empirical, theoretical) = pairwise_divergence_comparison(&group_sizes, samples);
        plot_comparison(&group_sizes, empirical, theoretical, "Pairwise divergence");
    }
}

fn pairwise_divergence_comparison(group_sizes: &Vec<usize>, samples: usize) -> (Vec<f64>, Vec<f64>) {
    let mut empirical_depths = Vec::new();
    let mut theoretical_depths = Vec::new();

    for group_size in group_sizes {
        let empirical_depth = (0..samples)
            .collect::<Vec<usize>>()
            .par_iter()
            .map(|_| {
                let coalescent = Coalescent::new(*group_size, thread_rng());
                let mut rng = thread_rng();
                coalescent.sample_genealogy(&mut rng).mean_pairwise_divergence()
            })
            .sum::<f64>() / samples as f64;

        fn theoretical_iter(x: f64, counter: usize) -> f64 {
            let n: f64 = counter as f64;
            let p: f64 = 2.0 / (n * (n - 1.0) );
            p * (2.0 * p) + (1.0 - p) * (2.0 * p + x)
        }
        let mut counter = 2;
        let mut x = 2.0;
        while counter < *group_size {
            x = theoretical_iter(x, counter);
            counter += 1;
        }
        let theoretical_depth = x;
        
        empirical_depths.push(empirical_depth);
        theoretical_depths.push(theoretical_depth);
    }
    
    (empirical_depths, theoretical_depths)
}

fn length_comparison(group_sizes: &Vec<usize>, samples: usize) -> (Vec<f64>, Vec<f64>) {
    let mut empirical_depths = Vec::new();
    let mut theoretical_depths = Vec::new();

    for group_size in group_sizes {
        let empirical_depth = (0..samples)
            .collect::<Vec<usize>>()
            .par_iter()
            .map(|_| {
                let coalescent = Coalescent::new(*group_size, thread_rng());
                let mut rng = thread_rng();
                coalescent.sample_genealogy(&mut rng).length()
            })
            .sum::<f64>() / samples as f64;

        let theoretical_depth = 2.0 * ((*group_size as f64 - 1.0).ln() + 0.57721 + 1.0 / (2 * (*group_size  - 1) ) as f64 );
        
        empirical_depths.push(empirical_depth);
        theoretical_depths.push(theoretical_depth);
    }
    
    (empirical_depths, theoretical_depths)
}

fn plot_comparison(group_sizes: &Vec<usize>, empirical: Vec<f64>, theoretical: Vec<f64>, topic: &str) {
    pre::process::Comparison::new(vec![
        (group_sizes, empirical)
            .preexplore()
            .title("empirical")
            .to_owned(),
        (group_sizes, theoretical)
            .preexplore()
            .title("theoretical")
            .to_owned(),
        
        ])
        .title(format!("{}: empirical mean vs expectation", topic))
        .labelx("initial group size")
        .labely("time")
        .logx(2)
        .plot(&format!("{}", topic))
        .unwrap();
}

/// # Output
/// 
/// (empirical, theoretical). 
fn depth_comparison(group_sizes: &Vec<usize>, samples: usize) -> (Vec<f64>, Vec<f64>) {
    let mut empirical_depths = Vec::new();
    let mut theoretical_depths = Vec::new();

    for group_size in group_sizes {
        let empirical_depth = (0..samples)
            .collect::<Vec<usize>>()
            .par_iter()
            .map(|_| {
                let coalescent = Coalescent::new(*group_size, thread_rng());
                let mut rng = thread_rng();
                coalescent.sample_genealogy(&mut rng).depth()
            })
            .sum::<f64>() / samples as f64;

        let theoretical_depth = 2.0 * (1.0 - 1.0 / *group_size as f64);
        
        empirical_depths.push(empirical_depth);
        theoretical_depths.push(theoretical_depth);
    }
    
    (empirical_depths, theoretical_depths)
}
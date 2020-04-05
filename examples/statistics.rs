// Crates
use preexplorer::prelude::*;
use rayon::prelude::*;

// Functions
use rand::thread_rng;

// Structs
use coalescence::Coalescent;
use average::Variance;
use rand_distr::Poisson;


// Traits
use rand::distributions::Distribution;



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
    if false {
        let power = 10;
        let group_sizes: Vec<usize> = (1..=power).map(|i| 2_usize.pow(i as u32)).collect();
        let samples = 1000;

        let (empirical, theoretical) = pairwise_divergence_comparison(&group_sizes, samples);
        plot_comparison(&group_sizes, empirical, theoretical, "Pairwise divergence");
    }


    // Compare depth: empirical vs theoretical
    if false {
        let power = 10;
        let group_sizes: Vec<usize> = (1..=power).map(|i| 2_usize.pow(i as u32)).collect();
        let samples = 1000;

        let (empirical, _) = depth_comparison(&group_sizes, samples);
        plot_emp_var(&group_sizes, empirical, "Depth");
    }

    // Compare length: empirical vs theoretical
    if false {
        let power = 10;
        let group_sizes: Vec<usize> = (1..=power).map(|i| 2_usize.pow(i as u32)).collect();
        let samples = 10;

        let (empirical, _) = length_comparison(&group_sizes, samples);
        plot_emp_var(&group_sizes, empirical, "Length");
    }

    // Compare pairwise divergence: empirical vs theoretical
    if false {
        let power = 10;
        let group_sizes: Vec<usize> = (1..=power).map(|i| 2_usize.pow(i as u32)).collect();
        let samples = 1000;

        let (empirical, _) = pairwise_divergence_comparison(&group_sizes, samples);
        plot_emp_var(&group_sizes, empirical, "Pairwise divergence");
    }

    if false {
        let power = 10;
        let group_sizes: Vec<usize> = (1..=power).map(|i| 2_usize.pow(i as u32)).collect();
        let samples = 1000;

        let empirical = divergence_sample(&group_sizes, samples);
        plot_emp_var(&group_sizes, empirical, "Pair divergence");
    }


    if false {
        let power = 10;
        let group_sizes: Vec<usize> = (1..=power).map(|i| 2_usize.pow(i as u32)).collect();
        let samples = 1000;

        let empirical = heterozygosity_sample(&group_sizes, samples);
        plot_with_error(&group_sizes, empirical, "Heterozygosity");
    }
}

fn plot_with_error(group_sizes: &Vec<usize>, empirical: Vec<Variance>, topic: &str) {
    let mut data = Vec::new();
    for i in 0..empirical.len() {
        data.push(group_sizes[i] as f64);
        data.push(empirical[i].mean());
        data.push(empirical[i].error());
    }
    let dim = 3;

    pre::Data::new(data, dim)
        .id(topic)
        .title(topic)
        .labelx("initial group size")
        .labely("time")
        .logx(2)
        // .write_plot_script()
        // .unwrap()
        .save()
        .unwrap();
}


fn heterozygosity_sample(group_sizes: &Vec<usize>, samples: usize) -> Vec<Variance> {
    let mut empirical = Vec::new();

    for group_size in group_sizes {
        let empirical_sample: Variance = (0..samples)
            .collect::<Vec<usize>>()
            .par_iter()
            .map(|_| {
                let coalescent = Coalescent::new(*group_size, thread_rng());
                let mut rng = thread_rng();
                let depth = coalescent.sample_genealogy(&mut rng).divergence(0, 1) / 2.0;
                let poi = Poisson::new(depth).unwrap();

                let mutations_1: u64 = poi.sample(&mut rng);
                let mutations_2: u64 = poi.sample(&mut rng);
                ((mutations_1 - mutations_2) % 2 == 0) as u32 as f64
            })
            .collect::<Vec<f64>>()
            .iter()
            .collect();
        
        empirical.push(empirical_sample);
    }
    
    empirical
}


fn divergence_sample(group_sizes: &Vec<usize>, samples: usize) -> Vec<Variance> {
    let mut empirical_depths = Vec::new();

    for group_size in group_sizes {
        let empirical_depth: Variance = (0..samples)
            .collect::<Vec<usize>>()
            .par_iter()
            .map(|_| {
                let coalescent = Coalescent::new(*group_size, thread_rng());
                let mut rng = thread_rng();
                coalescent.sample_genealogy(&mut rng).divergence(0, 1)
            })
            .collect::<Vec<f64>>()
            .iter()
            .collect();
        
        empirical_depths.push(empirical_depth);
    }
    
    empirical_depths
}

fn plot_emp_var(group_sizes: &Vec<usize>, empirical: Vec<Variance>, topic: &str) {
    (group_sizes, empirical.iter().map(|v| v.sample_variance()))
        .preexplore()
        .title(format!("{}: empirical variance", topic))
        .labelx("initial group size")
        .labely("time")
        .logx(2)
        .plot(&format!("{}", topic))
        .unwrap();
}



fn pairwise_divergence_comparison(group_sizes: &Vec<usize>, samples: usize) -> (Vec<Variance>, Vec<f64>) {
    let mut empirical_depths = Vec::new();
    let mut theoretical_depths = Vec::new();

    for group_size in group_sizes {
        let empirical_depth: Variance = (0..samples)
            .collect::<Vec<usize>>()
            .par_iter()
            .map(|_| {
                let coalescent = Coalescent::new(*group_size, thread_rng());
                let mut rng = thread_rng();
                coalescent.sample_genealogy(&mut rng).mean_pairwise_divergence()
            })
            .collect::<Vec<f64>>()
            .iter()
            .collect();
        let theoretical_depth = 2.0;
        
        empirical_depths.push(empirical_depth);
        theoretical_depths.push(theoretical_depth);
    }
    
    (empirical_depths, theoretical_depths)
}

fn length_comparison(group_sizes: &Vec<usize>, samples: usize) -> (Vec<Variance>, Vec<f64>) {
    let mut empirical_depths = Vec::new();
    let mut theoretical_depths = Vec::new();

    for group_size in group_sizes {
        let empirical_depth: Variance = (0..samples)
            .collect::<Vec<usize>>()
            .par_iter()
            .map(|_| {
                let coalescent = Coalescent::new(*group_size, thread_rng());
                let mut rng = thread_rng();
                coalescent.sample_genealogy(&mut rng).length()
            })
            .collect::<Vec<f64>>()
            .iter()
            .collect();

        let theoretical_depth = 2.0 * ((*group_size as f64 - 1.0).ln() + 0.57721 + 1.0 / (2 * (*group_size  - 1) ) as f64 );
        
        empirical_depths.push(empirical_depth);
        theoretical_depths.push(theoretical_depth);
    }
    
    (empirical_depths, theoretical_depths)
}

fn plot_comparison(group_sizes: &Vec<usize>, empirical: Vec<Variance>, theoretical: Vec<f64>, topic: &str) {
    pre::process::Comparison::new(vec![
        (group_sizes, empirical.iter().map(|v| v.mean()).collect::<Vec<f64>>())
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
fn depth_comparison(group_sizes: &Vec<usize>, samples: usize) -> (Vec<Variance>, Vec<f64>) {
    let mut empirical_depths = Vec::new();
    let mut theoretical_depths = Vec::new();

    for group_size in group_sizes {
        let empirical_depth: Variance = (0..samples)
            .collect::<Vec<usize>>()
            .par_iter()
            .map(|_| {
                let coalescent = Coalescent::new(*group_size, thread_rng());
                let mut rng = thread_rng();
                coalescent.sample_genealogy(&mut rng).depth()
            })
            .collect::<Vec<f64>>()
            .iter()
            .collect();

        let theoretical_depth = 2.0 * (1.0 - 1.0 / *group_size as f64);
        
        empirical_depths.push(empirical_depth);
        theoretical_depths.push(theoretical_depth);
    }
    
    (empirical_depths, theoretical_depths)
}
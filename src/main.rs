fn main() {
    let group_size = 5;
    let rng = rand::thread_rng();
    let coalescent = coalescence::Coalescent::new(group_size, rng);

    let mut rng = rand::thread_rng();
    println!("{:?}", coalescent.sample_genealogy(&mut rng) );
}
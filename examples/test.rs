extern crate rand;

fn main() {
let dissim = ndarray::arr2(&[[0,1,2,3],[1,0,4,5],[2,4,0,6],[3,5,6,0]]);
let mut meds = kmedoids::random_initialization(4, 2, &mut rand::thread_rng());
let (loss, assingment, n_iter, n_swap) = kmedoids::fasterpam(&dissim, &mut meds, 100);
println!("Loss is: {}", loss);
}

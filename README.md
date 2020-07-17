# ConsProb, which Estimates Average Posterior Probabilities on Sparse Global Pairwise RNA Structural Alignment
# Installation
This project has been written in Rust, a systems programming language.
You need to install the Rust components, which are rustc (the Rust compiler), cargo (the Rust package manager), and the Rust standard library.
Visit [the Rust website](https://www.rust-lang.org) to see more about the language.
You can install the components with one line as follows:
```bash
$ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
The installation is arranged by [Rustup](https://github.com/rust-lang-nursery/rustup.rs), which enables to switch easily a compiler in use.
You can install ConsProb: 
```bash
$ cargo install consprob
```
Check if the program has been installed properly:
```bash
$ consprob # Its available command options will be displayed.
```
You can run the program with a prepared test set of sampled tRNAs:
```bash
$ git clone https://github.com/heartsh/consprob && cd consprob
$ cargo test --release -- --nocapture
```

# Author
[Heartsh](https://github.com/heartsh)

# License
Copyright (c) 2018 Heartsh  
Licensed under [the MIT license](http://opensource.org/licenses/MIT).

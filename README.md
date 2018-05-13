# STRAP Program, Iterative-Bayesian-estimation Program for RNA Structural-alignment Probability Matrix Quadruple
This project provides an iterative-Bayesian-estimation binary for RNA structural-alignment probability matrix quadruples.

# Dependencies
This binary depends on [the ParasoR program](https://github.com/carushi/ParasoR) for estimation of base-pairing probability matrices.
Install this program following the installation instruction over this link.

# Installation
This project has been written in mainly Rust, a systems programming language.
So first you need to install the Rust compiler, Rust package manager, and Rust standard library. 
Visit [the Rust website](https://www.rust-lang.org) to see more about this language.
You can install these 3 components with 1 line as follows:
```bash
$ curl https://sh.rustup.rs -sSf | sh
```
The above installation is done by [Rustup](https://github.com/rust-lang-nursery/rustup.rs), so you can easily switch a compiler to use. 
Now you can install the STRAP program as follows: 
```bash
$ cargo install strap
```
Check if this program has been installed properly as follows:
```bash
$ strap
```
If you're interested in how much fast this program is, run the benchmark prepared for this program as follows:
```bash
$ git clone https://github.com/heartsh/strap && cd strap
$ cargo test --release -- --nocapture
```

# Author
[Heartsh](https://github.com/heartsh)

# License
Copyright (c) 2018 Heartsh  
Licensed under [the MIT license](http://opensource.org/licenses/MIT).

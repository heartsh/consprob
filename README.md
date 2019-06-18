# RNAfamProb Program, Program for Estimations of Posterior Probabilities on RNA Structural alignment
This project provides the RNAfamProb program, a program for estimations of posterior probabilities on RNA structural alignment.

# Installation
This project has been written in mainly Rust, a systems programming language.
So first, you need to install the Rust compiler, package manager, and standard library.
Visit [the Rust website](https://www.rust-lang.org) to see more about this language.
You can install these 3 components with 1 line as follows:
```bash
$ curl https://sh.rustup.rs -sSf | sh
```
The above installation is done by [Rustup](https://github.com/rust-lang-nursery/rustup.rs), so you can easily switch a compiler in use.
Now you can install the RNAfamProb program as follows: 
```bash
$ cargo install rnafamprob
```
Check if this program has been installed properly as follows:
```bash
$ rnafamprob
```
If you're interested in how much fast this program is, run the benchmark prepared for this program as follows:
```bash
$ git clone https://github.com/heartsh/rnafamprob && cd rnafamprob
$ cargo test --release -- --nocapture
```

# Author
[Heartsh](https://github.com/heartsh)

# License
Copyright (c) 2018 Heartsh  
Licensed under [the MIT license](http://opensource.org/licenses/MIT).

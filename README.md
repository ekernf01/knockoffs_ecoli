## *E. coli* experiments

Code to repeat our experiments on *E. coli* for our manuscript "Model-X knockoffs suggest fundamental limits on regulatory network identification". For more info, see the [central project repo](https://github.com/ekernf01/knockoffs_paper). Anyone should be able to repeat our experiments by using the same environment we used on AWS (Ubuntu 18.04); see `run_on_aws.sh` for details. Plots will be saved to a subdirectory `v28/figures`. Note that figure 2C uses entirely simulated data and is not generated in the ecoli repo but rather the [software demo repo](https://github.com/ekernf01/knockoffs_quick_demo).


```sh
docker pull    ekernf01/knockoffs_ecoli
docker run -it ekernf01/knockoffs_ecoli
```
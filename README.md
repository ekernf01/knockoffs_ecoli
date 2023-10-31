## *E. coli* experiments

Code to repeat our experiments on *E. coli* for our manuscript "Model-X knockoffs suggest fundamental limits on regulatory network identification". For more info, see the [central project repo](https://github.com/ekernf01/knockoffs_paper). Anyone should be able to repeat our experiments using a clean install of Ubuntu 20.04 (e.g. on AWS). You can reproduce results by calling `run_on_aws.sh`. Plots will be saved to a subdirectory `v35/figures`. Note that figure 2C uses entirely simulated data and is not generated in the ecoli repo but rather the [software demo repo](https://github.com/ekernf01/knockoffs_quick_demo).

We also have a docker image for better chance of long-term reproducibility. To start it, run these commands.

```sh
sudo docker pull    ekernf01/knockoffs_ecoli
sudo docker run -it --rm ekernf01/knockoffs_ecoli
```

This will land you in an interactive shell in the Docker container, and you can run our experiments by running [`run_in_docker.sh`](https://github.com/ekernf01/knockoffs_ecoli/blob/main/run_in_docker.sh).

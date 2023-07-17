# opflow
Optical Flow using Dynamic Programming

## Reproducibility using Guix

- Install [Guix](https://guix.gnu.org/en/download/), either as a package
manager on top of a linux system, or as an operating system.

- Set up the compilers using
[Guix time-machine](https://guix.gnu.org/manual/en/html_node/Invoking-guix-time_002dmachine.html)
```
guix time-machine --container channels.scm --environment --manifest manifest.scm 
```

- Run the scripts as indicatd in the [jeif98/README.md](jeif98/README.md)

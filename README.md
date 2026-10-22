# W-state T-count lower bound: exhaustive stabilizer-tensor search

Code accompanying the paper *Heralded W-State Preparation with a Single
Non-Clifford Operation: Protocol and a Four-T Lower Bound* (A. Mitra).

It checks whether any post-selected Clifford circuit fed with `t` magic states
`|T> = T|+>` can output a target state exactly. Every such circuit is described
by a stabilizer state on `t + d` qubits (`d` = number of output qubits), so the
script enumerates **all** stabilizer states on `t + d` qubits, using the
Dehaene–De Moor canonical form, and contracts each one with `|T>^t`.

## Requirements

- Python 3.8+
- NumPy

```
pip install numpy
```

## Usage

Validation runs, and the W-state search for t = 0, 1, 2 (a few seconds):

```
python gadget_search.py
```

Expected output:

```
validation t=0: states=6 solutions=0 best_p=0.0000
validation t=1: states=60 solutions=0 best_p=0.0000
validation t=2: states=1080 solutions=4 best_p=0.7500
W t=0: states=1080 solutions=0 distinct_p=[]
W t=1: states=36720 solutions=0 distinct_p=[]
W t=2: states=2423520 solutions=0 distinct_p=[]
```

The W-state search for t = 3, over all 315,057,600 six-qubit stabilizer states
(about 10 minutes on a laptop):

```
python run_t3.py .
```

Expected output:

```
states 315057600 solutions 0
```

## Result

No post-selected Clifford circuit using three or fewer T gates outputs the
three-qubit W state. Every heralded Clifford+T preparation of the W state
therefore uses at least four T gates per successful run.

## License

MIT (see `LICENSE`).

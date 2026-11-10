# Exact Clifford+T preparation of W states: code (v2.1)

Code for the paper "Exact Clifford+T Preparation of W States: Lower Bounds and
Heralded Constructions" (A. Mitra, 2026). Version 1.0 accompanied the
three-qubit paper. Version 2.0 adds the N-qubit material. Version 2.1 updates
`cost_table.py` to truncate the cheapest exact W_M with M > N, from any route,
as in Eq. (3) of the paper's arXiv v2.

Requires Python 3 and NumPy.

| File | What it does | Runtime |
|---|---|---|
| `gadget_search.py` | Enumerates all stabilizer states on t + d qubits and finds post-selected Clifford gadgets that turn \|T>^t into a target state. Run directly for the W3 search with t <= 2 and the validation cases. | seconds |
| `run_t3.py` | Run as `python3 run_t3.py .`. The t = 3 search for W3 (315,057,600 six-qubit tensors): no solution. | ~10 min |
| `split_search.py` | Searches for exact amplitude splits. `python3 split_search.py split 5` covers 1:2, 1:3, 1:4, 1:8, 1:16, 2:3, 3:4 and 3:5 with t <= 5. `python3 split_search.py w3 3` covers the three-qubit roots 1:2:4 etc. with t <= 3. | ~15 min each |
| `split_adaptive.py` | For one split and t, finds the best success probability counting every Clifford-correctable outcome, and saves the best gadget tensor as JSON. `python3 split_adaptive.py 4 1:8` | seconds (t = 4), ~15 min (t = 5) |
| `gadget_1to8_t4.json`, `gadget_1to2_t2.json` | Saved gadget tensors. The first t qubits are the \|T> inputs and the last qubit is the output. | |
| `cost_table.py` | Computes Table 2 and Fig. 2 of the paper: expected T-counts for N = 3..24 by steps, tree + steps, CH tree, RUS tree and truncation of the cheapest exact W_M (M up to 2N), with the best route in the last column. `python3 cost_table.py 24` | seconds |
| `verify_constructions.py` | Statevector checks of the step proposition (n = 2..10), the 1:8 split circuit and the two-AND W5 circuit. | seconds |

Expected results:

- `split_search.py split 5`: 1:2 first at t = 2, 1:4 and 1:8 first at t = 4, no solutions for 1:3, 2:3, 3:4, 3:5 or 1:16 with t <= 5.
- `split_search.py w3 3`: no solutions for 1:2:4, 1:2:8, 1:4:8, 1:2:16 or 1:4:16 with t <= 3.
- `split_adaptive.py`: 1:2 at t = 2 gives p = 3/4, 1:4 at t = 4 gives p = 5/8 (both match Paetznick and Svore), and 1:8 at t = 4 or 5 gives p = 9/32.
- `cost_table.py 24`: truncation is the best route for 12 of the 22 sizes, e.g. W9 23.6 (from W10), W17 48.5 and W19 43.4 (from W20). The largest ratio to the lower bound max(4, N-1) is about 3.0, at N = 17.

License: MIT.

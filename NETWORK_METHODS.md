# Network construction

## Input and distance

The CLI reads the haplotype table produced by `view -t table`, with metadata
columns `chrom,pos,ref,alt` followed by unique haplotype names. Each entry must
be 0, 1, 2, -1 or missing. A locus with -1 or a missing entry in any haplotype
is excluded from **all** distances. This complete-locus policy keeps the same
denominator and set of loci for every pair. The number excluded is logged, not
added to result columns. No complete loci is an error, including a one-node
input with no observed loci.

Distance is the integer count of unequal genotype states: 0 versus 2 counts
one, as does 0 versus 1. It is not allele dosage distance or an estimate of
evolutionary time. These unphased genotype profiles do not establish a phased
genealogy. Different named profiles can become identical on complete loci;
their nodes and zero-weight edges remain distinct. Complete-locus deletion
may discard much information when missingness is extensive; its result describes
only retained loci, not the discarded positions.

## MST and MSN

MST uses Kruskal's algorithm through NetworkX, with canonical node and edge
ordering to make ties reproducible under input-column permutation. A tree has
n-1 edges when connected; a singleton has one node and zero edges.

MSN is the **union of the edges in every minimum spanning tree**, or of the
minimum spanning forests if a caller supplies a disconnected distance graph.
It retains all equal-weight alternatives, consistent with `include.ties=TRUE`
in poppr's documented MSN concept. This is a distance network; it introduces
no inferred median/ancestral nodes, statistical test or significance level.

Implementation: sort edges by weight. For each equal-weight batch, retain
every edge joining different components formed using strictly lighter edges.
Only then union all endpoints in the batch. An edge whose endpoints are already
joined by lighter edges cannot occur in an MST (cycle property). Every other
edge can be included in a spanning forest of the equal-weight component graph
and hence in an MST. Testing before unioning the batch retains all such ties.
Hamming weights are integers, so exact equality is used, not a numerical tolerance.
Graph construction costs O(H²L); sorting its edges costs O(H² log H), with
O(H²) edge storage. This patch does not address large-network memory scaling.

The removed implementation compared a direct distance to the **sum** along an
MST path times `threshold_factor=1.2`. For Hamming distance, the triangle
inequality makes every pair pass this criterion, producing a complete graph.
The threshold parameter has been removed from both CLI and Python method.
An MSN may still legitimately be complete when all alternatives are tied.

## Output

CSV retains exactly `source,target,weight,type`; ordinary edges use the existing
`type=edge`. Isolated nodes use `type=node` and empty target/weight, without a
self-loop. GEXF/GML retain native nodes. TXT/edgelist retains ordinary tab-separated
edge rows and adds `#node\t<JSON string>` only for isolated nodes. The hastat
reader restores these; generic comment-skipping edge-list readers will not.
No status column or companion summary file is added. `FancyNetwork` uses the
same reader as the export round-trip checks.

## Display

The default group palette matches the existing hastat palette used in Pi and
other plots; sorted group names determine color order. Custom node colors remain
supported. The public example uses opaque sectors so population colors match
the Pi curves exactly. Group and size legends are packed outside the network.
Samples (n) shows up to three actual nonzero node counts: minimum, middle unique
count and maximum. Circle diameters use the same square-root-of-area scaling
as the plotted nodes. Labels are integers and row heights adapt to circle sizes.
When no sample assignments are supplied, unit-sized placeholder nodes do not
produce a sample-count legend. These are display changes, not network inference.

## References and validation

- Kruskal, J. B. (1956). On the shortest spanning subtree of a graph and the
  traveling salesman problem. Proceedings of the AMS 7:48–50.
  DOI: [10.1090/S0002-9939-1956-0078686-7](https://doi.org/10.1090/S0002-9939-1956-0078686-7).
- [NetworkX minimum_spanning_tree, version 3.4](https://networkx.org/documentation/networkx-3.4/reference/algorithms/generated/networkx.algorithms.tree.mst.minimum_spanning_tree.html).
- [poppr.msn documentation, include.ties](https://grunwaldlab.github.io/poppr/reference/poppr.msn.html).
  This supplies the tied-network interpretation; numerical validation uses an
  independent exhaustive oracle, not agreement with another package.

`tests/test_network.py` enumerates every spanning tree for all 729 assignments
of weights {0,1,2} on the six edges of K4 and 30 seeded K5 cases. It compares
the exact minimum-tree edge union and total MST weight, and checks missing
calls, ordering, disconnected graphs, singleton serialization and actual plot
loading. `validate_network.py` also reruns stat/population tests and compares
archived versus current outputs. Synthetic audit files live only in
`qa/20261031_network/`.

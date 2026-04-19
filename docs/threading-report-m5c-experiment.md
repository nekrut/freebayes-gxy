# Threading experiment — M5 Phase C-1 (hypothesis falsified)

**Branch:** `claude/rewrite-freebayes-oL5RQ` @ pending commit
**Result:** negative. Code reverted; only this report ships.

## Hypothesis

Phase A segfaulted because each Rayon worker called
`bam::IndexedReader::from_path` concurrently. The threading report for
Phase B (`threading-report-m5b.md`) catalogued four workarounds, of
which #1 was "sequential index load + parallel `.fetch`":

> If the hts_sys race is only at index-load time, pre-opening N
> readers on the main thread (sequentially) and handing each to a
> worker should sidestep it.

Cheap test; potentially unlocks real wall-clock speedup.

## Experiment

Rewrote `run_call_parallel` to pre-open one `bam::IndexedReader` per
tile on the main thread, zip them with the tile list, and iterate via
`rayon::into_par_iter` so each worker owns one ready reader and only
needs to call `.fetch` + `.read`.

Relevant code shape:

```rust
// Main thread: sequential opens.
let mut readers = Vec::with_capacity(tiles.len());
for _ in 0..tiles.len() {
    readers.push(bam::IndexedReader::from_path(&cli.bam)?);
}

// Workers: owned reader + tile.
tiles
    .into_iter()
    .zip(readers)
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|((tid, window), mut reader)| {
        reader.fetch((tid, window.start as i64, window.end as i64))?;
        // ... walk + clump + call ...
    })
    .collect();
```

## Observed

On the 100 kb / 19971-read fixture at `--threads 4 --tile-size
10000` (10 tiles, 4 workers):

```
run 1 exit=139  (SIGSEGV)
run 2 exit=139
run 3 exit=139
```

Reproducible across runs. Only 2 records in the output file before
crash, suggesting the crash happens during concurrent `.fetch` or
`.read`.

## Conclusion

**The race is deeper than `from_path`.** Pre-opening readers does
not sidestep it. Suspect culprit: the first `.fetch` call still
loads (parts of) the BAM index lazily, or `hts_sys`'s BGZF decoder
shares mutable state across reader instances. Either way, option 1
from the Phase B report is falsified.

## Remaining workarounds (from threading-report-m5b.md)

2. **Single-reader producer + parallel-consumer walk** — one thread
   reads records from the BAM into a channel, N workers consume
   records and run walk + clump. Depends on `walk_record` +
   `clump_observations` being the next biggest slice of runtime
   after `bam::Reader::read`. Likely viable; adds channel overhead
   and a lock-free pileup aggregator.
3. **Pre-shard BAM into per-region files offline.** Expensive for
   large BAMs; workable as a pipeline-level decision.
4. **Switch to `noodles`.** Pure Rust, designed for Send-across-
   threads BAM readers. Larger dep change; probably the right long
   term.

## What lands in this commit

Only this report. `run_call_parallel` stays at the Phase B design
(single-threaded pileup + parallel per-site call), as committed in
`44eda60`.

## Next test to run (if someone picks this up)

Option 2 (single-reader producer + parallel-consumer walk):
- Main thread: read records, push `bam::Record` clones to a
  crossbeam channel.
- N workers: consume records, run `walk_record` + `clump_observations`,
  emit observations into a `DashMap<(i32, i64), Vec<AlleleObs>>`.
- Main thread (after channel close): iterate the dashmap into
  `Pileup`, run the existing parallel call stage.

Profiling showed `walk_record` + `clump` likely account for a large
fraction of the "pileup" stage's 2.6s (we don't have a finer
breakdown). If ≥50% of pileup time is walk+clump, this unlocks
measurable speedup.

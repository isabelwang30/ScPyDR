# Benchmarking scPyDR against scanpy's pca and umap functions

## Timing

```
time scpydr benchmark/data -u -v

time [WRITE COMMAND TO RUN SCANPY PCA ON BENCHMARK]
```

## Memory usage

See: https://github.com/jhclark/memusg

```
memusg scpydr benchmark/data -u -v

memusg [WRITE COMMAND TO RUN SCANPY PCA ON BENCHMARK]
```
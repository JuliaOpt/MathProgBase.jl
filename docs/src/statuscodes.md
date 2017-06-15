```@meta
CurrentModule = MathProgBase
```


## Termination Status

The `TerminationStatus` attribute is meant to explain the reason why the solver stopped executing. The value of the attribute is of type `TerminationStatusCode`.

```@docs
TerminationStatus
TerminationStatusCode
```

## Result Status

The `PrimalStatus` and `DualStatus` attributes are meant to explain how to interpret the result returned by the solver. The value of the attributes are of type `ResultStatusCode`.

```@docs
PrimalStatus
DualStatus
ResultStatusCode
```

## Basis Status

TODO: attributes and status codes for LP basis status


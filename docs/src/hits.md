# Hits

The hit type hierarchy is currently being redesigned!

`NeRCA.jl` has three main types of hits: `AbstractDAQHit`, `AbstractMCHit` and
`CalibratedHit`.

| hit type                        | supertype        | guaranteed attributes                              |
|---------------------------------|------------------|----------------------------------------------------|
| `Hit` (in future `SnapshotHit`) | `AbstractDAQHit` | `channel_id`, `dom_id`, `t`, `tot`                 |
| `TriggeredHit`                  | `AbstractDAQHit` | `channel_id`, `dom_id`, `t`, `tot`, `trigger_mask` |
| `MCHit`                         | `AbstractMCHit`  | `a`, `origin`, `pmt_id`, `t`                       |

```@autodocs
Modules = [NeRCA]
Pages = ["hits.jl"]
```

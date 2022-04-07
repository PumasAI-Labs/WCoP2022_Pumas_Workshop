# Bioequivalence in Pumas

You should have in the folder two files:

- `be.jl` - Pumas code in Julia
- `data4be.csv` - dataset

## Data Dictionary for `data4be.csv`

|Variable|Label|Type|Code|
|--------|-----|----|----|
|id|SubjectID|NUM||
|sequence|Sequence|CHAR|RT – **R**eference **T**est & TR – **T**est **R**eference|
|period|Study period|NUM||
|AUC|Area under the curve|NUM||
|Cmax|Maximum concentration|NUM||

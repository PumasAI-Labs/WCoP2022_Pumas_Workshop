# NCA in Pumas

You should have in the folder two files:

- `nca_01.jl` - Pumas code in Julia
- `nca_02.jl` - Pumas code in Julia
- `data4be.csv` - dataset

## Data Dictionary for `data4NCA_sad.csv`

| Variable      | Label                       | Type  | Code                |
|---------------|-----------------------------|-------|---------------------|
| id            | Subject ID                  | NUM   |                     |
| time          | Time in hours               | NUM   |                     |
| concentration | CTMX concentration in μg/mL | NUM   |                     |
| amt           | CTMX dose in mg             | NUM   |                     |
| evid          | Dosing event                | NUM   | 0 = NO 1 = YES      |
| cmt           | Compartment                 | NUM   |                     |
| rate          | Dosing rate                 | NUM   |                     |
| age           | Age in years                | NUM   |                     |
| wt            | Weight in kg                | NUM   |                     |
| doselevel     | CTMX Dose in mg             | NUM   |                     |
| isPM          | Is a poor metaboliser       | CHAR  |                     |
| isfed         | Fed condition               | CHAR  |                     |
| sex           | sex                         | CHAR  |                     |
| route         | Route of administration     | CHAR  | ev = extravascular  |

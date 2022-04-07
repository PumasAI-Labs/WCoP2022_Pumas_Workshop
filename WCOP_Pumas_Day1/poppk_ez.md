# Population Pharmacokinetic (Pop PK) in Pumas

You should have in the folder two files:

- `poppk_iv.jl` - Pumas code in Julia
- `data4poppk_iv.csv` - dataset

## Data Dictionary for `data4poppk_iv.csv`

| Variable | Label                     | Type | Code                |
|----------|---------------------------|------|---------------------|
| ID       | Subject ID                | NUM  |                     |
| TIME     | Time in hrs               | NUM  |                     |
| CONC     | concentration             | NUM  |                     |
| AMT      | Dose in mg                | NUM  |                     |
| AGE      | Age in yrs                | NUM  |                     |
| WEIGHT   | Weight in kg              | NUM  |                     |
| SCR      | Serum Creatinine          | NUM  |                     |
| ISMALE   | Gender                    | NUM  | 0 = MALE 1 = FEMALE |
| eGFR     | Glomular filteration rate | NUM  |                     |
| EVID     | Dosing event              | NUM  | 0 = NO 1 = YES      |
| CMT      | Compartment               | NUM  |                     |

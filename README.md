# Model class
Includes:
 - A
 - Alpha

# Model family
Includes:
 - A00
 - A01
 - Alpha00
 - Alpha01

# Model numbers
 - Alpha01/Mdl00 (Prereged RSA)
 - Alpha01/Mdl01 (Prereged long-axis model [q])
 - Alpha01/Mdl02 (Saturated long-axis model [q])
 - Alpha01/Mdl03 (Saturated long-axis model [epsilon])

# Alpha00
This folder contains analysis results for Least Square Seperate (LSS) models using different modeling approaches and configurations.

## Directory Structure

The analysis results are organized in a hierarchical structure based on model identifiers, position types, integer values, and run numbers:

``` text
Analysis/
├── Alpha**00**/
│   ├── i0/
│   ├── i1/
│   ├── ...
│   ├── i5/
└── Alpha**01**/
    ├── a0/
    ├── a1/
    ├── ...
    ├── a5/
    ├── b0/
    ├── ...
    ├── b4/
    └── b5/
```

## Naming Convention

### Model Identifiers

- **Alpha00**: Native space no-ordering LSS models
  - Does not distinguish between a and b positions
  - Uses integer notation: `ix_Ry` where `x` is the integer value (0-5) and `y` is the run number

- **Alpha01**: Native space ordering LSS models
  - Distinguishes between a and b positions
  - Uses position notation: `[a|b]x_Ry` where position is either `a` or `b`, `x` is the integer value (0-5), and `y` is the run number
    Within **Alpha01**, several models are implemented:
      - `mdl00` is the preregistered mixed-effects model applied to the Fisher z-transformed statistic.
      - `mdl01` is the preregistered axis model.
      - `mdl02` is above model with interaction of colocation
      - `mdl03` is axis analysis using epsilon term
      - `mdl04` is a fully saturated analysis of the zTemplate searchlight 
      - `mdl05` looks at the simple effects highlighted by mdl04 
      - scripts `z0-02` generate the t-patterns used by all models 
      - scripts `z03-04` generate the searchlight used by mdl04-5
      - scripts within the model folders are numbered carrying on from the scripts which they depend on in the above folder  

### Directory Naming Pattern

```
Analysis/Alpha**MM**/[position]X_RY/
```

Where:
- `MM`: Two-digit model identifier (00, 01, etc.)
- `position`: Position identifier (`i` for Alpha00 models, `a` or `b` for Alpha01+ models)
- `X`: Integer value (0-5)
- `Y`: Run number (1, 2, 4, 5, etc.)

## Model Types

### Alpha00 Models
- **Approach**: Native space no-ordering
- **Position Handling**: No distinction between positions
- **Directory Format**: `Analysis/Alpha**00**/iX_RY/`
- **Example**: `Analysis/Alpha**00**/i3_R2/`

### Alpha01 Models
- **Approach**: Native space ordering
- **Position Handling**: Distinguishes between a and b positions
- **Directory Format**: `Analysis/Alpha**01**/[a|b]X_RY/`
- **Examples**: 
  - `Analysis/Alpha**01**/a2_R1/`
  - `Analysis/Alpha**01**/b4_R3/`

### Model Specification Files
  Each model is specified using an .mlx file. Decimal suffixes are used to indicate progressively more complex versions of a given base model, reflecting the incremental addition of terms.

  For example:
  - `mdl00.0` denotes the basic version of the preregistered mixed-effects model.
  - `mdl00.1` denotes an extended version of mdl00 with additional interaction terms.<br>
  Further decimal increments (e.g. `mdl00.2`) indicate subsequent model extensions.<br>
  This decimal-based convention applies consistently across all model families (e.g. mdl01, mdl02).

## Usage

Each directory contains the analysis results for a specific model configuration and run. The modeling approaches can be distinguished by their Alpha identifiers:

- Use `Alpha**00**` for analyses that don't require position ordering
- Use `Alpha**01**` for analyses that incorporate position ordering effects
- Additional Alpha identifiers can be added for future modeling approaches

## Data Organization

- Integer values range from 0 to 5
- Run numbers vary based on the specific analysis (commonly R1, R2, R4, R5)
- Each sub-directory contains the complete analysis output for that specific configuration

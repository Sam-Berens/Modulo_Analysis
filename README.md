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

## Usage

Each directory contains the analysis results for a specific model configuration and run. The modeling approaches can be distinguished by their Alpha identifiers:

- Use `Alpha**00**` for analyses that don't require position ordering
- Use `Alpha**01**` for analyses that incorporate position ordering effects
- Additional Alpha identifiers can be added for future modeling approaches

## Data Organization

- Integer values range from 0 to 5
- Run numbers vary based on the specific analysis (commonly R1, R2, R4, R5)
- Each sub-directory contains the complete analysis output for that specific configuration

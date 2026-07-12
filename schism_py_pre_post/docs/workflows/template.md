# Workflow title

## Purpose

State what the workflow accomplishes and when it should be used.

## Prerequisites

- Required Python environment and external programs.
- Required credentials or API keys.
- Required model configuration and data.

## Inputs

| Input | Format | Description |
|---|---|---|
| `path/to/input` | Format | What this input contains. |

## Procedure

### 1. First operation

Explain the operation, then provide a complete command:

```bash
python path/to/script.py \
  --input path/to/input \
  --output-dir path/to/output
```

### 2. Next operation

Continue with the next reproducible step.

## Outputs

| Output | Description |
|---|---|
| `output/example` | What downstream workflow consumes it. |

## Validation

Provide fast checks that catch incomplete or incorrectly ordered data:

```bash
test -s path/to/output
```

## Assumptions and limitations

- Record time-zone, datum, units, depth, and interpolation assumptions.
- State whether rerunning is safe and whether the output is resumable.

## Troubleshooting

Describe workflow-specific problems and link to
[general troubleshooting](../troubleshooting.md).


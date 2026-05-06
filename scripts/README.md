# MATLAB Parity Scripts

This folder contains scripts to assist in achieving exact parity between MATLAB and Python.

## 🔍 Execution Tracing (`extract_watershed_trace.m`)

This script (and the logic within) is used to generate JSONL execution traces from MATLAB that can be compared with Python traces using `dev/scripts/analysis/compare_execution_traces.py`.

### Integration into `get_edges_by_watershed.m`

To debug a full run, it is recommended to add the following snippets directly to the MATLAB source:

1.  **Open file at function start**:
    ```matlab
    fid = fopen('matlab_trace.jsonl', 'w');
    ```

2.  **Record iteration start** (inside `while ~isempty(available_locations)`):
    ```matlab
    event = struct('event', 'iteration_start', ...
                   'iteration', iteration, ...
                   'current_linear', current_location, ...
                   'current_energy', energy_map_temp(current_location));
    fprintf(fid, '%s\n', jsonencode(event));
    ```

3.  **Record seed selection** (inside `for seed_idx = seed_index_range`):
    ```matlab
    event = struct('event', 'seed_selected', ...
                   'seed_idx', seed_idx, ...
                   'selected_linear', next_location, ...
                   'selected_energy', current_strel_energies(strel_idx));
    fprintf(fid, '%s\n', jsonencode(event));
    ```

4.  **Close file at function end**:
    ```matlab
    fclose(fid);
    ```

### Comparison

Once you have `matlab_trace.jsonl` and `python_trace.jsonl`, use the Python comparison tool:
```powershell
python dev/scripts/analysis/compare_execution_traces.py matlab_trace.jsonl python_trace.jsonl
```

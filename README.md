## Project Overview
This research project focuses on a **parallel implementation** of proximity calculation for a set of points over multiple time steps. The goal is to identify groups of points that satisfy specific proximity criteria, making the problem computationally intensive due to the need to check distances between all point pairs at each time step.

### Problem Definition
Given a set of points in a 2D plane, a point is said to satisfy the **Proximity Criteria** if there exist at least `K` points with a distance less than a specified value `D` from it. The task is to evaluate this criterion for `tCount + 1` values of parameter `t`:

\[ t = \frac{2 * i}{\text{tCount}} - 1 \quad \text{where } i = 0, 1, 2, ..., \text{tCount} \]

The computation must identify at least 3 points satisfying the criteria for each `t` value, and if found, further evaluation for that specific `t` should stop.

### Input and Output Specifications
1. **Input File (`input.txt`)**:
   - **Line 1**: Contains four parameters - `N` (number of points), `K` (minimum number of points satisfying the criteria), `D` (distance threshold), and `TCount`.
   - **Next N lines**: Each line contains the parameters for each point: `id`, `x1`, `x2`, `a`, `b`.

   **Example Input:**
   ```
   4 2 1.23 100
   0 2.2 1.2 2 45.07
   1 -1 26.2 44 -3.3
   2 -43.3 12.2 4.7 20
   3 11.0 -6.6 12.5 23
   ```

2. **Output File (`output.txt`)**:
   - For each `t` value where 3 points satisfy the criteria, a line in the format:
     ```
     Points pointID1 pointID2 pointID3 satisfy Proximity Criteria at t = t1
     ```
   - If no such points are found, the output will be:
     ```
     There were no 3 points found for any t.
     ```

## Parallelization Approach
To optimize the computational process, the solution leverages a combination of **MPI**, **OpenMP**, and **CUDA** for efficient parallel computation:

### MPI (Message Passing Interface)
- **Purpose**: Distribute time step ranges across multiple processes.
- **Rationale**: Time steps are independent, making them ideal for distribution.
- **Architecture**: Implement a **Master-Worker Model**:
  - **Master Process**: Distributes the range of `t` values to each worker.
  - **Worker Processes**: Compute the proximity criteria in parallel for assigned `t` ranges and send results back to the master.

### OpenMP (Multi-threading)
- **Purpose**: Parallelize point distance calculations within each MPI process.
- **Rationale**: Each point pair can be processed independently across threads, minimizing execution time.

### CUDA (GPU Acceleration)
- **Purpose**: Offload the most computationally expensive parts, specifically the **X** and **Y** coordinate calculations, to the GPU.
- **Rationale**: GPUs excel at parallel arithmetic operations, significantly improving performance.

### Performance Consideration
The parallel solution should outperform the sequential implementation, demonstrating efficiency on a distributed system like **VLAB** (two computers from different pools when using MPI).

## How to Run the Project
1. Compile the code using the appropriate compilers for MPI, OpenMP, and CUDA.
2. Run the job using the following command:

   ```bash
   sbatch sbatch.sbatch
   ```

3. Ensure the `input.txt` file is present in the same directory as the executable.

## Expected Output
The output will list time steps and point IDs that satisfy the proximity criteria:

**Example Output:**
```
Points 8, 18, 9 satisfy Proximity Criteria at t = -1.000000
Points 8, 18, 9 satisfy Proximity Criteria at t = -0.993355
...
Points 8, 18, 9 satisfy Proximity Criteria at t = 1.000000
```

If no points meet the criteria for any `t`, the output will state:

```
There were no 3 points found for any t.
```

## Implementation Details
- The code uses **MPI** to distribute the range of `t` values across multiple processes.
- **OpenMP** is used to speed up distance calculations within each process.
- **CUDA** is used to accelerate coordinate calculations on the GPU.

### Key Equations
The x and y coordinates for each point are calculated using the following equations:

\[
x = \left( \frac{x2 - x1}{2} \right) * \sin \left( \frac{t * \pi}{2} \right) + \frac{x2 + x1}{2}
\]
\[
y = a * x + b
\]

### Load Balancing Consideration
Effective load balancing strategies are employed to ensure that the computational workload is evenly distributed across processes and threads, minimizing idle time and maximizing resource utilization.

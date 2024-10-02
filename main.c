#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include "myProto.h"

// Calculates the X corr
double calc_x(double x1, double x2, double t);

// Calculates the Y corr
double calc_y(double x, double a, double b);

// Calculates distance between (x1, y1) , (x2, y2)
double calc_distance(double x1, double y1, double x2, double y2);

// Processes points to check  if there is  groups of points that satisfy the proximity in each t value and  save the result
void process_points(Point points[], int N, int K, double D, int t_start, int t_end, char *result_buffer, int *result_length, int t_count);

// Read the file input data , (points , N, K, D, , tCount)
int read_file(const char *filename, Point **points, int *N, int *K, double *D, int *tCount);

// Distributes the work for the prossess acc to the size of the work for each prossess.
void distribute_job(int world_size, int tCount);

// getting  results from all processes, and putting them in a single buffer .
void get_results(int world_size, char *result_buffer, int *result_length);

// Write the result buffer to  output file
void write_output(char *result_buffer, int result_length);

// the processing receiving task ranges of worker process and and sending results back to the boss prossess after working on their tasks
void process_worker(int N, int K, double D, int tCount, Point *points);

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    Point *points;
    int N, K, tCount;
    double D;

    if (world_size < 2) {
        if (world_rank == 0) {
                fprintf(stderr, "Error: The program requires at least 3 processes.\n");
            }
            MPI_Abort(MPI_COMM_WORLD, 1); // Exit the MPI environment with error code 1
        return 1;
    }

    if (world_rank == 0 && !read_file("input.txt", &points, &N, &K, &D, &tCount))
    {
        printf("Error while reading the file.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&K, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&D, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (world_rank != 0)
        points = (Point *)malloc(N * sizeof(Point)); // create struct point
    MPI_Bcast(points, N * sizeof(Point), MPI_BYTE, 0, MPI_COMM_WORLD);

    if (world_rank == 0) // im the boss prossess
    {
        distribute_job(world_size, tCount); // sending jobs
        char result_buffer[1024 * 256];
        int result_length = 0;
        get_results(world_size, result_buffer, &result_length); // getting the result of the jobs from workers
        write_output(result_buffer, result_length);             // write the result to output file
    }
    else // im a prossesss worker
    {
        process_worker(N, K, D, tCount, points);
    }
    free(points);
    MPI_Finalize();
    return 0;
}

void distribute_job(int world_size, int tCount)
{
    int num_processes = world_size - 1; // Num of worker processes
    int t_per_process = tCount / num_processes;
    int rem = tCount % num_processes;

    for (int i = 1; i < world_size; i++)
    {
        int extra = (i <= rem) ? 1 : 0;
        int t_start = (i - 1) * t_per_process + (i - 1 < rem ? i - 1 : rem);
        int t_end = t_start + t_per_process + extra - 1;

        // make sure thate last process arrive to and prossess the  tCount
        if (i == world_size - 1)
        {
            t_end = tCount;
        }
        // Debugging
        // printf("Process %d: t_start = %d, t_end = %d\n", i, t_start, t_end);
        MPI_Send(&t_start, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        MPI_Send(&t_end, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
}

void get_results(int world_size, char *result_buffer, int *result_length)
{
    for (int i = 1; i < world_size; i++)
    {                                                                                  // loop on all prossess
        int worker_reslen;                                                             // the length of the result od the worker
        MPI_Recv(&worker_reslen, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // getting the result from workers
        MPI_Recv(result_buffer + *result_length, worker_reslen, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        *result_length += worker_reslen;
    }
}

void write_output(char *result_buffer, int result_length)
{
    MPI_File output;
    int file_exists;

    // Try to delete the file if it exists
    MPI_File_delete("output.txt", MPI_INFO_NULL);

    // Create and open a new file for writing
    MPI_File_open(MPI_COMM_SELF, "output.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &output);

    // Check if the first character of result_buffer is '\0', meaning no points were found
    if (result_buffer[0] == '\0')
    {
        const char *no_points_message = "There were no 3 points found for any t.\n";
        MPI_File_write(output, no_points_message, strlen(no_points_message), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    else
    {
        // Write the result to the output file
        MPI_File_write(output, result_buffer, result_length, MPI_CHAR, MPI_STATUS_IGNORE);
    }

    // Close the file
    MPI_File_close(&output);
}

void process_worker(int N, int K, double D, int tCount, Point *points)
{
    int t_start, t_end;
    MPI_Recv(&t_start, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recvining time of the start
    MPI_Recv(&t_end, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   // recvining time of the end
    char result_buffer[1024 * 256];
    int result_length = 0;
    process_points(points, N, K, D, t_start, t_end, result_buffer, &result_length, tCount); // prossess on points acc each t between t stat and t end
    MPI_Send(&result_length, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);                             // sending the res length and the the result
    MPI_Send(result_buffer, result_length, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
}

void process_points(Point points[], int N, int K, double D, int t_start, int t_end, char *result_buffer, int *result_length, int t_count)
{

    int num_iterations = t_end - t_start + 1;
    int temp_found_ids[num_iterations][3];
    int temp_found_count[num_iterations];

    char local_buffer[100 * num_iterations];
    local_buffer[0] = '\0'; // Ensure it's empty
    int offset = 0;

#pragma omp parallel for
    for (int i = t_start; i <= t_end; i++)
    {
        double t = (double)(2 * (double)i / (double)t_count) - 1;
        double x_val[N], y_val[N];
        int Id_points[3];
        int found = 0;

        ComputeXYValues(points, x_val, y_val, N, t);

    #pragma omp parallel for shared(found, Id_points)
        for (int j = 0; j < N; j++)
        {
            if (found >= 3)
                continue; // Skip if already found 3 points (non-breaking due to parallelism)

            int count = 0;
            for (int k = 0; k < N; k++)
            {
                if (j != k && calc_distance(x_val[j], y_val[j], x_val[k], y_val[k]) < D)
                {
                    count++;
                    if (count >= K)
                    {
                        // this critical section only for the inner loop to avoid race condition between the threads
                        #pragma omp critical
                        {
                            if (found < 3) // Re-check inside critical section
                                Id_points[found++] = points[j].id;
                        }
                        break;
                    }
                }
            }
        }

        if (found >= 3)
        {
            temp_found_ids[i - t_start][0] = Id_points[0];
            temp_found_ids[i - t_start][1] = Id_points[1];
            temp_found_ids[i - t_start][2] = Id_points[2];
            temp_found_count[i - t_start] = found;
        }
        else
        {
            temp_found_count[i - t_start] = 0;
        }
    }

    for (int i = 0; i < num_iterations; i++)
    {
        if (temp_found_count[i] >= 3)
        {
            double t = (double)(2 * (double)(i + t_start) / (double)t_count) - 1;
            int length = snprintf(local_buffer + offset, sizeof(local_buffer) - offset,
                                  "Points %d, %d, %d satisfy Proximity Criteria at t = %f\n",
                                  temp_found_ids[i][0], temp_found_ids[i][1], temp_found_ids[i][2], t);
            offset += length;
        }
    }

    memcpy(result_buffer, local_buffer, offset);
    *result_length = offset;
}

double calc_x(double x1, double x2, double t)
{
    return ((x2 - x1) / 2.0) * sin(t * 3.14 / 2.0) + (x2 + x1) / 2.0;
}

double calc_y(double x, double a, double b)
{
    return a * x + b;
}

double calc_distance(double x1, double y1, double x2, double y2)
{
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

int read_file(const char *filename, Point **points, int *N, int *K, double *D, int *tCount)
{
    FILE *input = fopen(filename, "r");
    if (input == NULL)
    {
        return 0;
    }

    fscanf(input, "%d %d %lf %d", N, K, D, tCount);

    *points = (Point *)malloc(*N * sizeof(Point));
    if (*points == NULL)
    {
        fclose(input);
        return 0;
    }

#pragma omp parallel for
    for (int i = 0; i < *N; i++)
    {
        fscanf(input, "%d %lf %lf %lf %lf", &(*points)[i].id, &(*points)[i].x1, &(*points)[i].x2, &(*points)[i].a, &(*points)[i].b);
        printf("%d %lf %lf %lf %lf\n", (*points)[i].id, (*points)[i].x1, (*points)[i].x2, (*points)[i].a, (*points)[i].b);
    }

    fclose(input);
    return 1;
}
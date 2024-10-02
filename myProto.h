#ifndef MYPROTO_H
#define MYPROTO_H

typedef struct {
    int id;
    double x1;
    double x2;
    double a;
    double b;
} Point;

void ComputeXYValues(Point* points, double* x_vals, double* y_vals, int N, double t);

#endif // MYPROTO_H
#include <iostream>
#include <cmath>
#include <limits>
#include <chrono>
#include <thread>
using namespace std;

static int g_threadCount = 1;
struct Point {
    double x, y;
};
double euclidean_distance(const Point &a, const Point &b) {
    double dx = a.x - b.x,
           dy = a.y - b.y;
    return sqrt(dx*dx + dy*dy);
}
void compute_dist_chunk(double** dist, Point* pts, int n, int startRow, int endRow) {
    for (int i = startRow; i < endRow; ++i) {
        dist[i][i] = 0.0;
        for (int j = i+1; j < n; ++j) {
            double d = euclidean_distance(pts[i], pts[j]);
            dist[i][j] = dist[j][i] = d;
        }
    }
}
double** compute_distance_matrix(Point* points, int n) {
    double** dist = new double*[n];
    for (int i = 0; i < n; ++i)
        dist[i] = new double[n];
    int T = min(g_threadCount, n);
    thread* thr = new thread[T];
    for (int t = 0; t < T; ++t) {
        int start = (n * t) / T;
        int end   = (n * (t+1)) / T;
        thr[t] = thread(compute_dist_chunk,
                        dist, points, n, start, end);
    }
    for (int t = 0; t < T; ++t)
        thr[t].join();
    delete[] thr;
    return dist;
}

int* greedy_tour(double** dist, int n) {
    bool* visited = new bool[n]();
    int* tour = new int[n];
    tour[0] = 0; visited[0] = true;
    int current = 0;
    for (int i = 1; i < n; ++i) {
        int next = -1;
        double best = numeric_limits<double>::max();
        for (int j = 0; j < n; ++j) {
            if (!visited[j] && dist[current][j] < best) {
                best = dist[current][j];
                next = j;
            }
        }
        tour[i] = next; visited[next] = true;
        current = next;
    }
    delete[] visited;
    return tour;
}
double compute_length(int* tour, double** dist, int n, bool cycle) {
    double L = 0;
    for (int i = 0; i+1 < n; ++i)
        L += dist[tour[i]][tour[i+1]];
    if (cycle) L += dist[tour[n-1]][tour[0]];
    return L;
}
void reverse_segment(int* a, int i, int j) {
    while (i < j) {
        swap(a[i], a[j]);
        ++i; --j;
    }
}
bool two_opt_swap(int* tour, double** dist, int n) {
    double oldL = compute_length(tour, dist, n, true);
    for (int i = 0; i+1 < n; ++i) {
        for (int k = i+1; k < n; ++k) {
            int* tmp = new int[n];
            memcpy(tmp, tour, n * sizeof(int));
            reverse_segment(tmp, i+1, k);
            double newL = compute_length(tmp, dist, n, true);
            if (newL < oldL) {
                memcpy(tour, tmp, n * sizeof(int));
                delete[] tmp;
                return true;
            }
            delete[] tmp;
        }
    }
    return false;
}
bool three_opt_swap(int* tour, double** dist, int n) {
    double oldL = compute_length(tour, dist, n, true);
    for (int i = 0; i < n-2; ++i) {
        for (int j = i+1; j < n-1; ++j) {
            for (int k = j+1; k < n; ++k) {
                for (int c = 1; c <= 4; ++c) {
                    int* tmp = new int[n];
                    int idx=0;
                    for (int m = 0; m <= i; ++m) tmp[idx++]=tour[m];
                    if (c==1) {
                      for (int m=j; m>i; --m) tmp[idx++]=tour[m];
                      for (int m=j+1; m<=k; ++m) tmp[idx++]=tour[m];
                    } else if (c==2) {
                      for (int m=i+1; m<=j; ++m) tmp[idx++]=tour[m];
                      for (int m=k; m>j; --m) tmp[idx++]=tour[m];
                    } else if (c==3) {
                      for (int m=j; m>i; --m) tmp[idx++]=tour[m];
                      for (int m=k; m>j; --m) tmp[idx++]=tour[m];
                    } else {
                      for (int m=j+1; m<=k; ++m) tmp[idx++]=tour[m];
                      for (int m=i+1; m<=j; ++m) tmp[idx++]=tour[m];
                    }
                    for (int m=k+1; m<n; ++m) tmp[idx++]=tour[m];
                    double newL = compute_length(tmp, dist, n, true);
                    if (newL < oldL) {
                        memcpy(tour, tmp, n*sizeof(int));
                        delete[] tmp;
                        return true;
                    }
                    delete[] tmp;
                }
            }
        }
    }
    return false;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout<<"Threads for distance matrix: ";
    cin>>g_threadCount;
    int n; cout<<"#points: "; cin>>n;
    Point* pts = new Point[n];
    for(int i=0;i<n;++i) cin>>pts[i].x>>pts[i].y;
    int mode; cout<<"0=path,1=cycle: "; cin>>mode;
    bool cycle=(mode==1);
    double** dist = compute_distance_matrix(pts,n);
    int* tour = greedy_tour(dist,n);
    cout<<"\nGreedy tour length="<<compute_length(tour,dist,n,cycle)<<"\n";
    while(two_opt_swap(tour,dist,n)){}
    cout<<"After 2-opt: "<<compute_length(tour,dist,n,cycle)<<"\n";
    while(three_opt_swap(tour,dist,n)){}
    cout<<"After 3-opt: "<<compute_length(tour,dist,n,cycle)<<"\n";
    delete[] tour;
    delete[] pts;
    for(int i=0;i<n;++i) delete[] dist[i];
    delete[] dist;
}

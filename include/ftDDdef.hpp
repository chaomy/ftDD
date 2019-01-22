#ifndef _pfDDdef_
#define _pfDDdef_

#include <string>
#include <unordered_map>
#include <vector>

#define PI 3.14159265359
#define INVPI 0.31830988618
#define Mc 0.1

using std::string;
using std::unordered_map;
using std::vector;

class DDdat;

// fitting parameters
struct FitVariable {
  double err, mob, alp, bta;
};

// Mat Constants
struct MatConstant {
  double mu, burg, mub, taylor;
};

class slip {
 public:
  string nm;
  vector<double> b;    // burgers
  vector<double> n;    // plane normal
  vector<double> dns;  // density
  vector<double> flx;  // flux
  slip(string nm, vector<double> b, vector<double> n) : nm(nm), b(b), n(n) {}
};

class DDdat {
 public:
  vector<slip> slips;
};

class DDfcc : public DDdat {
 public:
  DDfcc() {
    // For fluxtot_b1: burgers vector 1/2 [ 1  1  0]
    // 4         ( 1 -1  1), edge
    // 5         (-1  1  1), edge
    // 6         ( 1 -1  1), screw
    // 7         (-1  1  1), screw
    slips.push_back(slip("[110](1-11)e", {1, 1, 0}, {1, -1, 1}));
    slips.push_back(slip("[110](-111)e", {1, 1, 0}, {-1, 1, 1}));
    slips.push_back(slip("[110](-111)s", {1, 1, 0}, {1, -1, 1}));
    slips.push_back(slip("[110](-111)s", {1, 1, 0}, {-1, 1, 1}));

    // For fluxtot_b2: burgers vector 1/2 [ 1 -1  0]
    // 4         ( 1  1  1), edge
    // 5         ( 1  1 -1), edge
    // 6         ( 1  1  1), screw
    // 7         ( 1  1 -1), screw
    slips.push_back(slip("[1-10](111)e", {1, -1, 0}, {1, 1, 1}));
    slips.push_back(slip("[1-10](11-1)e", {1, -1, 0}, {1, 1, -1}));
    slips.push_back(slip("[1-10](111)s", {1, -1, 0}, {1, 1, 1}));
    slips.push_back(slip("[1-10](11-1)s", {1, -1, 0}, {1, 1, -1}));

    // For fluxtot_b3: burgers vector 1/2 [ 1  0  1]
    // 4         ( 1  1 -1), edge
    // 5         (-1  1  1), edge
    // 6         ( 1  1 -1), screw
    // 7         (-1  1  1), screw
    slips.push_back(slip("[101](11-1)e", {1, 0, 1}, {1, 1, -1}));
    slips.push_back(slip("[101](-111)e", {1, 0, 1}, {-1, 1, 1}));
    slips.push_back(slip("[101](11-1)s", {1, 0, 1}, {1, 1, -1}));
    slips.push_back(slip("[101](-111)s", {1, 0, 1}, {-1, 1, 1}));

    // For fluxtot_b4: burgers vector 1/2 [ 1  0 -1]
    // 4         ( 1  1  1), edge
    // 5         ( 1 -1  1), edge
    // 6         ( 1  1  1), screw
    // 7         ( 1 -1  1), screw
    slips.push_back(slip("[10-1](111)e", {1, 0, -1}, {1, 1, 1}));
    slips.push_back(slip("[10-1](1-11)e", {1, 0, -1}, {1, -1, 1}));
    slips.push_back(slip("[10-1](111)s", {1, 0, -1}, {1, 1, 1}));
    slips.push_back(slip("[10-1](1-11)s", {1, 0, -1}, {1, -1, 1}));

    // For fluxtot_b5: burgers vector 1/2 [ 0  1  1]
    // 4         ( 1  1 -1), edge
    // 5         ( 1 -1  1), edge
    // 6         ( 1  1 -1), screw
    // 7         ( 1 -1  1), screw
    slips.push_back(slip("[011](11-1)e", {0, 1, 1}, {1, 1, -1}));
    slips.push_back(slip("[011](1-11)e", {0, 1, 1}, {1, -1, 1}));
    slips.push_back(slip("[011](11-1)s", {0, 1, 1}, {1, 1, -1}));
    slips.push_back(slip("[011](1-11)s", {0, 1, 1}, {1, -1, 1}));

    // For fluxtot_b6: burgers vector 1/2 [ 0  1 -1]
    // 4         ( 1  1  1), edge
    // 5         (-1  1  1), edge
    // 6         ( 1  1  1), screw
    // 7         (-1  1  1), screw
    slips.push_back(slip("[01-1](111)e", {0, 1, -1}, {1, 1, 1}));
    slips.push_back(slip("[01-1](-111)e", {0, 1, -1}, {-1, 1, 1}));
    slips.push_back(slip("[01-1](111)s", {0, 1, -1}, {1, 1, 1}));
    slips.push_back(slip("[01-1](-111)s", {0, 1, -1}, {-1, 1, 1}));
  }
};

template <class T>
inline void outProd(const vector<T>& a, const vector<T>& b,
                    vector<vector<T>>& M) {
  const int N = int(a.size());
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) M[i][j] = a[i] * b[j];
}

inline double derror(const double& a) {
  double r = fabs(a);
  return r < Mc ? r * r : Mc * (2 * r - Mc);
}

inline void split(const string& s, const char* delim, vector<string>& v) {
  // duplicate original string, return a char pointer and free  memories
  char* dup = strdup(s.c_str());
  char* token = strtok(dup, delim);
  while (token != NULL) {
    v.push_back(string(token));
    // the call is treated as a subsequent calls to strtok:
    // the function continues from where it left in previous invocation
    token = strtok(NULL, delim);
  }
  free(dup);
}

  // ------    -------------------------------------------
  // 1         Plastic strain
  // 2         Strain
  // 3         Flux due to climb
  // 8         Simulation timestep duration

#endif